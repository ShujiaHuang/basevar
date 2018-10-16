"""
This is a Process module for BaseType by BAM/CRAM

"""
import sys
import time
import multiprocessing

import pysam

from . import utils
from . import bam
from .basetypeprocess import basetypeprocess


class BaseVarSingleProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """
    def __init__(self, ref_file, aligne_files, in_popgroup_file, out_vcf_file,
                 out_cvg_file, regions, samples, cmm=None):
        """
        Store input file, options and output file name.

        Parameters:
        ===========

            samples: list like
                A list of sample id
        """
        self.ref_file_hd = pysam.FastaFile(ref_file)
        self.aligne_files = aligne_files
        self.out_vcf_file = out_vcf_file
        self.out_cvg_file = out_cvg_file
        self.samples = samples
        self.cmm = cmm
        self.regions = {}

        # store the region into a dict
        for chrid, start, end in regions:

            if chrid not in self.regions:
                self.regions[chrid] = []

            self.regions[chrid].append([start, end])

        # Cache a batch of aligne file handle
        self.ali_files_hd = []
        for f in self.aligne_files:

            try:
                bf = pysam.AlignmentFile(f)

            except ValueError:
                sys.stderr.write('[ERROR] Input file: %s is not BAM nor CRAM.\n' % f)
                self._close_aligne_file()
                sys.exit(1)

            self.ali_files_hd.append(bf)

            # loading population group
            # group_id => [a list samples_index]
            self.popgroup = {}
            if in_popgroup_file and len(in_popgroup_file):
                self.popgroup = utils.load_popgroup_info(self.samples, in_popgroup_file)

    def _close_aligne_file(self):

        self.ref_file_hd.close()
        for bf in self.ali_files_hd:
            bf.close()

        return

    def run(self):
        """
        Run the process of calling variant and output
        """
        vcf_header = utils.vcf_header_define()
        group = []  # Just for the header of CVG file
        if self.popgroup:
            for g in self.popgroup.keys():
                g_id = g.split('_')[0]  # ignore '_AF'
                group.append(g_id)
                vcf_header.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="Allele '
                                  'frequency in the %s populations calculated base on LRT, in '
                                  'the range (0,1)">' % (g_id, g_id))

        with open(self.out_vcf_file, 'w') as VCF, open(self.out_cvg_file, 'w') as CVG:

            CVG.write('##fileformat=CVGv1.0\n')
            CVG.write('##Group_info is the depth of A:C:G:T:Indel\n')
            CVG.write('\t'.join(['#CHROM', 'POS', 'REF', 'Depth'] + self.cmm.BASE +
                                ['Indel', 'FS', 'SOR', 'Strand_Coverage(REF_FWD,'
                                 'REF_REV,ALT_FWD,ALT_REV)\t%s\n' % '\t'.join(group)]))

            # set header
            VCF.write('\n'.join(vcf_header) + '\n')
            VCF.write('\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
                                 'INFO\tFORMAT'] + self.samples) + '\n')

            for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):
                # ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
                # ``iter_tokes`` is a list of iterator for each sample's input file
                tmp_region = []
                for p in regions:  # covert to 1d-array
                    tmp_region.extend(p)

                tmp_region = sorted(tmp_region)

                start, end = tmp_region[0], tmp_region[-1]
                iter_tokes = []
                for i, bf in enumerate(self.ali_files_hd):
                    try:
                        # 0-base
                        iter_tokes.append(bf.pileup(chrid, start-1, end))
                    except ValueError:
                        if self.cmm.debug:
                            sys.stderr.write("# [WARMING] Empty region %s:%d-%d in %s" %
                                             (chrid, start-1, end, self.aligne_files[i]))
                        iter_tokes.append('')

                # get sequence of chrid from reference fasta
                fa = self.ref_file_hd.fetch(chrid)

                # Set iteration marker: 1->iterate; 0->Do not iterate or hit the end
                n = 0
                sample_info = [utils.fetch_next(it) for it in iter_tokes]
                for start, end in regions:

                    sys.stderr.write('[INFO] Fetching info from %d samples in region %s'
                                     ' at %s\n' % (len(iter_tokes),
                                                   chrid + ":" + str(start) + "-" + str(end),
                                                   time.asctime())
                                     )

                    for position in xrange(start, end + 1):

                        if n % 100000 == 0:
                            sys.stderr.write("[INFO] loading lines %d at position %s:%d\t%s\n" %
                                             (n+1, chrid, position, time.asctime()))
                        n += 1

                        (sample_bases, sample_base_quals, strands, mapqs, read_pos_rank,
                         indels) = bam.fetch_base_by_position(
                            position - 1,  # postion for pysam is 0-base
                            sample_info,
                            iter_tokes,
                            fa,  # Fa sequence for indel sequence
                            is_scan_indel=True
                        )

                        ref_base = fa[position-1]

                        # ignore positions if coverage=0 or ref base is 'N' base
                        if (not sample_bases) or (ref_base.upper() not in ['A', 'C', 'G', 'T']):
                            continue

                        # Calling varaints by Basetypes and output VCF and Coverage files.
                        basetypeprocess(chrid, position, ref_base, sample_bases, sample_base_quals,
                                        mapqs, strands, indels, read_pos_rank, self.popgroup,
                                        self.cmm, CVG, VCF)

        self._close_aligne_file()


###############################################################################
class BaseVarMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.
    """
    def __init__(self, ref_in_file, aligne_files, pop_group_file, out_vcf_file,
                 out_cvg_file, regions, samples_id, cmm=None):
        """
        Constructor.

        regions: 2d-array like, required
                It's region info , format like: [[chrid, start, end], ...]
        """
        multiprocessing.Process.__init__(self)

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.single_process = BaseVarSingleProcess(ref_in_file,
                                                   aligne_files,
                                                   pop_group_file,
                                                   out_vcf_file,
                                                   out_cvg_file,
                                                   regions,
                                                   samples_id,
                                                   cmm=cmm)

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()
