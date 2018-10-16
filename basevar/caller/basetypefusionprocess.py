"""
This is a process module of BaseType for fusion file
"""
import sys
import time
import multiprocessing

import pysam

from . import utils
from .basetypeprocess import basetypeprocess


class BaseTypeFusionSingleProcess(object):
    """
    simple class to repesent a single process.
    """
    def __init__(self, in_ref_file, in_fusion_files, in_popgroup_file,
                 out_vcf_file, out_cvg_file, regions, samples, cmm=None):
        """
        Store input file, options and output file name.

        Parameters:
            ``in_fusion_files``: Array like
                fusion file is 0-base coordinate system
        ===========

            samples: list like
                A list of sample id
        """
        self.ref_file_hd = pysam.FastaFile(in_ref_file)
        self.in_fusion_files = in_fusion_files
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

        # Cache a batch of file handle which index by tabix
        self.tb_files = [pysam.TabixFile(f) for f in self.in_fusion_files]

        # loading population group
        # group_id => [a list samples_index]
        self.popgroup = {}
        if in_popgroup_file and len(in_popgroup_file):
            self.popgroup = utils.load_popgroup_info(self.samples, in_popgroup_file)

    def _close_file(self):

        # close the reference file
        self.ref_file_hd.close()

        # close all of the fusion files
        for tb in self.tb_files:
            tb.close()

        return

    def run(self):
        """
        Run the process of calling variant and output
        """
        vcf_header = utils.vcf_header_define()
        group = [] # Just for the header of CVG file
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
                # ``regions`` is a 2-D array : [[start1,end1], [start2,end2], ...]
                # ``iter_tokes`` is a list of iterator for each sample's input file
                tmp_region = []
                for p in regions:  # covert to 1d-array
                    tmp_region.extend(p)

                tmp_region = sorted(tmp_region)
                start, end = tmp_region[0], tmp_region[-1]
                iter_tokes = []

                # All the fusion file is 0-base system
                for i, bt in enumerate(self.tb_files):
                    try:
                        iter_tokes.append(bt.fetch(chrid, start-1, end))
                    except ValueError:
                        if self.cmm.debug:
                            sys.stderr.write("# [WARMING] Empty region %s:%d-%d in %s" %
                                             (chrid, start-1, end, self.in_fusion_files[i]))
                        iter_tokes.append('')

                # get sequence of chrid from reference fasta
                fa = self.ref_file_hd.fetch(chrid)

                # Init the sample information
                n = 0
                sample_info = [utils.fetch_next(it) for it in iter_tokes]
                for start, end in regions:

                    sys.stderr.write('[INFO] Fetching info from %d samples in region %s'
                                     ' at %s\n' % (len(iter_tokes),
                                                   chrid + ":" + str(start) + "-" + str(end),
                                                   time.asctime()))

                    for position in xrange(start, end + 1):

                        if n % 100000 == 0:
                            sys.stderr.write("[INFO] loading lines %d at position %s:%d\t%s\n" %
                                             (n+1, chrid, position, time.asctime()))

                        n += 1
                        # sample_base, sample_base_qual, strands, mapqs and
                        # read_pos_rank are listed the same orde with each other.
                        (sample_bases, sample_base_quals, strands, mapqs, read_pos_rank,
                         indels) = self.fetch_base_by_position(
                            fa,  # Fa sequence for indel sequence
                            position - 1,  # postion for fusion file is 0-base
                            sample_info,
                            iter_tokes
                        )

                        ref_base = fa[position-1]

                        # ignore positions if coverage=0 or ref base is 'N' base
                        if (not sample_bases) or (ref_base.upper() not in ['A', 'C', 'G', 'T']):
                            continue

                        basetypeprocess(chrid, position, ref_base, sample_bases, sample_base_quals,
                                        mapqs, strands, indels, read_pos_rank, self.popgroup,
                                        self.cmm, CVG, VCF)

                        # These two lines just for debug.
                        # VCF.write("\t".join([chrid, str(position), ".", ref_base, ".\t.\t.\t.\t."]+[":".join(map(str, [a,b,c])) for a,b,c in zip(sample_bases,sample_base_quals,strands)]))
                        # sys.exit(1)

        self._close_file()

    def fetch_base_by_position(self, fa, position, sample_info, iter_tokes):
        """
        sample_info example: 0-base
            '#CHROM  START  END  TYPE  MAPQ  SO  Read_POS   BASE_QUAL'
            'chr3	1002728	1002776	<NON_REF>	60	+	16	G'
        """

        base_quals = []
        bases = []
        strands = []
        mapqs = []
        read_pos_rank = []
        indels = []

        for i, sample_pos_line in enumerate(sample_info):
            # 'chr3	1002728	1002776	<NON_REF>	60	+	16	G'

            bs, qs, strand, mapq, rpr, indel, sample_info[i] = (
                self.seek_position(fa, position, sample_pos_line, iter_tokes[i])
            )

            bases.append(bs)
            base_quals.append(qs)
            strands.append(strand)
            mapqs.append(mapq)
            read_pos_rank.append(rpr)

            if indel:
                indels.append(indel.upper())
            else:
                indels.append("")

        return bases, base_quals, strands, mapqs, read_pos_rank, indels

    def seek_position(self, fa, target_pos, sample_line, sample_iter):
        """
        Get mapping info for specific position.
            `fa`: Just for scan indel
        """
        base, strand, indel, qual, rpr, mapq = 'N', '.', '', '!', 0, 0  # Init
        if sample_line:
            # 'chr3	1002728	1002776	<NON_REF>	60	+	16	G'
            tok = sample_line.strip().split()
            start = int(tok[1])
            end = int(tok[2])

            # target_pos must be in the [``start``, ``end``), and not allow == ``end``
            # cause this it's the 0-base coordinate system
            if target_pos >= end:

                while target_pos >= end:

                    sample_line = utils.fetch_next(sample_iter)
                    if sample_line:

                        tok = sample_line.strip().split()
                        start = int(tok[1])
                        end = int(tok[2])
                    else:
                        # The end of file. Break the loop.
                        break

            # In case of sample_line may hit the end of file
            if sample_line and target_pos >= start:

                mapq = int(tok[4])
                strand = tok[5]
                rpr = target_pos - start + int(tok[6])
                qual = tok[7]

                base = tok[3]
                if base == '<NON_REF>':
                    # reference base.
                    # target_pos is 0-base
                    base = fa[target_pos].upper()
                elif base == '.':
                    # deletion
                    indel = '-' + fa[start:end]

                elif len(base) > 1:
                    # insertion
                    indel = '+' + base

        qual = ord(qual) - 33
        return base, qual, strand, mapq, rpr, indel, sample_line


###############################################################################
class BaseVarFusionMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.
    """
    def __init__(self, ref_in_file, fusion_files, pop_group_file,
                 out_vcf_file, out_cvg_file, regions,samples, cmm=None):
        """
        Constructor.

        regions: 2d-array like, required
                It's region info , format like: [[chrid, start, end], ...]
        """
        multiprocessing.Process.__init__(self)

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.single_process = BaseTypeFusionSingleProcess(ref_in_file,
                                                          fusion_files,
                                                          pop_group_file,
                                                          out_vcf_file,
                                                          out_cvg_file,
                                                          regions,
                                                          samples,
                                                          cmm=cmm)

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()

