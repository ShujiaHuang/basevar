"""
This is a Process module for BaseType by BAM/CRAM

"""
import os
import sys
import time
import multiprocessing

from pysam import FastaFile

from . import utils
from . import bam
from .basetypeprocess import variantsdiscovery

REMOVE_BATCH_FILE = False


class BaseVarSingleProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """

    def __init__(self, ref_file, aligne_files, in_popgroup_file, regions, samples, mapq=10, batchcount=1000,
                 out_vcf_file=None, out_cvg_file=None, rerun=False, cmm=None):
        """
        Constructor.

        Store input file, options and output file name.

        Parameters:
        ===========
            samples: list like
                A list of sample id

            regions: 2d-array like, required
                    It's region info , format like: [[chrid, start, end], ...]
        """
        self.ref_file_path = ref_file
        self.ref_file_hd = FastaFile(ref_file)
        self.aligne_files = aligne_files
        self.out_vcf_file = out_vcf_file
        self.out_cvg_file = out_cvg_file
        self.batch_count = batchcount
        self.samples = samples
        self.smart_rerun = rerun
        self.mapq = mapq
        self.cmm = cmm
        self.regions = {}

        # store the region into a dict
        for chrid, start, end in regions:

            if chrid not in self.regions:
                self.regions[chrid] = []

            self.regions[chrid].append([start, end])

        # loading population group
        # group_id => [a list samples_index]
        self.popgroup = {}
        if in_popgroup_file and len(in_popgroup_file):
            self.popgroup = utils.load_popgroup_info(self.samples, in_popgroup_file)

    def run(self):
        """
        Run the process of calling variant and output files.

        """
        vcf_header = utils.vcf_header_define(self.ref_file_path)
        group = []  # Just for the header of CVG file
        if self.popgroup:
            for g in self.popgroup.keys():
                g_id = g.split('_')[0]  # ignore '_AF'
                group.append(g_id)
                vcf_header.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="Allele frequency in the %s '
                                  'populations calculated base on LRT, in the range (0,1)">' % (g_id, g_id))

        CVG = open(self.out_cvg_file, 'w')
        CVG.write('##fileformat=CVGv1.0\n')
        CVG.write('##Group information is the depth of A:C:G:T:Indel\n')
        CVG.write('\t'.join(['#CHROM', 'POS', 'REF', 'Depth'] + self.cmm.BASE +
                            ['Indel', 'FS', 'SOR',
                             'Strand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)\t%s\n' % '\t'.join(group)]))

        VCF = open(self.out_vcf_file, 'w') if self.out_vcf_file else None
        if VCF:  # set header if VCF is not None
            VCF.write('\n'.join(vcf_header) + '\n')
            VCF.write('\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
                                 'INFO\tFORMAT'] + self.samples) + '\n')

        is_empty = True
        tmpd, name = os.path.split(os.path.realpath(self.out_cvg_file))
        batchfile_dir = utils.safe_makedir(tmpd + "/batchfiles.%s" % name.split(".cvg")[0])
        for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):

            # create batch file for variant discovery
            batchfiles = bam.create_batchfiles_for_regions(chrid,
                                                           regions,
                                                           self.batch_count,
                                                           self.aligne_files,
                                                           self.ref_file_hd.fetch(chrid),
                                                           self.mapq,
                                                           batchfile_dir,
                                                           sample_ids=self.samples,
                                                           is_smart_rerun=self.smart_rerun)

            # Process of variants discovery
            fa = self.ref_file_hd.fetch(chrid)  # get sequence of chrid from reference fasta
            _is_empty = variantsdiscovery(chrid, batchfiles, fa, self.popgroup, self.cmm, CVG, VCF)

            if not _is_empty:
                is_empty = False

            if REMOVE_BATCH_FILE:
                for f in batchfiles:
                    os.remove(f)

        CVG.close()
        if VCF:
            VCF.close()

        self.ref_file_hd.close()

        if is_empty:
            sys.stderr.write("\n************************************************************************\n"
                             "[WARNING] No reads are satisfy with the mapping quality (>=%d) in all your "
                             "bamfiles.\nWe get nothing in %s " % (self.mapq, self.out_cvg_file))
            if VCF:
                sys.stderr.write("and %s " % self.out_vcf_file)

            sys.stderr.write("Now the time is %s\n" % time.asctime())

        return


###############################################################################
class BaseVarMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.

    This class is much benefit than using ``BaseVarSingleProcess`` as a ``multiprocessing.Process`` directly
    It's a shield for ``BaseVarSingleProcess``
    """

    def __init__(self, ref_in_file, aligne_files, pop_group_file, regions, samples_id, mapq=10, batchcount=1000,
                 out_vcf_file=None, out_cvg_file=None, rerun=False, cmm=None):
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
                                                   regions,
                                                   samples_id,
                                                   mapq=mapq,
                                                   batchcount=batchcount,
                                                   out_cvg_file=out_cvg_file,
                                                   out_vcf_file=out_vcf_file,
                                                   rerun=rerun,
                                                   cmm=cmm)

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()
