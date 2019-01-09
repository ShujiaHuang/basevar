"""
This is a Process module for BaseType by BAM/CRAM

"""
import sys
import time
import multiprocessing

from pysam import FastaFile

from . import utils
from .basetypeprocess import batchfile_variants_discovery


class BaseVarBatchProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """

    def __init__(self, ref_file, batch_files, in_popgroup_file, regions, samples,
                 out_vcf_file=None, out_cvg_file=None, cmm=None):
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
        self.fa_file_hd = FastaFile(ref_file)
        self.batch_files = batch_files
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

        # loading population group
        # group_id => [a list samples_index]
        self.popgroup = {}
        if in_popgroup_file and len(in_popgroup_file):
            self.popgroup = utils.load_popgroup_info(self.samples, in_popgroup_file)

    def run(self):
        """
        Run the process of calling variant and output files.

        """
        info, group = [], []
        if self.popgroup:
            for g in self.popgroup.keys():
                g_id = g.split('_AF')[0]  # ignore '_AF'
                group.append(g_id)
                info.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="Allele frequency in the %s '
                               'populations calculated base on LRT, in the range (0,1)">' % (g_id, g_id))

        VCF = open(self.out_vcf_file, 'w') if self.out_vcf_file else None
        if VCF:
            vcf_header = utils.vcf_header_define(self.fa_file_hd.filename, info="\n".join(info), samples=self.samples)
            VCF.write("%s\n" % "\n".join(vcf_header))

        CVG = open(self.out_cvg_file, 'w')
        CVG.write('%s\n' % "\n".join(utils.cvg_header_define(group)))

        is_empty = True
        for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):

            # Process of variants discovery
            fa = self.fa_file_hd.fetch(chrid)
            _is_empty = batchfile_variants_discovery(chrid, regions, fa, self.batch_files, self.popgroup,
                                                     self.cmm, CVG, VCF)

            if not _is_empty:
                is_empty = False

        CVG.close()
        if VCF:
            VCF.close()

        self.fa_file_hd.close()

        if is_empty:
            sys.stderr.write("\n************************************************************************\n"
                             "[WARNING] No reads are satisfy with the mapping quality in all your\n"
                             "input files.\n We get nothing in %s " % (self.out_cvg_file))
            if VCF:
                sys.stderr.write("and %s " % self.out_vcf_file)

        sys.stderr.write("%s\n" % time.asctime())
        return


###############################################################################
class BaseVarBatchMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.

    This class is much benefit than using ``BaseVarSingleProcess`` as a ``multiprocessing.Process`` directly
    It's a shield for ``BaseVarSingleProcess``
    """

    def __init__(self, ref_in_file, align_files, pop_group_file, regions, samples_id,
                 out_vcf_file=None, out_cvg_file=None, cmm=None):
        """
        Constructor.

        regions: 2d-array like, required
                It's region info , format like: [[chrid, start, end], ... ]
        """
        multiprocessing.Process.__init__(self)

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.single_process = BaseVarBatchProcess(ref_in_file,
                                                  align_files,
                                                  pop_group_file,
                                                  regions,
                                                  samples_id,
                                                  out_cvg_file=out_cvg_file,
                                                  out_vcf_file=out_vcf_file,
                                                  cmm=cmm)

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()
