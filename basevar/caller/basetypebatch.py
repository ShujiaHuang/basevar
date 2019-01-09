"""
This is a Process module for batch BaseType by BAM/CRAM

"""
import sys

from pysam import FastaFile

from . import utils
from .basetypeprocess import output_header, batchfile_variants_discovery


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

        # store the region into a dict
        self.regions = utils.regions2dict(regions)

        # loading population group
        # group_id => [a list samples_index]
        self.popgroup = {}
        if in_popgroup_file and len(in_popgroup_file):
            self.popgroup = utils.load_popgroup_info(self.samples, in_popgroup_file)

    def run(self):
        """
        Run the process of calling variant and output files.

        """
        VCF = open(self.out_vcf_file, 'w') if self.out_vcf_file else None
        CVG = open(self.out_cvg_file, 'w')
        output_header(self.fa_file_hd.filename, self.samples, self.popgroup, CVG, out_vcf_handle=VCF)

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
                             "[WARNING] No reads are satisfy with the mapping quality in all your input\n\n"
                             "We get nothing in %s " % self.out_cvg_file)
            if VCF:
                sys.stderr.write("and %s " % self.out_vcf_file)

        return
