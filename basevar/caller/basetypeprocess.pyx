# cython: profile=True
"""
This is a Process module for BaseType by BAM/CRAM

"""
import sys
import time

from basevar.io.fasta import FastaFile

from basevar.log import logger
from basevar import utils

from basevar.caller.variantcaller cimport variant_discovery_in_regions


cdef class BaseVarProcess:
    """
    simple class to repesent a single BaseVar process.
    """
    def __cinit__(self, samples, align_files, ref_file, regions, out_vcf_file=None,
                  out_cvg_file=None, options=None):
        """Constructor.

        Store input file, options and output file name.

        Parameters:
        ===========
            samples: list like
                A list of sample id

            regions: 2d-array like, required
                    It's region info , format like: [[chrid, start, end], ...]
        """
        self.samples = samples
        self.align_files = align_files
        self.fa_file_hd = FastaFile(ref_file, ref_file + ".fai")
        self.out_vcf_file = out_vcf_file
        self.out_cvg_file = out_cvg_file
        self.options = options

        self.regions = regions
        self.dict_regions = utils.regions2dict(regions)

        # loading population group
        # group_id => [a list samples_index]
        self.popgroup = {}
        if options.pop_group_file and len(options.pop_group_file):
            self.popgroup = utils.load_popgroup_info(self.samples, options.pop_group_file)

    def run(self):
        self.run_variant_discovery_in_regions()  # do not create batch files
        return

    cdef void run_variant_discovery_in_regions(self):
        """Run the process of calling variant without creating batch files.
        This function will hit A BIG IO problem when we need to read huge number of BAM files.
        """
        start_time = time.time()
        cdef bint is_empty
        try:
            is_empty = variant_discovery_in_regions(
                self.fa_file_hd,
                self.align_files,
                self.regions,
                self.samples,
                self.popgroup,
                self.out_cvg_file,
                self.out_vcf_file,
                self.options
            )

        except Exception, e:
            logger.error("Error happen in run_variant_discovery_in_regions(): %s" % e)
            sys.exit(1)

        logger.info("Running variant_discovery_in_regions done, %d seconds elapsed.\n" % (time.time() - start_time))
        if is_empty:
            logger.warning("\n***************************************************************************\n"
                           "[WARNING] No reads are satisfy with the mapping quality (>=%d) in all of your\n"
                           "input files. We get nothing in %s \n\n" % (self.options.mapq, self.out_cvg_file))
            if self.out_vcf_file:
                logger.warning("and %s " % self.out_vcf_file)

        return
