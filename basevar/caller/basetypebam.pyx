"""
This is a Process module for BaseType by BAM/CRAM

"""
from __future__ import division

import os
import sys

from basevar.io.fasta import FastaFile

from basevar import utils
from basevar.caller.batchcaller cimport create_batchfiles_in_regions
from .basetypeprocess import output_header, variants_discovery

cdef bint REMOVE_BATCH_FILE = 0


class BaseVarProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """

    def __init__(self, samples, align_files, ref_file, regions, out_vcf_file=None,
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

        self.smart_rerun = True if options.smartrerun else False
        self.options = options

        # store the region into a dict
        self.regions = utils.regions2dict(regions)

        # loading population group
        # group_id => [a list samples_index]
        self.popgroup = {}
        if options.pop_group_file and len(options.pop_group_file):
            self.popgroup = utils.load_popgroup_info(self.samples, options.pop_group_file)

    def run(self):
        """Run the process of calling variant and output files.
        """
        VCF = open(self.out_vcf_file, "w") if self.out_vcf_file else None
        CVG = open(self.out_cvg_file, "w")
        output_header(self.fa_file_hd.filename, self.samples, self.popgroup, CVG, out_vcf_handle=VCF)

        is_empty = True
        tmpd, name = os.path.split(os.path.realpath(self.out_cvg_file))
        cache_dir = utils.safe_makedir(tmpd + "/Batchfiles.%s.WillBeDeletedWhenJobsFinish" % name)

        if self.smart_rerun:
            # remove the last modification file
            utils.safe_remove(utils.get_last_modification_file(cache_dir))

        for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):
            # refseq = self.fa_file_hd.fetch(chrid)
            # create batch file for variant discovery
            batchfiles = create_batchfiles_in_regions(chrid,
                                                      regions,
                                                      self.align_files,
                                                      self.fa_file_hd,
                                                      self.samples,
                                                      cache_dir,
                                                      self.options,
                                                      self.smart_rerun)

            # Process of variants discovery
            # _is_empty = variants_discovery(chrid, refseq, batchfiles, self.popgroup, self.options.min_af, CVG, VCF)
            #
            # if not _is_empty:
            #     is_empty = False

            if REMOVE_BATCH_FILE:
                for f in batchfiles:
                    os.remove(f)

        CVG.close()
        if VCF:
            VCF.close()

        self.fa_file_hd.close()

        if is_empty:
            sys.stderr.write("\n***************************************************************************\n"
                             "[WARNING] No reads are satisfy with the mapping quality (>=%d) in all of your\n"
                             "input files. We get nothing in %s " % (self.options.mapq, self.out_cvg_file))
            if VCF:
                sys.stderr.write("and %s " % self.out_vcf_file)

            sys.stderr.write("\n\n")

        if REMOVE_BATCH_FILE:
            try:
                os.removedirs(cache_dir)
            except OSError:
                sys.stderr.write("[WARNING] Directory not empty: %s, please delete it by yourself\n" % cache_dir)

        return
