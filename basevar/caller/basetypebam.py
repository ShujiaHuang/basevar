"""
This is a Process module for BaseType by BAM/CRAM

"""
from __future__ import division

import os
import sys

from pysam import FastaFile

from . import bam
from . import utils
from .basetypeprocess import output_header, variants_discovery

REMOVE_BATCH_FILE = False


class BaseVarProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """

    def __init__(self, ref_file, align_files, in_popgroup_file, regions, samples, mapq=10, batchcount=50,
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
        self.fa_file_hd = FastaFile(ref_file)
        self.align_files = align_files
        self.out_vcf_file = out_vcf_file
        self.out_cvg_file = out_cvg_file

        self.batch_count = batchcount
        self.samples = samples
        self.smart_rerun = rerun
        self.mapq = mapq
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
        tmpd, name = os.path.split(os.path.realpath(self.out_cvg_file))
        cache_dir = utils.safe_makedir(tmpd + "/Batchfiles.%s.WillBeDeletedWhenJobsFinish" % name)

        if self.smart_rerun:
            # remove the last modification file
            utils.safe_remove(utils.get_last_modification_file(cache_dir))

        for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):

            # get fasta sequence by chrid
            fa = self.fa_file_hd.fetch(chrid)

            # create batch file for variant discovery
            batchfiles = bam.create_batchfiles_for_regions(chrid,
                                                           regions,
                                                           self.batch_count,
                                                           self.align_files,
                                                           fa,
                                                           self.mapq,
                                                           cache_dir,
                                                           sample_ids=self.samples,
                                                           is_smart_rerun=self.smart_rerun)

            # Process of variants discovery
            _is_empty = variants_discovery(chrid, fa, batchfiles, self.popgroup, self.cmm, CVG, VCF)

            if not _is_empty:
                is_empty = False

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
                             "input files.\n We get nothing in %s " % (self.mapq, self.out_cvg_file))
            if VCF:
                sys.stderr.write("and %s " % self.out_vcf_file)

        if REMOVE_BATCH_FILE:
            try:
                os.removedirs(cache_dir)
            except OSError:
                sys.stderr.write("[WARNING] Directory not empty: %s, please delete it by yourself\n" % cache_dir)

        return
