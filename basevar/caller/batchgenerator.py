"""
This is a Process module for batchgenerator from BAM/CRAM

"""
from __future__ import division

import os
import sys

from pysam import FastaFile

from . import bam
from . import utils

REMOVE_BATCH_FILE = True


class BatchProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """

    def __init__(self, ref_file, align_files, regions, samples, mapq=10, batchcount=50,
                 out_batch_file=None, rerun=False):
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
        self.align_files = align_files
        self.out_batch_file = out_batch_file

        self.batch_count = batchcount
        self.samples = samples
        self.mapq = mapq
        self.smart_rerun = rerun

        # store the region into a dict
        self.regions = utils.regions2dict(regions)

    def run(self):
        """
        Run the process of calling variant and output files.

        """
        suffix_group = {}
        tmppd, name = os.path.split(os.path.realpath(self.out_batch_file))
        cache_dir = utils.safe_makedir(tmppd + "/Batchtemp.%s.WillBeDeletedWhenrJobsFinish" % name)

        if self.smart_rerun:
            # remove the last modification file
            utils.safe_remove(utils.get_last_modification_file(cache_dir))

        for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):

            # get fasta sequence of chrid
            fa = self.ref_file_hd.fetch(chrid)

            # create batch file for variant discovery
            region_batch_files = bam.create_batchfiles_for_regions(chrid,
                                                                   regions,
                                                                   self.batch_count,
                                                                   self.align_files,
                                                                   fa,
                                                                   self.mapq,
                                                                   cache_dir,
                                                                   sample_ids=self.samples,
                                                                   is_smart_rerun=self.smart_rerun)

            # The same group of samples has the same suffix in filename
            for filename in region_batch_files:
                # BaseVar.chrxx.xxx.xxx.2_10.batch.gz
                suffix = filename.split(".")[-3]  # get "2_10"
                if suffix not in suffix_group:
                    suffix_group[suffix] = []

                suffix_group[suffix].append(filename)

        # close reference file handle
        self.ref_file_hd.close()

        regions_batch_files = []
        for k, batch_files in suffix_group.items():
            # Merge different regions but the same samples
            region_batch_file = "%s/temp.all.%s.%s" % (cache_dir, k, name)
            regions_batch_files.append(region_batch_file)
            utils.merge_files(batch_files, region_batch_file, is_del_raw_file=True)

        # Merge the same region but different samples
        utils.merge_batch_files(regions_batch_files, self.out_batch_file, is_del_raw_file=True)

        if REMOVE_BATCH_FILE:
            try:
                os.removedirs(cache_dir)
            except OSError:
                sys.stderr.write("[WARNING] Directory not empty: %s, please delete it by yourself\n" % cache_dir)

        return
