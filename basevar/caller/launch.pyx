# cython: profile=True
"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``runner.py`` is "Sauron", and 'launch.pyx' module could just be called by it.
"""
from __future__ import division

import os
import sys
import time

from basevar.log import logger
from basevar import utils
from basevar.utils cimport generate_regions_by_process_num

from basevar.io.bam cimport get_sample_names
from basevar.caller.do import CallerProcess, process_runner
from basevar.caller.basetypeprocess cimport BaseVarProcess


class BaseTypeRunner(object):

    def __init__(self, args):
        """init function
        """
        # setting parameters
        self.alignfiles = args.input
        if args.infilelist:
            self.alignfiles += utils.load_file_list(args.infilelist)

        self.nCPU = args.nCPU
        self.reference_file = args.referencefile
        self.outvcf = args.outvcf if args.outvcf else None
        self.outcvg = args.outcvg
        self.options = args

        # setting the resolution of MAF
        self.options.min_af = utils.set_minaf(len(self.alignfiles)) if (args.min_af is None) else args.min_af
        logger.info("Finish loading arguments and we have %d BAM/CRAM files for "
                    "variants calling." % len(self.alignfiles))

        # Loading positions if not been provided we'll load all the genome
        regions = utils.load_target_position(self.reference_file, args.positions, args.regions)
        self.regions_for_each_process = generate_regions_by_process_num(regions, process_num=self.nCPU,
                                                                        convert_to_2d=False)

        # ``samples_id`` has the same size and order as ``aligne_files``
        self.sample_id = get_sample_names(self.alignfiles, True if args.filename_has_samplename else False)

    def basevar_caller(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call variants by BaseType ... %s\n' % time.asctime())

        out_vcf_names = []
        out_cvg_names = []

        processes = []
        # Always create process manager even if nCPU==1, so that we can
        # listen signals from main thread
        for i in range(self.nCPU):
            sub_cvg_file = self.outcvg + '_temp_%s' % i
            out_cvg_names.append(sub_cvg_file)

            if self.outvcf:
                sub_vcf_file = self.outvcf + '_temp_%s' % i
                out_vcf_names.append(sub_vcf_file)
            else:
                sub_vcf_file = None

            tmpd, name = os.path.split(os.path.realpath(sub_cvg_file))
            cache_dir = tmpd + "/Batchfiles.%s.WillBeDeletedWhenJobsFinish" % name

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s, %s]\n' % (i + 1, self.nCPU, sub_vcf_file, sub_cvg_file))

            if self.options.smartrerun and os.path.isfile(sub_cvg_file) and (not os.path.exists(cache_dir)):
                # if `cache_dir` is not exist and `sub_cvg_file` is exists means
                # `sub_cvg_file and sub_vcf_file` has been finish successfully.
                continue

            cache_dir = utils.safe_makedir(cache_dir)
            processes.append(CallerProcess(BaseVarProcess,
                                           self.sample_id,
                                           self.alignfiles,
                                           self.reference_file,
                                           self.regions_for_each_process[i],
                                           out_cvg_file=sub_cvg_file,
                                           out_vcf_file=sub_vcf_file,
                                           cache_dir=cache_dir,
                                           options=self.options))

        process_runner(processes)

        # Final output
        utils.output_cvg_and_vcf(out_cvg_names, out_vcf_names, self.outcvg, outvcf=self.outvcf)
        return processes

    def basevar_caller_singleprocess(self):
        """
        Run variant caller --------- Just for Testting, when we done, please delete this function!!!!!!
        """
        sys.stderr.write('[INFO] Start call variants by BaseType ... %s\n' % time.asctime())

        out_vcf_names = []
        out_cvg_names = []

        processes = []
        # Always create process manager even if nCPU==1, so that we can
        # listen signals from main thread
        for i in range(self.nCPU):
            sub_cvg_file = self.outcvg + '_temp_%s' % i
            out_cvg_names.append(sub_cvg_file)

            if self.outvcf:
                sub_vcf_file = self.outvcf + '_temp_%s' % i
                out_vcf_names.append(sub_vcf_file)
            else:
                sub_vcf_file = None

            tmpd, name = os.path.split(os.path.realpath(sub_cvg_file))
            cache_dir = tmpd + "/Batchfiles.%s.WillBeDeletedWhenJobsFinish" % name

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s, %s]\n' % (i + 1, self.nCPU, sub_vcf_file, sub_cvg_file))

            if self.options.smartrerun and os.path.isfile(sub_cvg_file) and (not os.path.exists(cache_dir)):
                # if `cache_dir` is not exist and `sub_cvg_file` is exists means
                # `sub_cvg_file and sub_vcf_file` has been finish successfully.
                continue

            cache_dir = utils.safe_makedir(cache_dir)
            bp = BaseVarProcess(self.sample_id,
                                self.alignfiles,
                                self.reference_file,
                                self.regions_for_each_process[i],
                                out_cvg_file=sub_cvg_file,
                                out_vcf_file=sub_vcf_file,
                                cache_dir=cache_dir,
                                options=self.options)
            bp.run()

        # Final output
        utils.output_cvg_and_vcf(out_cvg_names, out_vcf_names, self.outcvg, outvcf=self.outvcf)
        return processes


class MergeRunner(object):
    """Runner for merging files"""

    def __init__(self, args):
        """init function"""

        self.inputfiles = args.input
        if args.infilelist:
            self.inputfiles += utils.load_file_list(args.infilelist)

        self.outputfile = args.outputfile

    def run(self):
        utils.output_file(self.inputfiles, self.outputfile)
        return


class NearbyIndelRunner(object):
    """Add Nearby Indel density and type information for each variants of VCF"""

    def __init__(self, args):
        """init function"""
        args.nearby_dis_around_indel = int(args.nearby_dis_around_indel)
        self.in_vcf_file = args.in_vcf_file
        self.in_cvg_file = args.in_cvg_file
        self.nearby_dis_around_indel = args.nearby_dis_around_indel

        sys.stderr.write('[INFO] basevar NearbyIndel'
                         '\n\t-i %s'
                         '\n\t-c %s'
                         '\n\t-d %d' % (args.in_vcf_file,
                                        args.in_cvg_file,
                                        args.nearby_dis_around_indel)
                         )

    def run(self):
        from .other import NearbyIndel

        nbi = NearbyIndel(self.in_vcf_file, self.in_cvg_file, nearby_distance=self.nearby_dis_around_indel)
        nbi.run()

        return
