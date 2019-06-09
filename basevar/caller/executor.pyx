"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``BaseVar.py`` is "Sauron", and 'executor.py' module could just be called by it.
"""
from __future__ import division

import sys
import time

from pysam import tabix_index

from basevar.log import logger
from basevar import utils
from basevar.io.bam cimport get_sample_names
from . import CallerProcess, process_runner

from .basetypebam import BaseVarProcess
# from .coverageprocess import CvgSingleProcess


def _generate_regions_for_each_process(regions, process_num=1):
    """create regions for each process"""

    regions_for_each_process = [[] for _ in range(process_num)]

    # calculate region size for each processes
    total_regions_size = sum([e - s + 1 for _, s, e in regions])
    delta = int(total_regions_size / process_num) + 1 if total_regions_size % process_num else \
        int(total_regions_size / process_num)

    # reset
    total_regions_size = 0
    for chrid, start, end in regions:

        sub_region_size = end - start + 1
        # find index of regions_for_each_process by the total_regions_size
        i = int(total_regions_size / delta)
        total_regions_size += sub_region_size

        pre_size = 0
        if len(regions_for_each_process[i]):
            pre_size = sum([e - s + 1 for _, s, e in regions_for_each_process[i]])

        if sub_region_size + pre_size > delta:
            d = delta - pre_size
            regions_for_each_process[i].append([chrid, start, start + d - 1])

            i += 1
            for k, pos in enumerate(range(start + d - 1, end, delta)):
                s = pos + 1 if pos + 1 < end else end
                e = pos + delta if pos + delta < end else end

                regions_for_each_process[i + k].append([chrid, s, e])

        else:
            regions_for_each_process[i].append([chrid, start, end])

    return regions_for_each_process


def _output_cvg_and_vcf(sub_cvg_files, sub_vcf_files, outcvg, outvcf=None):
    """CVG file and VCF file could use the same tabix strategy."""
    for out_final_file, sub_file_list in zip([outcvg, outvcf], [sub_cvg_files, sub_vcf_files]):

        if out_final_file:
            _output_file(sub_file_list, out_final_file, del_raw_file=True)

    return


def _output_file(sub_files, out_file_name, del_raw_file=False):
    if out_file_name.endswith(".gz"):
        utils.merge_files(sub_files, out_file_name, output_isbgz=True, is_del_raw_file=del_raw_file)

        # Column indices are 0-based. Note: this is different from the tabix command line
        # utility where column indices start at 1.
        tabix_index(out_file_name, force=True, seq_col=0, start_col=1, end_col=1)
    else:
        utils.merge_files(sub_files, out_file_name, is_del_raw_file=del_raw_file)

    return


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
        self.outcvg = args.outcvg if args.outcvg else None
        self.options = args

        # setting the resolution of MAF
        self.options.min_af = utils.set_minaf(len(self.alignfiles)) if (args.min_af is None) else args.min_af
        logger.info("Finish loading arguments and we have %d BAM/CRAM files for "
                    "variants calling." % len(self.alignfiles))

        # Loading positions or load all the genome regions
        regions = utils.load_target_position(self.reference_file, args.positions, args.regions)
        self.regions_for_each_process = _generate_regions_for_each_process(regions, process_num=self.nCPU)

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

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s, %s]\n' % (i + 1, self.nCPU, sub_vcf_file, sub_cvg_file))

            processes.append(CallerProcess(BaseVarProcess,
                                           self.sample_id,
                                           self.alignfiles,
                                           self.reference_file,
                                           self.regions_for_each_process[i],
                                           out_cvg_file=sub_cvg_file,
                                           out_vcf_file=sub_vcf_file,
                                           options=self.options))

        process_runner(processes)

        # Final output
        _output_cvg_and_vcf(out_cvg_names, out_vcf_names, self.outcvg, outvcf=self.outvcf)

        return processes


class CoverageRunner(object):

    def __init__(self, args, cmm=utils.CommonParameter()):
        """init function
        """
        self.referencefile = args.referencefile
        self.nCPU = args.nCPU
        self.outputfile = args.outputfile
        self.cmm = cmm

        # Loading positions if not provid we'll load all the genome
        self.regions = utils.load_target_position(self.referencefile, args.positions, args.regions)

        # Get all the input alignement files
        self.alignefiles = utils.load_file_list(args.infilelist)
        sys.stderr.write('[INFO] Finish loading parameters and input file list %s\n' % time.asctime())

    def run(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call varaintis by BaseType ... %s\n' % time.asctime())

        # Always create process manager even if nCPU==1, so that we can
        # listen for signals from main thread
        regions_for_each_process = [[] for _ in range(self.nCPU)]
        if len(self.regions) < self.nCPU:
            # We cut the region evenly to fit nCPU if regions < nCPU
            for chrid, start, end in self.regions:
                delta = int((end - start + 1) / self.nCPU)
                if delta == 0:
                    delta = 1

                for i, pos in enumerate(range(start - 1, end, delta)):
                    s = pos + 1 if pos + 1 < end else end
                    e = pos + delta if pos + delta < end else end

                    regions_for_each_process[i % self.nCPU].append([chrid, s, e])

        else:
            for i, region in enumerate(self.regions):
                regions_for_each_process[i % self.nCPU].append(region)

        out_cvg_names = set()
        processes = []
        for i in range(self.nCPU):
            sub_cvg_file = self.outputfile + '_temp_%s' % i

            out_cvg_names.add(sub_cvg_file)
            # processes.append(CallerProcess(CvgSingleProcess,
            #                                self.referencefile,
            #                                self.alignefiles,
            #                                sub_cvg_file,
            #                                regions_for_each_process[i],
            #                                cmm=self.cmm))

        process_runner(processes)

        # Final output file name
        utils.merge_files(out_cvg_names, self.outputfile, is_del_raw_file=True)

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
        _output_file(self.inputfiles, self.outputfile)
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

        return self
