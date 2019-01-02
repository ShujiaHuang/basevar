"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``BaseVar.py`` is "Sauron", and 'executor.py' module could just be called by it.
"""
from __future__ import division

import os
import sys
import time

from pysam import AlignmentFile

from . import utils
from .basetypebam import BaseVarMultiProcess
from .coverageprocess import CvgMultiProcess


# from .vqsr import vqsr


class BaseTypeBamRunner(object):

    def __init__(self, args, cmm=utils.CommonParameter()):
        """init function
        """
        # setting parameters
        self.nCPU = args.nCPU
        self.mapq = args.mapq
        self.referencefile = args.referencefile
        self.pop_group_file = args.pop_group_file
        self.batchcount = args.batchcount

        self.outvcf = args.outvcf if args.outvcf else None
        self.outcvg = args.outcvg

        self.smartrerun = True if args.smartrerun else False
        if self.smartrerun:
            sys.stderr.write("************************************************\n"
                             "******************* WARNING ********************\n"
                             "************************************************\n"
                             ">>>>>>>> You have setted `smart rerun` <<<<<<<<<\n"
                             "Please make sure that all the parameters are the\n"
                             "same with your previous commands.\n"
                             "************************************************\n\n")

        # Loading positions or load all the genome regions
        self.regions = utils.load_target_position(self.referencefile, args.positions, args.regions)

        # Make sure you have input at least one bamfile.
        if not args.input and not args.infilelist:
            sys.stderr.write("[ERROR] Missing input BAM/CRAM files.\n\n")
            sys.exit(1)

        self.alignefiles = args.input
        if args.infilelist:
            self.alignefiles += utils.load_file_list(args.infilelist)

        # setting the resolution of MAF
        self.cmm = cmm
        if args.min_af is None:
            args.min_af = min(100.0 / len(self.alignefiles), 0.001, self.cmm.MINAF)

        self.cmm.MINAF = args.min_af
        sys.stderr.write('[INFO] Finish loading parameters and we have %d BAM/CRAM files '
                         'for variants calling at %s\n' % (len(self.alignefiles), time.asctime()))

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.sample_id = self._load_sample_id_from_bam(filename_has_samplename=args.filename_has_samplename)

    def _load_sample_id_from_bam(self, filename_has_samplename=True):
        """loading sample id in BAM/CRMA files from RG tag"""

        sys.stderr.write('[INFO] Start loading all samples\' id from alignment files\n')
        if filename_has_samplename:
            sys.stderr.write('[INFO] loading samples\' id from filename because you '
                             'set "--filename-has-samplename"\n')
        sample_id = []
        for i, al in enumerate(self.alignefiles):

            if i % 1000 == 0:
                sys.stderr.write("[INFO] loading %d/%d alignment files ... %s\n" %
                                 (i + 1, len(self.alignefiles), time.asctime()))

            if filename_has_samplename:
                filename = os.path.basename(al)

                # sample id should be the first element separate by ".",
                # e.g: "CL100045504_L02_61.sorted.rmdup.realign.BQSR.bam", "CL100045504_L02_61" is sample id.
                sample_id.append(filename.split(".")[0])

            else:

                # This situation will take a very long time to get sampleID from BAM header.
                bf = AlignmentFile(al)
                if 'RG' not in bf.header:
                    sys.stderr.write('[ERROR] Bam file format error: missing @RG in the header.\n')
                    bf.close()
                    sys.exit(1)

                sample_id.append(bf.header['RG'][0]['SM'])
                bf.close()

        sys.stderr.write('[INFO] Finish load all %d samples\' ID '
                         'from RG tag\n\n' % len(sample_id))
        return sample_id

    def run(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call variants by BaseType ... %s\n' % time.asctime())

        # Always create process manager even if nCPU==1, so that we can
        # listen signals from main thread
        regions_for_each_process = [[] for _ in range(self.nCPU)]
        if len(self.regions) < self.nCPU:
            # We cut the region into pieces to fit nCPU if regions < nCPU
            for chrid, start, end in self.regions:
                delta = int((end - start + 1) / self.nCPU)
                if delta == 0:
                    delta = 1

                for i, pos in enumerate(xrange(start - 1, end, delta)):
                    s = pos + 1 if pos + 1 < end else end
                    e = pos + delta if pos + delta < end else end

                    regions_for_each_process[i % self.nCPU].append([chrid, s, e])

        else:
            for i, region in enumerate(self.regions):
                regions_for_each_process[i % self.nCPU].append(region)

        out_vcf_names = set()
        out_cvg_names = set()

        processes = []
        for i in range(self.nCPU):
            sub_cvg_file = self.outcvg + '_temp_%s' % i
            out_cvg_names.add(sub_cvg_file)

            if self.outvcf:
                sub_vcf_file = self.outvcf + '_temp_%s' % i
                out_vcf_names.add(sub_vcf_file)
            else:
                sub_vcf_file = None

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s, %s]\n' % (i + 1, self.nCPU, sub_vcf_file,
                                             sub_cvg_file))

            processes.append(BaseVarMultiProcess(self.referencefile,
                                                 self.alignefiles,
                                                 self.pop_group_file,
                                                 regions_for_each_process[i],
                                                 self.sample_id,
                                                 mapq=self.mapq,
                                                 batchcount=self.batchcount,
                                                 out_cvg_file=sub_cvg_file,
                                                 out_vcf_file=sub_vcf_file,
                                                 rerun=self.smartrerun,
                                                 cmm=self.cmm))

        for p in processes:
            p.start()

        # listen for signal while any process is alive
        while True in [p.is_alive() for p in processes]:
            try:
                time.sleep(1)

            except KeyboardInterrupt:
                sys.stderr.write('KeyboardInterrupt detected, terminating '
                                 'all processes...\n')
                for p in processes:
                    p.terminate()

                sys.exit(1)

        # Make sure all process are finished
        for p in processes:
            p.join()

        # Final output file name
        utils.merge_files(out_cvg_names, self.outcvg, is_del_raw_file=True)
        if self.outvcf:
            utils.merge_files(out_vcf_names, self.outvcf, is_del_raw_file=True)

        return processes


# class VQSRRuner(object):
#     """Runner for VQSR"""
#     def __init__(self):
#         """Init function"""
#         self.vqsr = vqsr
#         return
#
#     def run(self):
#         self.vqsr.main(self.vqsr.cmdopts())
#
#         return


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
        sys.stderr.write('[INFO] Finish loading parameters and input file '
                         'list %s\n' % time.asctime())

    def run(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call varaintis by BaseType ... %s\n' %
                         time.asctime())

        # Always create process manager even if nCPU==1, so that we can
        # listen for signals from main thread
        regions_for_each_process = [[] for _ in range(self.nCPU)]
        if len(self.regions) < self.nCPU:
            # We cut the region evenly to fit nCPU if regions < nCPU
            for chrid, start, end in self.regions:
                delta = int((end - start + 1) / self.nCPU)
                if delta == 0:
                    delta = 1

                for i, pos in enumerate(xrange(start - 1, end, delta)):
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
            processes.append(CvgMultiProcess(self.referencefile,
                                             self.alignefiles,
                                             sub_cvg_file,
                                             regions_for_each_process[i],
                                             cmm=self.cmm))

        for p in processes:
            p.start()

        # listen for signal while any process is alive
        while True in [p.is_alive() for p in processes]:
            try:
                time.sleep(1)

            except KeyboardInterrupt:
                sys.stderr.write('KeyboardInterrupt detected, terminating '
                                 'all processes...\n')
                for p in processes:
                    p.terminate()

                sys.exit(1)

        # make sure all process are finished
        for p in processes:
            p.join()

        # Final output file name
        utils.merge_files(out_cvg_names, self.outputfile, is_del_raw_file=True)

        return processes


class MergeRunner(object):
    """Runner for merging files"""

    def __init__(self, args):
        """init function"""

        self.outputfile = args.outputfile
        # Load all files
        self.inputfiles = utils.load_file_list(args.infilelist)

    def run(self):
        utils.merge_files(self.inputfiles, self.outputfile)
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
