"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``BaseVar.py`` is "Sauron", and 'executor.py' module could just be called by it.
"""
from __future__ import division

import os
import sys
import argparse
import time

from pysam import AlignmentFile, TabixFile

from . import utils
from .fusion import Fusion
from .basetypebam import BaseVarMultiProcess as BamBaseVarMultiProcess
from .basetypefusion import BaseVarFusionMultiProcess
from .coverageprocess import CvgMultiProcess
# from .vqsr import vqsr


class BaseTypeBamRunner(object):

    def __init__(self, args, cmm=utils.CommonParameter()):
        """init function
        """
        # setting parameters
        self.referencefile = args.referencefile
        self.nCPU = args.nCPU
        self.pop_group_file = args.pop_group_file
        self.mapq = args.mapq
        self.batchcount = args.batchcount

        self.outvcf = args.outvcf if args.outvcf else None
        self.outcvg = args.outcvg

        self.smartrerun = True if args.smartrerun else False
        if self.smartrerun:
            sys.stderr.write("***********************************************\n"
                             "******************* WARNING *******************\n"
                             "***********************************************\n"
                             ">>>>>>>> You have setted `smart rerun` <<<<<<<<\n"
                             "Please make sure that all the parameters are the "
                             "same with your previous commands.\n"
                             "***********************************************\n\n")

        # Loading positions or load all the genome regions
        self.regions = utils.load_target_position(self.referencefile, args.positions, args.regions)

        # Get all the input alignement files
        if not args.input and not args.infilelist:
            sys.stderr.write("[ERROR] Missing input BAM/CRAM files.\n\n")
            sys.exit(1)

        self.alignefiles = args.input
        if args.infilelist:
            self.alignefiles += utils.load_file_list(args.infilelist)

        # setting the resolution of MAF
        self.cmm = cmm
        if args.min_af is None:
            args.min_af = min(100.0/len(self.alignefiles), 0.001, self.cmm.MINAF)

        self.cmm.MINAF = args.min_af
        sys.stderr.write('[INFO] Finish loading parameters and we have %d BAM/CRAM files '
                         'for variants calling %s\n' % (len(self.alignefiles), time.asctime()))

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

            processes.append(BamBaseVarMultiProcess(self.referencefile,
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

        return


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

        return


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


class BaseTypeFusionRunner(object):
    def __init__(self, cmm=utils.CommonParameter()):
        """init function
        """
        optp = argparse.ArgumentParser()
        optp.add_argument('basetypefusion')
        optp.add_argument('-I', '--fusion-file-list', dest='infilelist', metavar='FILE',
                          help='Fusion file list, one line per file.', default='')
        optp.add_argument('-R', '--reference', dest='referencefile', metavar='FILE',
                          help='Input reference fasta file.', default='')
        optp.add_argument('-O', '--outprefix', dest='outprefix', metavar='FILE',
                          default='out', help='The prefix of output files. [out]')

        optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                          help='skip unlisted positions (chr pos). [None]', default='')
        optp.add_argument('--region', metavar='chr:start-end', dest='region',
                          help='Skip position which not in these regions. Comma delimited '
                               'list of regions (chr:start-end). Could be a file contain the '
                               'regions.', default='')

        optp.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int,
                          help='Number of processer to use. [1]', default=1)
        optp.add_argument('-m', '--min_af', dest='min_af', type=float, metavar='MINAF',
                          help='By setting min AF to skip uneffective caller positions '
                               'to accelerate program speed. Usually you can set it to '
                               'be min(0.001, 100/x), x is the size of your population.'
                               '[min(0.001, 100/x)]')

        # special parameter for calculating specific population allele frequence
        optp.add_argument('--pop-group', dest='pop_group_file', metavar='FILE', type=str,
                          help='Calculating the allele frequency for specific population.')

        # special parameter to limit the function of BaseType
        optp.add_argument('--justdepth', dest='justdepth', type=bool,
                          help='Just output the depth information for each position [False]',
                          default=False)

        opt = optp.parse_args()
        self.opt = opt

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] Missing bamfile.\n')

        if len(opt.referencefile) == 0:
            optp.error('[ERROR] Missing reference fasta file.\n')

        # Loading positions if not provid we'll load all the genome
        self.regions = utils.load_target_position(opt.referencefile, opt.positions,
                                                  opt.region)

        # Get all the input align fusion files
        self.fusionfiles = utils.load_file_list(opt.infilelist)

        self.cmm = cmm
        if self.opt.min_af is None:
            self.opt.min_af = min(100.0 / len(self.fusionfiles), 0.001, self.cmm.MINAF)

        # reset threshold of min allele frequence threshold by sample size
        self.cmm.MINAF = self.opt.min_af

        sys.stderr.write('[INFO] Finish loading parameters and input file '
                         'list %s\n' % time.asctime())

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.sample_id = self._load_sample_id()

    def _load_sample_id(self):
        """loading sample id in BAM/CRMA files from RG tag"""

        sys.stderr.write('[INFO] Start loading all samples\' id from alignment files\n')
        # loading sample'id from the header of fusion files
        sample_id = []
        for i, f in enumerate(self.fusionfiles):

            tf = TabixFile(f)
            try:
                # get sample ID: '##RG\tSM:SAMPLE_ID'
                header = [h for h in tf.header if h.startswith('##RG\tSM:')][0]
                sample_id.append(header.split(':')[-1])
                tf.close()

            except IndexError:

                sys.stderr.write('[ERROR] File header has no sample tag mark '
                                 'by "SM:", Please check %s!' % f)
                tf.close()
                sys.exit(1)

        sys.stderr.write('[INFO] Finish load all %d samples\' ID '
                         'from RG tag\n\n' % len(sample_id))
        return sample_id

    def run(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call variants by BaseType ... %s\n' %
                         time.asctime())

        # Always create process manager even if nCPU==1, so that we can
        # listen for signals from main thread
        regions_for_each_process = [[] for _ in range(self.opt.nCPU)]
        if len(self.regions) < self.opt.nCPU:
            # We cut the region into pieces to fit nCPU if regions < nCPU
            for chrid, start, end in self.regions:
                delta = int((end - start + 1) / self.opt.nCPU)
                if delta == 0:
                    delta = 1

                for i, pos in enumerate(xrange(start - 1, end, delta)):
                    s = pos + 1 if pos + 1 < end else end
                    e = pos + delta if pos + delta < end else end

                    regions_for_each_process[i % self.opt.nCPU].append([chrid, s, e])

        else:
            for i, region in enumerate(self.regions):
                regions_for_each_process[i % self.opt.nCPU].append(region)

        out_vcf_names = set()
        out_cvg_names = set()

        processes = []
        for i in range(self.opt.nCPU):
            sub_cvg_file = self.opt.outprefix + '_temp_%s' % i + '.cvg.tsv'
            out_cvg_names.add(sub_cvg_file)

            if not self.opt.justdepth:
                sub_vcf_file = self.opt.outprefix + '_temp_%s' % i + '.vcf'
                out_vcf_names.add(sub_vcf_file)
            else:
                sub_vcf_file = None

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s, %s]\n' % (i + 1, self.opt.nCPU, sub_vcf_file,
                                             sub_cvg_file))

            processes.append(BaseVarFusionMultiProcess(self.opt.referencefile,
                                                       self.fusionfiles,
                                                       self.opt.pop_group_file,
                                                       regions_for_each_process[i],
                                                       self.sample_id,
                                                       out_cvg_file=sub_cvg_file,
                                                       out_vcf_file=sub_vcf_file,
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
        out_cvg_file = self.opt.outprefix + '.cvg.tsv'  # position coverage
        utils.merge_files(out_cvg_names, out_cvg_file, is_del_raw_file=True)

        if not self.opt.justdepth:
            out_vcf_file = self.opt.outprefix + '.vcf'
            utils.merge_files(out_vcf_names, out_vcf_file, is_del_raw_file=True)

        return


class FusionRunner(object):

    def __init__(self):
        """
        init function
        """
        optp = argparse.ArgumentParser()

        optp.add_argument('fusion')
        optp.add_argument('-I', '--inputfile', dest='inbamfile', metavar='FILE',
                          help='BAM/CRAM file list, one line per file.', default='')

        optp.add_argument('-R', '--reference', dest='referencefile', metavar='FILE',
                          help='Input reference fasta file.', default='')

        optp.add_argument('-O', '--outputfile', dest='outfile', metavar='FILE',
                          default='out.fusion', help='Output fusion file. [out.fusion]')

        opt = optp.parse_args()
        self.opt = opt

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] Missing input BAM/CRAM file.\n')

        if len(opt.referencefile) == 0:
            optp.error('[ERROR] Missing reference fasta file.\n')

    def run(self):

        # Get alignment sample ID
        bf = AlignmentFile(self.opt.inbamfile)
        sample_id = bf.header['RG'][0]['SM']
        bf.close()

        callfusion = Fusion(self.opt.referencefile, self.opt.inbamfile)
        with open(self.opt.outfile, 'w') as OUT:
            OUT.write('##fileformat=Fusion_v1.0 and the coordinate is 0-base system\n')
            OUT.write('##RG\tSM:%s\n' % sample_id)
            OUT.write('\t'.join(['#CHROM', 'START', 'END', 'TYPE', 'MAPQ',
                                 'SO', 'Read_POS', 'BASE_QUAL']) + '\n')
            for fusion in callfusion.generate_fusion():
                info = '\t'.join(map(str, [fusion.chrid,
                                           fusion.start,
                                           fusion.end,
                                           fusion.alt,
                                           fusion.mapq,
                                           fusion.strand_orientation,
                                           fusion.read_first_position,
                                           fusion.base_quality]
                                     )
                                 )

                OUT.write(info + '\n')
