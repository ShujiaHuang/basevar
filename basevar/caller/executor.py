"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``BaseVar.py`` is "Sauron", and this module could just be called by it.
"""
from __future__ import division

import sys
import argparse
import time

from pysam import FastaFile

from . import utils
from .basetypeprocess import BaseVarMultiProcess
from .basetypebamprocess import BaseVarMultiProcess as BamBaseVarMultiProcess
from .coverageprocess import BaseVarMultiProcess as CvgBaseVarMultiProcess
from .vqsr import vqsr


class BaseTypeBamRunner(object):

    def __init__(self, cmm=utils.CommonParameter()):
        """init function
        """
        optp = argparse.ArgumentParser()
        optp.add_argument('basetypebam')
        optp.add_argument('-o', '--outprefix', dest='outprefix', metavar='FILE',
                          default='out', help='The prefix of output files. [out]')
        optp.add_argument('-l', '--aligne-file-list', dest='infilelist', metavar='FILE',
                          help='Input alignmernt file list.', default='')
        optp.add_argument('-r', '--reference', dest='referencefile', metavar='FILE',
                          help='Input reference fasta file.', default='')

        optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                          help='skip unlisted positions (chr pos)', default='')
        optp.add_argument('-R', '--regions', metavar='chr:start-end', dest='regions',
                          help='skip positions not in (chr:start-end)', default='')

        optp.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int,
                          help='Number of processer to use. [1]', default=1)
        optp.add_argument('-m', '--min_af', dest='min_af', type=float, metavar='MINAF',
                          default=0.001, help='By setting min AF to skip uneffective '
                                              'caller positions to accelerate program '
                                              'speed. [0.001]')

        ## special parameter to limit the function of BaseType
        optp.add_argument('--justdepth', dest='justdepth', type=bool,
                          help='Just calculate the depth distribution [False]',
                          default=False)

        opt = optp.parse_args()
        self.opt = opt
        self.cmm = cmm

        # reset threshold of init min allele frequence by read depth
        self.cmm.MINAF = self.opt.min_af

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] At least input one mpileup file.\n')

        if len(opt.referencefile) == 0:
            optp.error('[ERROR] Missing reference fasta file.\n')

        # Loading positions if not provid we'll load all the genome
        self.regions = self._loading_position(opt.positions, opt.regions)

        # Get all the input alignement files
        self.alignefiles = utils.load_file_list(opt.infilelist)
        sys.stderr.write('[INFO] Finish loading parameters and input file '
                         'list %s\n' % time.asctime())

    def _loading_position(self, position, region):

        # Loading positions
        _sites = utils.get_list_position(position) if position else {}

        if len(region):
            chrid, reg = region.strip().split(':')
            reg = map(int, reg.split('-'))
            if chrid not in _sites:
                _sites[chrid] = []

            _sites[chrid].append([reg[0], reg[1]])

        # merge and sorted the regions
        # [[chrid1, start1, end1], [chrid2, start2, end2], ...]
        regions = []
        for chrid, v in sorted(_sites.items(), key=lambda x: x[0]):
            for start, end in utils.merge_region(v):
                regions.append([chrid, start, end])

        # load all the genome if no position or regions provide
        if not regions:

            sys.stderr.write('[WARNINGS] Program will load all the genome cause '
                             'there is not any positions and regions provided.\n')
            fa = FastaFile(self.opt.referencefile)
            regions = [[ci, 1, fa.get_reference_length(ci)]
                       for ci in fa.references]
            fa.close()

        return regions

    def run(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call varaintis by BaseType ... %s\n' %
                         time.asctime())

        # Always create process manager even if nCPU==1, so that we can
        # listen for signals from main thread
        regions_for_each_process = [[] for _ in range(self.opt.nCPU)]
        if len(self.regions) < self.opt.nCPU:
            # We cut the region evenly to fit nCPU if regions < nCPU
            for chrid, start, end in self.regions:
                delta = int((end-start+1) / self.opt.nCPU)
                if delta == 0:
                    delta = 1

                for i, pos in enumerate(xrange(start-1, end, delta)):
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
            sub_vcf_file = self.opt.outprefix + '_temp_%s' % i + '.vcf'
            sub_cvg_file = self.opt.outprefix + '_temp_%s' % i + '.cvg.tsv'

            out_vcf_names.add(sub_vcf_file)
            out_cvg_names.add(sub_cvg_file)

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s, %s]\n' % (i+1, self.opt.nCPU, sub_vcf_file,
                                             sub_cvg_file))

            processes.append(BamBaseVarMultiProcess(self.opt.referencefile,
                                                    self.alignefiles,
                                                    sub_vcf_file,
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
        out_vcf_file = self.opt.outprefix + '.vcf'
        out_cvg_file = self.opt.outprefix + '.cvg.tsv'  # position coverage

        utils.merge_files(out_vcf_names, out_vcf_file, is_del_raw_file=True)
        utils.merge_files(out_cvg_names, out_cvg_file, is_del_raw_file=True)

        return


class BaseTypeRunner(object):

    def __init__(self, cmm=utils.CommonParameter()):
        """init function
        """
        optp = argparse.ArgumentParser()
        optp.add_argument('basetype')
        optp.add_argument('-o', '--outprefix', dest='outprefix', metavar='FILE',
                          default='out', help='The prefix of output files. [out]')
        optp.add_argument('-l', '--mpileup-list', dest='infilelist', metavar='FILE',
                          help='The input mpileup file list.', default='')
        optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                          help='skip unlisted positions (chr pos)', default='')
        optp.add_argument('-R', '--regions', metavar='chr:start-end', dest='regions',
                          help='skip positions not in (chr:start-end)', default='')
        optp.add_argument('-s', '--sample-list', dest='samplelistfile', metavar='FILE',
                          help='The sample list.')
        optp.add_argument('-S', '--subsample-list', dest='subsample', metavar='FILE',
                          help='Skip samples not in subsample-list, one sample per row.')
        optp.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int,
                          help='Number of processer to use. [1]', default=1)
        optp.add_argument('-m', '--min_af', dest='min_af', type=float, metavar='MINAF',
                          default=0.001, help='By setting min AF to skip uneffective '
                                              'caller positions to accelerate program '
                                              'speed. [0.001]')

        opt = optp.parse_args()
        self.opt = opt

        self.cmm = cmm
        # reset threshold of init min allele frequence by read depth
        self.cmm.MINAF = self.opt.min_af

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] At least input one mpileup file\n')

        if len(opt.samplelistfile) == 0:
            optp.error('[ERROR] Must input the sample\'s ID list file by (-s)')

        if len(opt.positions) == 0 and len(opt.regions) == 0:
            optp.error('[ERROR] The list of position (-L or -R) is required.\n')

        # Loading positions
        self.regions = self._loading_position(opt)

        # Load all the input files
        self.mpileupfiles = utils.load_file_list(opt.infilelist)
        sys.stderr.write('[INFO] Finish loading parameters and mpileup '
                         'list %s\n' % time.asctime())

    def _loading_position(self, opt):

        # Loading positions
        _sites = utils.get_list_position(opt.positions) if opt.positions else {}

        if len(opt.regions):
            chrid, reg = opt.regions.strip().split(':')
            reg = map(int, reg.split('-'))
            if chrid not in _sites: _sites[chrid] = []

            _sites[chrid].append([reg[0], reg[1]])

        # merge and sorted the regions
        # [[chrid1, start1, end1], [chrid2, start2, end2], ...]
        regions = []
        for chrid, v in sorted(_sites.items(), key=lambda x: x[0]):
            for start, end in utils.merge_region(v):
                regions.append([chrid, start, end])

        return regions

    def run(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call varaintis by BaseType ... %s' %
                         time.asctime())

        # Always create process manager even if nCPU==1, so that we can
        # listen for signals from main thread
        regions_for_each_process = [[] for _ in range(self.opt.nCPU)]
        if len(self.regions) < self.opt.nCPU:
            # We cut the region evenly to fit nCPU if regions < nCPU
            for chrid, start, end in self.regions:
                delta = int((end-start+1) / self.opt.nCPU)
                if delta == 0:
                    delta = 1

                for i, pos in enumerate(xrange(start-1, end, delta)):
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
            sub_vcf_file = self.opt.outprefix + '_temp_%s' % i + '.vcf'
            sub_cvg_file = self.opt.outprefix + '_temp_%s' % i + '.cvg.tsv'

            out_vcf_names.add(sub_vcf_file)
            out_cvg_names.add(sub_cvg_file)
            processes.append(BaseVarMultiProcess(self.mpileupfiles,
                                                 sub_vcf_file,
                                                 sub_cvg_file,
                                                 regions_for_each_process[i],
                                                 self.opt,
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
        out_vcf_file = self.opt.outprefix + '.vcf'
        out_cvg_file = self.opt.outprefix + '.cvg.tsv'  # position coverage

        utils.merge_files(out_vcf_names, out_vcf_file, is_del_raw_file=True)
        utils.merge_files(out_cvg_names, out_cvg_file, is_del_raw_file=True)

        return


class VQSRRuner(object):
    """Runner for VQSR"""
    def __init__(self):
        """Init function"""
        self.vqsr = vqsr
        return

    def run(self):
        self.vqsr.main(self.vqsr.cmdopts())

        return


class CoverageRunner(object):

    def __init__(self, cmm=utils.CommonParameter()):
        """init function
        """
        optp = argparse.ArgumentParser()
        optp.add_argument('coverage')
        optp.add_argument('-o', '--outprefix', dest='outprefix', metavar='FILE',
                          default='out', help='The prefix of output files. [out]')
        optp.add_argument('-l', '--aligne-file-list', dest='infilelist', metavar='FILE',
                          help='Input alignmernt file list.', default='')
        optp.add_argument('-r', '--reference', dest='referencefile', metavar='FILE',
                          help='Input reference fasta file.', default='')

        optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                          help='skip unlisted positions (chr pos)', default='')
        optp.add_argument('-R', '--regions', metavar='chr:start-end', dest='regions',
                          help='skip positions not in (chr:start-end)', default='')

        optp.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int,
                          help='Number of processer to use. [1]', default=1)

        opt = optp.parse_args()
        self.opt = opt
        self.cmm = cmm

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] At least input one mpileup file.\n')

        if len(opt.referencefile) == 0:
            optp.error('[ERROR] Missing reference fasta file.\n')

        # Loading positions if not provid we'll load all the genome
        self.regions = self._loading_position(opt.positions, opt.regions)

        # Get all the input alignement files
        self.alignefiles = utils.load_file_list(opt.infilelist)
        sys.stderr.write('[INFO] Finish loading parameters and input file '
                         'list %s\n' % time.asctime())

    def _loading_position(self, position, region):

        # Loading positions
        _sites = utils.get_list_position(position) if position else {}

        if len(region):
            chrid, reg = region.strip().split(':')
            reg = map(int, reg.split('-'))
            if chrid not in _sites:
                _sites[chrid] = []

            _sites[chrid].append([reg[0], reg[1]])

        # merge and sorted the regions
        # [[chrid1, start1, end1], [chrid2, start2, end2], ...]
        regions = []
        for chrid, v in sorted(_sites.items(), key=lambda x: x[0]):
            for start, end in utils.merge_region(v):
                regions.append([chrid, start, end])

        # load all the genome if no position or regions provide
        if not regions:

            sys.stderr.write('[WARNINGS] Program will load all the genome cause '
                             'there is not any positions and regions provided.\n')
            fa = FastaFile(self.opt.referencefile)
            regions = [[ci, 1, fa.get_reference_length(ci)]
                       for ci in fa.references]
            fa.close()

        return regions

    def run(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call varaintis by BaseType ... %s\n' %
                         time.asctime())

        # Always create process manager even if nCPU==1, so that we can
        # listen for signals from main thread
        regions_for_each_process = [[] for _ in range(self.opt.nCPU)]
        if len(self.regions) < self.opt.nCPU:
            # We cut the region evenly to fit nCPU if regions < nCPU
            for chrid, start, end in self.regions:
                delta = int((end-start+1) / self.opt.nCPU)
                if delta == 0:
                    delta = 1

                for i, pos in enumerate(xrange(start-1, end, delta)):
                    s = pos + 1 if pos + 1 < end else end
                    e = pos + delta if pos + delta < end else end

                    regions_for_each_process[i % self.opt.nCPU].append([chrid, s, e])

        else:
            for i, region in enumerate(self.regions):
                regions_for_each_process[i % self.opt.nCPU].append(region)

        out_cvg_names = set()
        processes = []
        for i in range(self.opt.nCPU):
            sub_cvg_file = self.opt.outprefix + '_temp_%s' % i + '.cvg.tsv'

            out_cvg_names.add(sub_cvg_file)
            processes.append(CvgBaseVarMultiProcess(self.opt.referencefile,
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
        out_cvg_file = self.opt.outprefix + '.cvg.tsv'  # position coverage

        utils.merge_files(out_cvg_names, out_cvg_file, is_del_raw_file=True)

        return


class MergeRunner(object):
    """Runner for merging files"""
    def __init__(self):
        """init function"""
        optp = argparse.ArgumentParser()
        optp.add_argument('merge')
        optp.add_argument('-l', '--file-list', dest='infilelist', metavar='FILE',
                          help='The input files\' list.', default='')
        optp.add_argument('-o', '--outfile', dest='outfile', metavar='FILE', default='out',
                          help='The prefix of output files. [out]')

        opt = optp.parse_args()
        self.opt = opt

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] At least one input file.\n')

        # Load all files
        self.files = []
        if opt.infilelist:
            self.files.extend(utils.load_file_list(opt.infilelist))

    def run(self):
        utils.merge_files(self.files, self.opt.outfile)

        return


class NearbyIndelRunner(object):

    def __init__(self):
        """init function"""
        optp = argparse.ArgumentParser()
        optp.add_argument('nbi')
        optp.add_argument('-i', '--in-vcf-file', dest='in_vcf_file', metavar='FILE',
                          help='The input vcf files', default='')
        optp.add_argument('-c', '--in-cvg-file', dest='in_cvg_file', metavar='FILE',
                          help='Input coverage file which has indel information',
                          default='')
        optp.add_argument('-d', '--nearby-distance-around-indel',
                          dest='nearby_dis_around_indel',
                          help='The distance around indel. [16]', default=16)

        opt = optp.parse_args()
        opt.nearby_dis_around_indel = int(opt.nearby_dis_around_indel)
        self.opt = opt

        if len(sys.argv) == 2 and (len(opt.in_vcf_file) == 0 or len(opt.in_cvg_file) == 0):
            optp.error('[ERROR] At least one input file.\n')

        sys.stderr.write('[INFO] Parameters: python %s nbi' 
                         '\n\t-i %s'
                         '\n\t-c %s'
                         '\n\t-d %d' % (sys.argv[0],
                                        opt.in_vcf_file,
                                        opt.in_cvg_file,
                                        opt.nearby_dis_around_indel)
                         )

    def run(self):

        from .other import NearbyIndel

        nbi = NearbyIndel(self.opt.in_vcf_file,
                          self.opt.in_cvg_file,
                          nearby_distance=self.opt.nearby_dis_around_indel)
        nbi.run()

        return self
