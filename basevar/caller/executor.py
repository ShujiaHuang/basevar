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

from . import utils
from .variantcaller import BaseVarMultiProcess


class BaseTypeRunner(object):

    def __init__(self, cmm=utils.CommonParameter()):
        """init function
        """
        optp = argparse.ArgumentParser()
        optp.add_argument('basetype')
        optp.add_argument('-o', '--outprefix', dest='outprefix', metavar='FILE', default='out',
                          help='The prefix of output files. [out]')
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
        optp.add_argument('-m', '--min_af', dest='min_af', type=float, metavar='MINAF', default=0.001,
                          help='By setting min AF to skip uneffective caller positions to accelerate program speed. [0.001]')

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
        _sites = utils.get_list_position(opt.positions) if opt.positions else {}
        if len(opt.regions):
            chrid, reg = opt.regions.strip().split(':')
            reg = map(int, reg.split('-'))
            if chrid not in _sites:
                _sites[chrid] = []

            _sites[chrid].append([reg[0], reg[1]])

        # merge and sorted the regions
        # [[chrid1, start1, end1], [chrid2, start2, end2], ...]
        self.regions = []
        for chrid, v in sorted(_sites.items(), key=lambda x: x[0]):
            for start, end in utils.merge_region(v):
                self.regions.append([chrid, start, end])

        # Load all the mpileup files
        self.mpileupfiles = [f for f in sys.argv if '.mpileup.gz' in f]
        if opt.infilelist:
            self.mpileupfiles.extend(utils.load_file_list(opt.infilelist))

        print >> sys.stderr, '[INFO] Finish loading parameters and mpileup list %s'%time.asctime()

    def run(self):
        """
        Run variant caller
        """

        print >> sys.stderr, '[INFO] Start call varaintis by BaseType ... %s'%time.asctime()
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
            sub_vcf_file = self.opt.outprefix + '_temp_%s'%i + '.vcf'
            sub_cvg_file = self.opt.outprefix + '_temp_%s'%i + '.cvg.tsv'

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
                print >> sys.stderr, 'KeyboardInterrupt detected, terminating all processes...'
                for p in processes:
                    p.terminate()

                sys.exit(1)

        # make sure all process are finished
        for p in processes:
            p.join()

        # Final output file name
        out_vcf_file = self.opt.outprefix + '.vcf'
        out_cvg_file = self.opt.outprefix + '.cvg.tsv'  # position coverage

        utils.merge_files(out_vcf_names, out_vcf_file)
        utils.merge_files(out_cvg_names, out_cvg_file)

        return


class MergeRunner(object):
    """Runner for merging files"""
    def __init__(self):
        """init function"""
        optp = argparse.ArgumentParser()
        optp.add_argument('merge')
        optp.add_argument('-l', '--file-list', dest='infilelist', metavar='FILE',
                          help='The input mpileup file list.', default='')
        optp.add_argument('-o', '--outfile', dest='outfile', metavar='FILE', default='out',
                          help='The prefix of output files. [out]')

        opt = optp.parse_args()
        self.opt = opt

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] At least input one mpileup file\n')

        # Load all files
        self.files = []
        if opt.infilelist:
            self.files.extend(utils.load_file_list(opt.infilelist))


    def run(self):
        utils.merge_files(self.files, self.opt.outfile)

        return
