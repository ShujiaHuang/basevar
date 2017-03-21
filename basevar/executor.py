"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``BaseVar.py`` is "Sauron", and this module could just be called by it.
"""

import sys
import argparse
import pysam
import numpy as np
from scipy.stats import fisher_exact

import utils
import mpileup
from variantcaller import BaseType

class RunBaseType(object):

    def __init__(self, cmm=utils.CommonParameter()):
        """init function
        """
        self.cmm = cmm

        optp = argparse.ArgumentParser()
        optp.add_argument('basetype')
        optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                          help='skip unlisted positions (chr pos)', default='')
        optp.add_argument('-R', '--regions',
                          metavar='chr:start-end', dest='regions',
                          help='skip positions not in (chr:start-end)', default='')
        optp.add_argument('-l', '--mpileup-list', dest='infilelist', metavar='FILE',
                          help='The input mpileup file list.', default='')
        optp.add_argument('-s', '--sample-list', dest='samplelistfile',
                          metavar='FILE', help='The sample list.')
        optp.add_argument('-o', '--outprefix', dest='outprefix',
                          metavar='FILE', default='out',
                          help='The prefix of output files. [out]')
        optp.add_argument('-q', '--basequality', dest='basequality',
                          metavar='INT', default=5,
                          help='The minine base quality threshold [5]')

        opt = optp.parse_args()
        self.basequality_threshold = int(opt.basequality)
        self.out_vcf_file = opt.outprefix + '.vcf'
        self.out_cvg_file = opt.outprefix + '.cvg.tsv' # position coverage

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] At least one mpileup to input\n')

        if len(opt.samplelistfile) == 0:
            optp.error('[ERROR] Must input the sample\'s ID list file by (-s)')

        if len(opt.positions) == 0 and len(opt.regions) == 0:
            optp.error('[ERROR] The list of position (-L or -R) is required.\n')

        # Loading positions
        _sites = utils.get_list_position(opt.positions) if opt.positions else {}
        if len(opt.regions):
            chrid, reg = opt.regions.strip().split(':')
            reg = map(int, reg.split('-'))
            if chrid not in _sites: _sites[chrid] = []
            _sites[chrid].append([reg[0], reg[1]])

        self.sites = {k:utils.merge_region(v) for k, v in _sites.items()}

        # Load all the mpileup files
        self.files = [f for f in sys.argv if '.mpileup.gz' in f]
        if opt.infilelist:
            self.files.extend(utils.load_file_list(opt.infilelist))

        # open mpileup files which store a batch of sample info by tabix
        self.tb_files = [pysam.TabixFile(f) for f in self.files]

        # load all the samples, 2D array
        self.sample_id = []
        with open(opt.samplelistfile) as I:

            samplefiles = []
            for r in I:
                if r[0] == '#': continue
                samplefiles.append(r.strip().split()[0])

            for f in samplefiles:
                with open(f) as I:
                    self.sample_id.append([s.strip().split()[0] for s in I])

    def run(self):

        total_sample = []
        _ = [total_sample.extend(s) for s in self.sample_id]
        vcf_header = utils.vcf_header_define()

        with open(self.out_vcf_file, 'w') as VCF, open(self.out_cvg_file, 'w') as CVG:

            CVG.write('\t'.join(['#CHROM','POS','Depth'] +
                                self.cmm.BASE) + '\n')

            VCF.write('\n'.join(vcf_header) + '\n')
            VCF.write('\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
                                 'INFO\tFORMAT'] + total_sample) + '\n')

            for chrid, regions in sorted(self.sites.items(), key = lambda x:x[0]):
                # ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
                # fetch the position data from each mpileup files
                # `iter_tokes` is a list of iterator for each sample's mpileup
                tmp_region = []
                for p in regions: tmp_region.extend(p)
                tmp_region = sorted(tmp_region)

                start, end = tmp_region[0], tmp_region[-1]
                iter_tokes = []
                for i, tb in enumerate(self.tb_files):
                    try:
                        iter_tokes.append(tb.fetch(chrid, start-1, end))
                    except ValueError:
                        if self.cmm.debug:
                            print >> sys.stderr, ("# [WARMING] Empty region",
                                                  chrid, start-1, end,
                                                  self.files[i])
                        iter_tokes.append('')

                # Set iteration marker: 1->iterate; 0->donot
                # iterate or hit the end
                go_iter = [1] * len(iter_tokes)
                for start, end in regions:
                    for position in xrange(start, end+1):

                        sample_info = [mpileup.fetch_next(iter_tokes[i])
                                       if g else sample_info[i]
                                       for i, g in enumerate(go_iter)]

                        sample_base_qual = []
                        sample_base = []
                        strands = []
                        ref_base = ''
                        for i, sample_line in enumerate(sample_info):

                            sample_info[i], ref_base_t, bs, qs, strand, go_iter[i] = (
                                mpileup.seek_position(position, sample_line,
                                                      len(self.sample_id[i]),
                                                      iter_tokes[i])
                            )

                            sample_base.extend(bs)
                            strands.extend(strand)
                            sample_base_qual.extend([ord(q) - 33 for q in qs])

                            if not ref_base:
                                ref_base = ref_base_t

                        bt = BaseType(ref_base.upper(), sample_base,
                                      sample_base_qual)
                        bt.lrt()

                        # ACGT count and mark the refbase
                        CVG.write('\t'.join(
                            [chrid, str(position), str(int(bt.total_depth))] +
                            [str(bt.depth[b])
                             if b!=ref_base.upper() else str(bt.depth[b]) + '*'
                             for b in self.cmm.BASE]) + '\n')

                        if len(bt.alt_bases()) > 0:
                            #VCF.write('\n')
                            self._out_vcf_line(chrid, position, ref_base,
                                               sample_base, strands, bt, VCF)

        self._close_tabix()

        return

    def _out_vcf_line(self, chrid, position, ref_base, sample_base,
                      strands, bt, out_file_handle):
        #  
        alt_gt = {b:'./'+str(k+1) for k,b in enumerate(bt.alt_bases())}
        samples = []

        ref_fwd, ref_rev = 0, 0
        alt_fwd, alt_rev = 0, 0

        for k, b in enumerate(sample_base):

            # For sample FORMAT
            if b != 'N':
                # For the base which not in bt.alt_bases()
                if b not in alt_gt: alt_gt[b] = './.'
                gt = '0/.' if b==ref_base.upper() else alt_gt[b]

                ## bt.qual_pvalue[k] is the base quality of 'b'
                samples.append(gt+':'+b+':'+strands[k]+':'+
                               str(round(bt.qual_pvalue[k], 6)))
            else:
                samples.append('./.') ## 'N' base

            # For strand bias
            if b == 'N': continue
            if strands[k] == '+':
                if b == ref_base.upper():
                    ref_fwd += 1
                else:
                    alt_fwd += 1

            elif strands[k] == '-':
                if b == ref_base.upper():
                    ref_rev += 1
                else:
                    alt_rev += 1

        # Strand bias by fisher exact test
        # Normally you remove any SNP with FS > 60.0 and an indel with FS > 200.0
        fs = round(-10 * np.log10(fisher_exact(
                   [[ref_fwd, ref_rev],[alt_fwd, alt_rev]])[1]), 3)

        # base=>[AF, allele depth]
        af = {b:['%f' % round(bt.depth[b]/float(bt.total_depth), 6),
                 bt.depth[b]] for b in bt.alt_bases()}

        info = {'CM_DP': str(int(bt.total_depth)),
                'CM_AC': ','.join(map(str, [af[b][1] for b in bt.alt_bases()])),
                'CM_AF': ','.join(map(str, [af[b][0] for b in bt.alt_bases()])),
                'CM_EAF': ','.join(map(str, [bt.eaf[b] for b in bt.alt_bases()])),
                'FS': str(fs),
                'SB_REF': str(ref_fwd)+','+str(ref_rev),
                'SB_ALT': str(alt_fwd)+','+str(alt_rev)}

        out_file_handle.write('\t'.join([chrid, str(position), '.', ref_base,
                         ','.join(bt.alt_bases()), str(bt.var_qual()),
                         '.' if bt.var_qual() > self.cmm.QUAL_THRESHOLD else 'LowQual',
                         ';'.join([k+'='+v for k, v in sorted(
                            info.items(), key=lambda x:x[0])]),
                            'GT:AB:SO:BP'] + samples) + '\n')
        return

    def _close_tabix(self):
        for tb in self.tb_files:
            tb.close()


