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
from variantcaller import BaseType


def basetype(cmm=utils.CommonParameter()):

    optp = argparse.ArgumentParser()
    optp.add_argument('basetype')
    optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                      help='skip unlisted positions (chr pos)', default='')
    optp.add_argument('-R', '--regions', metavar='Region', dest='regions',
                      help='skip positions not in (chr:start-end)', default='')
    optp.add_argument('-l', '--mpileup-list', dest='infilelist', metavar='FILE',
                      help='The input mpileup file list.', default='')
    optp.add_argument('-s', '--sample-list', dest='samplelistfile',
                      metavar='FILE', help='The sample list.')

    opt = optp.parse_args()

    # Loading positions
    _sites = utils.get_list_position(opt.positions) if opt.positions else {}
    if len(opt.regions):
        chrid, reg = opt.regions.strip().split(':')
        reg = map(int, reg.split('-'))
        if chrid not in _sites: _sites[chrid] = []
        _sites[chrid].append([reg[0], reg[1]])

    sites = {k:utils.merge_region(v) for k, v in _sites.items()}

    # Load all the mpileup files
    files = [f for f in sys.argv if '.mpileup.gz' in f]
    if opt.infilelist:
        files.extend(utils.load_file_list(opt.infilelist))

    # open mpileup files which store a batch of sample info by tabix
    tb_files = [pysam.TabixFile(f) for f in files]

    # load all the samples
    sample_id = []
    with open(opt.samplelistfile) as I:

        samplefiles = []
        for r in I:
            if r[0] == '#': continue
            samplefiles.append(r.strip().split()[0])

        for f in samplefiles:
            with open(f) as I:
                sample_id.append([s.strip().split()[0] for s in I])

    total_sample = []
    _ = [total_sample.extend(s) for s in sample_id]
    print >> sys.stderr, '\t'.join(['#CHROM','POS','Depth'] + cmm.BASE)

    vcf_header = utils.vcf_header_define()
    print '\n'.join(vcf_header)
    print '\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'] +
                    total_sample)

#print "# init done", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "\n"
    for chrid, regions in sorted(sites.items(), key = lambda x:x[0]):
        # ``positions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
        # fetch the position data from each mpileup files
        # `iter_tokes` is a list of iterator for each sample's mpileup
        tmp_region = []
        for p in regions: tmp_region.extend(p)
        tmp_region = sorted(tmp_region)

        start, end = tmp_region[0], tmp_region[-1]
#start, end = positions[0], positions[-1]
        iter_tokes = []
        for i, tb in enumerate(tb_files):
            try:
                iter_tokes.append(tb.fetch(chrid, start-1, end))
            except ValueError:
                print >> sys.stderr, ("# [WARMING] Empty region",
                                      chrid, start-1, end, files[i])
                iter_tokes.append('')

        # Set iteration marker: 1->iterate; 0->donot iterate or hit the end
        go_iter = [1] * len(iter_tokes)
        for start, end in regions:

            for position in xrange(start, end+1):

                # just fetch one line each time for initail.
                sample_info = [utils.fetch_next(iter_tokes[i]) if g else sample_info[i]
                               for i, g in enumerate(go_iter)]

                sample_base_qual = []
                sample_base = []
                strands = []
                ref_base = ''
                for i, sample_line in enumerate(sample_info):

                    sample_info[i], ref_base_tmp, bases, quals, strand, go_iter[i] = utils.seek_position(
                        position, sample_line, len(sample_id[i]), iter_tokes[i])

                    sample_base.extend(bases)
                    strands.extend(strand)
                    sample_base_qual.extend([ord(qual) - 33 for qual in quals])

                    if not ref_base:
                        ref_base = ref_base_tmp

                bt = BaseType(ref_base.upper(), sample_base, sample_base_qual)
                bt.lrt()

                print >> sys.stderr, '\t'.join([chrid, str(position), str(int(bt.total_depth))] +
                                               [str(bt.depth[b]) for b in cmm.BASE]) # ACGT count
                if len(bt.alt_bases()) > 0:

                    samples = []
                    for k, b in enumerate(sample_base):
                        if b != 'N':
                            ## bt.qual_pvalue[k] is the base quality of 'b'
                            samples.append('./.:'+b+':'+strands[k]+':'+
                                           str(round(bt.qual_pvalue[k], 6)))
                        else:
                            samples.append('./.') ## 'N' base

                    # base=>[AF, allele depth]
                    af = {b:['%f' % round(bt.depth[b]/float(bt.total_depth), 6),
                             bt.depth[b]] for b in bt.alt_bases()}

                    ref_fwd, ref_rev = 0, 0
                    alt_fwd, alt_rev = 0, 0
                    for k, s in enumerate(strands):

                        if sample_base[k] == 'N': continue

                        if s == '+':
                            if sample_base[k] == ref_base.upper():
                                ref_fwd += 1
                            else:
                                alt_fwd += 1

                        elif s == '-':
                            if sample_base[k] == ref_base.upper():
                                ref_rev += 1
                            else:
                                alt_rev += 1

                    fs = round(-10 * np.log10(fisher_exact(
                               [[ref_fwd, ref_rev],[alt_fwd, alt_rev]])[1]), 3)

                    info = {'CM_DP': str(int(bt.total_depth)),
                            'CM_AC': ','.join(map(str, [af[b][1] for b in bt.alt_bases()])),
                            'CM_AF': ','.join(map(str, [af[b][0] for b in bt.alt_bases()])),
                            'CM_EAF': ','.join(map(str, [bt.eaf[b] for b in bt.alt_bases()])),
                            'FS': str(fs),
                            'SB_REF': str(ref_fwd)+','+str(ref_rev),
                            'SB_ALT': str(alt_fwd)+','+str(alt_rev)}

                    print '\t'.join([chrid, str(position), '.', ref_base,
                                     ','.join(bt.alt_bases()), str(bt.var_qual()),
                                     '.' if bt.var_qual() > cmm.QUAL_THRESHOLD else 'LowQual',
                                     ';'.join([k+'='+v for k, v in sorted(
                                        info.items(), key=lambda x:x[0])]),
                                        'GT:AB:SO:BP'] + samples)

    for tb in tb_files:
        tb.close()

    return














