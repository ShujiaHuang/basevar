"""
This module contain functions for calculating the genotype probability.
"""
import sys
import argparse
import pysam
import numpy as np

import utils


def BQprob():
    optp = argparse.ArgumentParser()
    optp.add_argument('bqprob')
    optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                      help='skip unlisted positions (chr pos)', default='')
    optp.add_argument('-l', '--mpileup-list', dest='infilelist',
                      metavar='FILE',
                      help='The input mpileup file list.', default='')
    optp.add_argument('-s', '--sample-list', dest='samplelistfile',
                      metavar='FILE',help='The sample list.')

    opt = optp.parse_args()

    # Loading positions
    sites = utils.get_list_position(opt.positions)

    # Load all the mpileup files
    files = [f for f in sys.argv if '.mpileup.gz' in f]
    if opt.infilelist:
        files.extend(utils.load_file_list(opt.infilelist))

    # open all the mpileup index by tabix and contant all the samples' id
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
    print >> sys.stderr, '\t'.join(['#CHROM','POS','Depth','A','C','G','T'])

    vcf_header = utils.vcf_header_define()
    print '\n'.join(vcf_header)
    print '\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'] +
                    total_sample)
    for chrid, positions in sorted(sites.items(), key = lambda x:x[0]):
        # fetch the position data from each mpileup files
        # `iter_tokes` is a list of iterator for each sample's mpileup
        start, end = positions[0], positions[-1]
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

        for start in positions:

            # fetch one line for initail.
            sample_info = [utils.fetch_next(iter_tokes[i]) if g else sample_info[i]
                           for i, g in enumerate(go_iter)]

            base_list = [chrid, str(start)]
            pvalue = []
            sample_base = []
            strands = []
            ref_base = ''
            for i, sample_line in enumerate(sample_info):

                sample_info[i], ref_base_tmp, bases, quals, strand, go_iter[i] = utils.seek_position(
                    start, sample_line, len(sample_id[i]), iter_tokes[i])

                _ = [sample_base.append(base) for base in bases]
                _ = [strands.append(s) for s in strand]

                q = np.array([ord(qual) - 33 for qual in quals])
                p_error = 10 ** (q/-10.0)
                pvalue.extend(1.0 - p_error)

                if not ref_base:
                    ref_base = ref_base_tmp.upper()

            # calculate the major and minor
            minor, major, base_depth_tuple = utils.get_minor_major(sample_base)
            if major == 'N':

                print >> sys.stderr, '\t'.join([chrid, str(start),
                    '0', '0', '0', '0', '0'])
            else:
                # ignore all 'N' positions
                base_depth_dict = dict(base_depth_tuple)
                print >> sys.stderr, '\t'.join([chrid, str(start),
                    str(sum([d[1] for d in base_depth_tuple]))] +
                    [str(base_depth_dict[b] if b in base_depth_dict else 0)
                     for b in ['A', 'C', 'G', 'T']]) # ACGT count

                mj, dj = major, base_depth_tuple[0][1]
                mi, di = mj, dj
                if len(base_depth_tuple) > 1:
                    mi, di = base_depth_tuple[1]

                # Finde alternate allele
                alt = []
                if mj != ref_base: alt.append(mj)
                if mi != ref_base: alt.append(mi)

                # ignore the reference line
                if not alt: continue

                samples = []
                af = {}
                n = 0
                for k, b in enumerate(sample_base):
                    samples.append('./.:'+b+':'+strands[k]+':'+str(round(pvalue[k], 6)))
                    if b != 'N':
                        n += 1.0
                        for ale in alt:
                            if b not in af: af[b] = [0.0, 0.0]  # expected allele count, allele depth
                            if b == ale:
                                pv = pvalue[k]
                                af[b][1] += 1
                            else:
                                pv = (1.0-pvalue[k]) / 3.0
                            af[b][0] += pv

                # base=>[expected count, AF, allele depth]
                af = {k:[round(f, 2), '%f' % round(f/n, 6), int(c)]
                      for k, (f, c) in af.items()}

                ref_fwd, ref_res = 0, 0
                alt_fwd, alt_res = 0, 0
                for k, s in enumerate(strands):
                    if s == 'F':
                        if sample_base[k] == ref_base.upper():
                            ref_fwd += 1
                        else:
                            alt_fwd += 1

                    elif s == 'R':
                        if sample_base[k] == ref_base.upper():
                            ref_res += 1
                        else:
                            alt_res += 1

                info = {'CM_DP': str(sum([d[1] for d in base_depth_tuple])),
                        'EAC': ','.join(map(str, [af[a][0] for a in alt])),
                        'CM_AC': ','.join(map(str, [af[a][2] for a in alt])),
                        'CM_AF': ','.join(map(str, [af[a][1] for a in alt])),
                        'SB_REF': str(ref_fwd)+','+str(ref_res),
                        'SB_ALT': str(alt_fwd)+','+str(alt_res)}

                print '\t'.join([chrid, str(start), '.', ref_base,
                                 ','.join(alt), '100', 'PASS',
                                 ';'.join([k+'='+v for k, v in sorted(
                                    info.items(), key=lambda x:x[0])]),
                                 'GT:AB:SO:BP'] + samples)

    for tb in tb_files:
        tb.close()

    return





