"""
Author: Shujia Huang
Date: 2016-12-06
"""
import sys
import argparse

def main():
    optp = argparse.ArgumentParser()
    optp.add_argument('-s', '--sample2province', metavar='FILE',
                      dest='samplelist',
                      help='The input mpileup file list.')
    optp.add_argument('-i', '--input', metavar='FILE', dest='infile',
                      help='The input files.')

    opt = optp.parse_args()

    sample2prov = {}
    prov_set = set()
    with open(opt.samplelist) as I:
        for r in I:
            # xk54r1w51s      Guangdong
            if r[0] == '#': continue
            col = r.strip().split()
            sample2prov[col[0]] = col[1]

            if col[1] != 'null':
                prov_set.add(col[1])

    idx2sample = {}
    print '\t'.join(['#CHROM', 'POS', 'REF', 'ALT'] +
                    [pr for pr in sorted(prov_set, key=lambda x:x)])

    with open(opt.infile) as I:
        for r in I:
            if r[:2] == '##': continue

            col = r.strip().split()
            if col[0] == '#CHROM':
                idx2sample = {i:s for i, s in enumerate(col[9:])}
            else:
                k = col[0] + ':' + col[1]

                alt = col[4]  # ALT
                pprob = {}
                for i, s in enumerate(col[9:]):  ## loop the sample

                    _, ale, _, prob = s.split(':')  # GT:AB:SO:BP
                    if ale == 'N': continue

                    prob = float(prob)
                    prov = sample2prov[idx2sample[i]] #+ '-' + k

                    if prov == 'null': continue
                    if prov not in pprob: pprob[prov] = [0.0, 0.0]
#pprob[prov][0] += prob if ale == alt else 0.0 # (1.0 - prob)/3.0
                    pprob[prov][0] += prob if ale == alt else (1.0 - prob)/3.0
                    pprob[prov][1] += 1.0

                print '\t'.join([col[0], col[1], col[3], col[4]] +
                                [str('%f' % (round(pprob[pr][0]/pprob[pr][1], 6)))
                                 if pr in pprob else '0.0'
                                 for pr in sorted(prov_set, key=lambda x:x)])


if __name__ == '__main__':
    main()




