"""
Calculation the selection factor of geographic

Author: Shujia Huang
Date: 2017-11-10
"""
from __future__ import division

import sys
import time
import gzip
import argparse

import numpy as np


from rpy2 import robjects
from rpy2 import rinterface


R_IntVector = robjects.IntVector
R_MATRIX = robjects.r['matrix']
R_Fisher_Test = robjects.r['fisher.test']

def get_pass_positon(in_file):

    sites = set()
    with open (in_file) as I:
        for r in I:
            tok = r.strip().split()
            sites.add(tok[0]+':'+tok[1])

    return sites

def get_list_position(in_site_file):

    sites = {}
    with open(in_site_file) as f:
        for r in f:
            # #CHROM  POS  ID  REF   ALT   DISEASE
            tok = r.strip().split()

            ## skip none SNP positions
            if len(tok[3]) > 1 or any(len(alt) > 1 for alt in tok[4].split(',')):
                continue

            pk = tok[0]+':'+tok[1]

            tok[1] = int(tok[1])
            tok[4] = tok[4].upper()  # convert allele to upper

            sites[pk] = tok

    return sites


def lookup_frequence_nearby_allele(freq, g_idx, nbf_slide_window, data,
                                   max_size=100000):

    low_freq = 0.9 * freq
    high_freq = 1.1 * freq

    tmp_nbf_ale = []
    if nbf_slide_window:
        first_bt_idx = len(nbf_slide_window)
        for i, d in enumerate(nbf_slide_window):
            if d[1] >= low_freq:
                first_bt_idx = i
                break

        # reset the data of nbf_slide_window
        tmp_nbf_ale = nbf_slide_window[first_bt_idx:]

    for i in xrange(g_idx, len(data)):

        if low_freq > data[i][1]: continue
        if high_freq < data[i][1]:
            g_idx = i  # record the last one
            break

        tmp_nbf_ale.append(data[i])

    idx = range(len(tmp_nbf_ale))
    np.random.shuffle(idx)  # random shuffling the index of nbf_slide_window

    nbf = [tmp_nbf_ale[i] for i in idx[0:max_size]] \
        if len(tmp_nbf_ale) > max_size else tmp_nbf_ale

    return g_idx, nbf, tmp_nbf_ale


def calculate_significant(nbf_data, have_fisher_test_res):
    """
    Fisher Exact test of nearby frequence positions
    """

    res = []
    for pos_key, f, alt_b, n, c, s in nbf_data:

        if pos_key not in have_fisher_test_res:
            north = map(int, map(float, n.split(':')[:2]))  # first element is REF base
            centr = map(int, map(float, c.split(':')[:2]))  # first element is REF base
            south = map(int, map(float, s.split(':')[:2]))  # first element is REF base

            depth = sum(north + centr + south)
            if depth:
                # User FisherExact to do Fisher exact test for mxn contingency table
                m = R_MATRIX(R_IntVector(north + centr + south), nrow=2)

                try:
                    pvalue = R_Fisher_Test(m, workspace=2e6)[0][0]
                    have_fisher_test_res[pos_key] = pvalue
                except rinterface.RRuntimeError:
                    sys.stderr.write('[SKIP] [workspace not enough(N:C:S)] %s\n' %
                                      '\t'.join([pos_key]+map(str, north + centr + south)))
                    continue

            else:
                have_fisher_test_res[pos_key] = 1.0

        res.append([pos_key, alt_b, have_fisher_test_res[pos_key]])

    return res


def get_rank(data, pos_key, alt_base):

    if not data:
        return 'NA', 'NA', 'NA'

    # sorted by pvalue
    data.sort(key=lambda x: x[1])

    n = 1
    for i, d in enumerate(data):

        # avoid multiple allele problem
        if d[0] == pos_key and d[1] == alt_base:
            n = i + 1
            break

    # pvalue, percentile_pvalue, percentile_rank
    return data[n-1][-1], round(float(n)/len(data), 6), "%d/%d" % (n, len(data))


if __name__ == '__main__':

    sys.stderr.write('[INFO] Program starting .. %s\n' % time.asctime())

    # program's parameters
    optp = argparse.ArgumentParser()
    optp.add_argument('-i', '--in_vcf_file', dest='in_file', metavar='FILE',
                      default='', help='Input files, required.')
    optp.add_argument('-n', '--nearfreqnum', metavar='NUM', dest='num', type=int,
                      help='Max number of allleles in one block', default=100000)
    optp.add_argument('-l', '--positions', metavar='FILE', dest='positions',
                      help='skip unlisted positions (chr pos)', default='')
    optp.add_argument('-p', '--pass_pos', metavar='FILE', dest='pass_pos',
                      help='skip all the position unlist here', default='')

    opt = optp.parse_args()

    if not opt.in_file:
        optp.error('[ERROR] The input file (-i) is required.\n')

    if not opt.positions:
        optp.error('[ERROR] The list of position (-l) is required.\n')

    sys.stderr.write('[INFO] Parameters: python %s '
                     '\n\t-i %s'
                     '\n\t-n %d'
                     '\n\t-p %s'
                     '\n\t-l %s\n\n' % (
                         sys.argv[0], opt.in_file, opt.num,
                         opt.pass_pos, opt.positions))

    target_sites = get_list_position(opt.positions)
    pass_sites = get_pass_positon(opt.pass_pos)

    # main function
    data = []
    have_fisher_test_res = {}
    with gzip.open(opt.in_file) \
            if opt.in_file.endswith(".gz") else open(opt.in_file) as I:

        for r in I:
            # chr22  34510633  C  A,G  0.0078,0.564  424.7:8.9:510.1  262.3:7.0:322.2  187.6:3.0:165.4  15.8:0:23.9
            if r.startswith('#'): continue
            tok = r.strip().split()

            if len(tok) != 9:
                sys.stderr.write('[ValueError] %s\n' % '\t'.join(tok))
                continue

            pos_key = tok[0] + ':' + tok[1]
            if pos_key not in pass_sites:
                continue

            alt_base = tok[3].split(',')
            alt_freq = map(float, tok[4].split(','))
            centr = map(int, map(float, tok[5].split(':'))) # fist one is REF base
            north = map(int, map(float, tok[6].split(':'))) # fist one is REF base
            south = map(int, map(float, tok[7].split(':'))) # fist one is REF base

            for i, (alt_f, alt_b) in enumerate(zip(alt_freq, alt_base)):

                # [pos_key, alele_freq, alt_base, north, central, south]
                data.append([pos_key,
                             alt_f,
                             alt_b.upper(),
                             str(north[0]) + ':' + str(north[i+1]) + ':' + str(round(float(north[1])/sum(north), 5)) if sum(north)>0 else 'NA',
                             str(centr[0]) + ':' + str(centr[i+1]) + ':' + str(round(float(centr[1])/sum(centr), 5)) if sum(centr)>0 else 'NA',
                             str(south[0]) + ':' + str(south[i+1]) + ':' + str(round(float(south[1])/sum(south), 5)) if sum(south)>0 else 'NA'])

    # sorted the whole array data from lowest to highest by allele_frequence
    data.sort(key=lambda x:x[1])

    # Now we search all our target position by frequence data
    g_idx = 0
    nbf_slide_window = []

    out_res = []

    print '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'DISEASE', 'ALT_freq',
                     'FisherExactTest:p-value', 'Percentile:p-value',
                     'Percentile:Rank','North(REF:ALT:AF)', 'Central(REF:ALT:AF)',
                     'South(REF:ALT:AF)'])

    sys.stderr.write('[INFO] ** Data is all loaded .. %s\n' % time.asctime())

    num = 0
    for d in data:
        # loop frequence from lowest to highest

        if num % 100000 == 0:
            sys.stderr.write('[INFO] ** Dealing lines %d %s\n' % (num, time.asctime()))

        num += 1

        pk, freq, alt = d[0:3]
        if pk not in target_sites:
            continue

        # skip un match allele postions
        if alt not in target_sites[pk][4]:
            continue

        ## catch the postione
        g_idx, nbf, nbf_slide_window = lookup_frequence_nearby_allele(
            freq, g_idx, nbf_slide_window, data, max_size=opt.num
        )

        pdata = calculate_significant(nbf, have_fisher_test_res)
        pvalue, percentile_pvalue, percentile_pos = get_rank(pdata, pk, alt)

        chr_id, pos = pk.split(':')
        pos = int(pos)

        out_res.append(target_sites[pk][:4] + [alt, target_sites[pk][5]] +
                       [freq, pvalue, percentile_pvalue, percentile_pos] + d[3:])

        sys.stderr.write('[OUT] %s\n' % '\t'.join(map(str, out_res[-1])))

    sys.stderr.write('[INFO] ** All lines %d done and outputting result .. %s\n' %
                     (num, time.asctime()))

    out_res.sort(key=lambda x: (x[0], x[1]))

    for d in out_res:
        print '\t'.join(map(str, d))












