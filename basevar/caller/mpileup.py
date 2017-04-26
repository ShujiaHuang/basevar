"""
Package for parse mpileup
Author: Shujia Huang
Date : 2016-07-19 14:14:21
"""
import re
import numpy as np
import sys

from . import utils

REOBJ_RE_START_END = re.compile('\^\S|\$')
REOBJ_RE_INDEL = re.compile('[-+]\d+[ACGTacgtNn]+')

def rmStartEnd(bases):
    """
    remove start(`^`) and end(`$`) character

    Examples
    --------

    Base example

    >>> import mpileup
    >>> bases="...,$.$.$A,..A...,,,.,,...+4AGGC...-5GTCGG......,a,^F,^].^F,"
    >>> mpileup.rmStartEnd(bases)
    ... ...,..A,..A...,,,.,,...+4AGGC...-5GTCGG......,a,,.,
    """
    return REOBJ_RE_START_END.sub('', bases)


def rmIndel(bases):
    """
    remove indels in pileup string

    Examples
    --------

    >>> import mpileup
    >>> bases="...,$.$.$A,..A...,,,.,,...+4AGGC...-5GTCGG......,a,^F,^].^F,"
    >>> mpileup.removeIndel(bases)
    ... ...,$.$.$A,..A...,,,.,,............,a,^F,^].^F,

    """
    return REOBJ_RE_INDEL.sub('', bases)


def fetch_next(iter_fh):
    """
    re-define the next funtion in fetch function of pysam TabixFile()
    prevent throunghing the 'StopIteration'
    """
    return utils.fetch_next(iter_fh)


def fetch_base_by_position(position, sample_id, sample_info, go_iter, iter_tokes,
                            is_scan_indel=False):

    sample_base_qual = []
    sample_base = []
    strands = []
    indels = []
    ref_base = ''
    for i, tb_sample_line in enumerate(sample_info):

        sample_info[i], ref_base_t, bs, qs, strand, go_iter[i], indel = (
            seek_position(position, tb_sample_line, len(sample_id[i]),
                          iter_tokes[i], is_scan_indel=is_scan_indel)
        )

        sample_base.extend(bs)
        strands.extend(strand)
        sample_base_qual.extend([ord(q) - 33 for q in qs])
        indels.extend(indel)

        if not ref_base:
            ref_base = ref_base_t

    return ref_base, sample_base, sample_base_qual, strands, indels


def seek_position(target_pos, sample_line, sample_num, sample_tb_iter,
                  is_scan_indel=False):

    ref_base = ''
    indels = []
    bases, quals, strand = [], [], []

    for _ in xrange(sample_num):
        bases.append('N')
        quals.append('!')
        strand.append('.')

    go_iter_mark = 0  # 1->iterate; 0->donot iterate or hit the end
    if sample_line:
        # chr2    181748  c       2       .,      EA
        tmp = sample_line.split('\t')  ## Must be '\t'
        pos = int(tmp[1])
        if pos == target_pos: # The same position

            ref_base = tmp[2]
            go_iter_mark = 1  # keep iterate
            for i in xrange(sample_num):

                try:
                    if tmp[3*(i+1)] != '0' and tmp[3*(i+1)+1] != '*':
                       strand[i], bases[i], quals[i], indel = first_base(
                           tmp[2], tmp[3*(i+1)+1], tmp[3*(i+1)+2],
                           is_scan_indel=is_scan_indel)

                       indels.extend(indel)

                except IndexError:
                    print >> sys.stderr, "[WARNING] IndexError. SampleID:", i+1, sample_num, len(tmp)
                    print >> sys.stderr, sample_line, "\n", tmp

        elif pos < target_pos:

            while pos < target_pos:

                sample_line = utils.fetch_next(sample_tb_iter)
                if sample_line:
                    tmp = sample_line.split('\t')  ## Must be '\t'
                    pos = int(tmp[1])
                else:
                    # The end of file. Break the loop.
                    break

            if pos == target_pos:

                ref_base = tmp[2]
                go_iter_mark = 1
                for i in xrange(sample_num):

                    try:
                        if tmp[3*(i+1)] != '0' and tmp[3*(i+1)+1] != '*':
                            strand[i], bases[i], quals[i], indel = first_base(
                                tmp[2], tmp[3*(i+1)+1], tmp[3*(i+1)+2],
                                is_scan_indel=is_scan_indel)

                            indels.extend(indel)

                    except IndexError:
                        print >> sys.stderr, "[WARNING] IndexError. SampleID:", i + 1, sample_num, len(tmp)
                        print >> sys.stderr, sample_line, "\n", tmp

        else:
            # pos > target_pos
            go_iter_mark = 0

    return sample_line, ref_base, bases, quals, strand, go_iter_mark, indels


def scan_indel(bases):
    """For scanning indel from ``bases`` info in mpileup column

    :param bases:
    :return: array-like. Indel
    """
    indel = []
    search_pos = 0
    while True:
        match = REOBJ_RE_INDEL.search(bases, pos=search_pos)
        if not match:
            break

        search_pos = match.end()  # jump to next match start position
        indel.append(match.group(0))  # catch indel

    return indel


def first_base(ref_base, bases, quality, is_scan_indel=False):
    """Just get the first base for each sample.
    """
    indel = scan_indel(bases) if is_scan_indel else []

    # ignore the indels, '^' or '$'
    b = rmIndel(rmStartEnd(bases))
    idx = 0

    ret_base = ref_base if b[idx] in [',', '.'] else b[idx]

    # Forwarstrand => +; reverseStrand => -.
    strand = '-' if (b[idx] == ',' or b[idx].islower()) else '+'
    return strand, ret_base.upper(), quality[idx], indel


def best_base(ref_base, bases, quality, is_scan_indel=False):
    """Just get the best quality base for each sample.
    """
    indel = scan_indel(bases) if is_scan_indel else []

    # ignore the indels, '^' or '$'
    b = rmIndel(rmStartEnd(bases))
    idx = np.argmax(quality, axis=0).eval() # get the best quality index

    ret_base = ref_base if b[idx] in [',', '.'] else b[idx]

    # Forwarstrand => +; reverseStrand => -.
    strand = '-' if (b[idx] == ',' or b[idx].islower()) else '+'
    return strand, ret_base.upper(), quality[idx], indel


def shuffle_base(ref_base, bases, quality, is_scan_indel=False):
    """
    """
    indel = scan_indel(bases) if is_scan_indel else []

    # ignore the indels, '^' or '$'
    b = rmIndel(rmStartEnd(bases))
    idx = range(len(b))
    np.random.shuffle(idx)  # shuffle the index

    ret_base = ref_base if b[idx[0]] in [',', '.'] else b[idx[0]]

    # Forwarstrand => +; reverseStrand => -.
    strand = '-' if (b[idx[0]] == ',' or b[idx[0]].islower()) else '+'

    return strand, ret_base.upper(), quality[idx[0]], indel

