"""
This module contain some main algorithms of BaseVar
"""
import numpy as np
from scipy.stats import fisher_exact


def strand_bias(ref_base, alt_base, sample_base, strands):
    """
    A method for calculating the strand bias of REF_BASE and ALT_BASE

    :param ref_base:  char, required
        The reference base

    :param alt_base: array like, required
        A list of alt bases.

    :param sample_base: array-like, required
        A list of bases which cover this position

    :param strands: array-like, equired
        '+' or '-' strand for each base in ``sample_base``

    :return: list-like
        FS, ref_fwd, ref_rev, alt_fwd, alt_rev
    """
    ref_fwd, ref_rev = 0, 0
    alt_fwd, alt_rev = 0, 0
    for k, b in enumerate(sample_base):

        # For strand bias
        if b == 'N': continue
        if strands[k] == '+':
            if b == ref_base.upper():
                ref_fwd += 1
            elif b in alt_base:
                alt_fwd += 1

        elif strands[k] == '-':
            if b == ref_base.upper():
                ref_rev += 1
            elif b in alt_base:
                alt_rev += 1

    # Strand bias by fisher exact test
    # Normally you remove any SNP with FS > 60.0 and an indel with FS > 200.0
    fs = round(-10 * np.log10(
        fisher_exact([[ref_fwd, ref_rev], [alt_fwd, alt_rev]])[1]), 3)

    return fs, ref_fwd, ref_rev, alt_fwd, alt_rev
