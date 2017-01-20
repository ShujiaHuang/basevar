"""
Package for parse mpileup
Author: Shujia Huang
Date : 2016-07-19 14:14:21
"""
import re

def rmStartEnd(bases):
    """
    remove start(`^`) and end(`$`) character

    Examples
    --------

    Base example

    >>> import mpileup
    >>> bases="...,$.$.$A,..A...,,,.,,...+5AGGC...-8GTCGG......,a,^F,^].^F,"
    >>> mpileup.clip(bases)
    ... ...,..A,..A...,,,.,,...+5AGGC...-8GTCGG......,a,,.,
    """
    return re.sub('\^\S|\$', '', bases)


def rmIndel(bases):
    """
    remove indels in pileup string

    Examples
    --------

    >>> import mpileup
    >>> bases="...,$.$.$A,..A...,,,.,,...+5AGGC...-8GTCGG......,a,^F,^].^F,"
    >>> mpileup.removeIndel(bases)
    ... ...,$.$.$A,..A...,,,.,,............,a,^F,^].^F,

    """
    return re.sub('[-+]\d+[ACGTacgtNn]+', '', bases)


