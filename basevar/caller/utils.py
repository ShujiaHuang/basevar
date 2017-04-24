"""
"""
import sys
import os
import gzip


class CommonParameter(object):
    """
    defined some globle common parameters
    """
    def __init__(self):
        self.LRT_THRESHOLD = 24   ## 24 corresponding to a chi-pvalue of 10^-6
        self.QUAL_THRESHOLD = 60  ## -10 * lg(10^-6)
        self.MLN10TO10 = -0.23025850929940458 # -np.log(10)/10
        self.BASE = ['A', 'C', 'G', 'T']
        self.BASE2IDX = {'A':0, 'C':1, 'G':2, 'T':3}
        self.debug = False
        self.MINAF = 0.001  # The effective base freqence threshold for 140k sample size


def vcf_header_define():
    header=['##fileformat=VCFv4.2',
            '##FILTER=<ID=LowQual,Description="Low quality">',
            ('##INFO=<ID=CM_EAF,Number=.,Type=Float,Description="An ordered, '
             'comma delimited list of Expectation Maxmization Estimated allele frequency">'),
            ('##INFO=<ID=CM_AF,Number=.,Type=Float,Description='
            '"An ordered, comma delimited list of allele frequencies base on read count">'),
            '##INFO=<ID=CM_AC,Number=.,Type=Float,Description="An ordered, comma delimited allele depth in CMDB">',
            '##INFO=<ID=CM_DP,Number=.,Type=Float,Description="Total Depth in CMDB">',
            '##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">',
            '##INFO=<ID=SB_REF,Number=.,Type=Integer,Description="Read number support REF: Forward,Reverse">',
            '##INFO=<ID=SB_ALT,Number=.,Type=Integer,Description="Read number support ALT: Forward,Reverse">']
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.append('##FORMAT=<ID=AB,Number=1,Type=String,Description="Allele Base">')
    header.append('##FORMAT=<ID=BP,Number=1,Type=String,Description="Base Probability which calculate by base quality">')
    header.append('##FORMAT=<ID=SO,Number=1,Type=String,Description="Strand orientation of the mapping base. Marked as + or -">')

    return header


def fetch_next(iter_fh):
    """
    re-define the next funtion in fetch function of pysam TabixFile()
    prevent throunghing the 'StopIteration'
    """

    if iter_fh == '': return ''

    try:
        line = iter_fh.next()
    except StopIteration:
        line = ''

    return line


def load_file_list(in_file):

    with open(in_file) as fh:
        files = [r.strip().split()[0] for r in fh if r[0] != '#']

    return files


def get_list_position(in_site_file):
    sites = {}
    with open(in_site_file) as f:
        for r in f:
            tok = r.strip().split()
            if tok[0] not in sites: sites[tok[0]] = []
            if len(tok) < 3:
                sites[tok[0]].append([int(tok[1]), int(tok[1])])
            else:
                sites[tok[0]].append([int(tok[1]), int(tok[2])])

    return sites


def merge_region(position_region, delta=1):
    """Merge a batch of sorted region

    Parameters
    ----------
    ``position_region``: a list like, required
        A regions (2D) array, format like: [[start1,end1], [start2,end2], ...]

    ``delta``: Integer, optinal

    Example
    -------
    ...
    >>> from basevar import utils
    >>> utils.merge_region([[1,1], [2,3], [4,6], [4,5], [8, 20], [9, 12]])
    ... [[1, 6], [8, 20]]

    """

    # sorted
    position_region.sort(key=lambda x:x[0])

    m_region = []
    prepos, start, end = '', '', ''
    flag = False
    for s, e in position_region:
        if s > e:
            print >> sys.stderr,('[ERROR]Your region start > end.'
                                 ' It is not allow when call Merge function\n')
            sys.exit(1)

        if prepos == '':
            # # The light is on => Get the region!
            if flag: m_region.append([start, end])
            start, end = s, e
            flag = True

        else:
            if prepos > s:
                print >> sys.stderr,('[ERROR]The array hasn\'t been sorted.\n')
                sys.exit(1)

            if delta + end >= s:
                end = e if e > end else end

            else:
                m_region.append([start, end])
                start, end = s, e

        prepos = s

    if flag:
        m_region.append([start, end])

    return m_region


def get_minor_major(base):
    """
    ``base`` is a [ATCG] string
    """
    bc = {}
    for b in base:
        b = b.upper()
        if b != 'N':
            bc[b] = bc.get(b, 0) + 1

    # It's a 2D list of tuple after sorting
    sorted_bc = sorted(bc.items(), key = lambda x:x[1], reverse=True)

    if len(sorted_bc):
        mi, mj = sorted_bc[-1][0], sorted_bc[0][0]
    else:
        mi, mj = 'N', 'N'

    # minor and major
    return mi, mj, sorted_bc


def expandedOpen(path, mode):
    try:
        return open(path, mode)
    except IOError:
        return open(os.path.expanduser(path), mode)


def Open(fileName, mode, compressLevel=9):
    """
    Function that allows transparent usage of dictzip, gzip and
    ordinary files
    """
    if fileName.endswith(".gz") or fileName.endswith(".GZ"):
        fileDir = os.path.dirname(fileName)
        if os.path.exists(fileDir):
            return gzip.GzipFile(fileName, mode, compressLevel)
        else:
            return gzip.GzipFile(os.path.expanduser(fileName), mode,
                                 compressLevel)
    else:
        return expandedOpen(fileName, mode)


