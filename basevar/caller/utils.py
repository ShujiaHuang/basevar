"""
"""
import sys
import os
import heapq
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
    >>> from basevar.caller import utils
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


class FileForQueueing(object):
    """
    """
    def __init__(self, the_file, line):
        """
        Store the file, and initialise the current value
        """
        self.the_file = the_file
        self.finishedReadingFile = False
        self.heap = []

        line = line
        cols = line.strip().split("\t")
        chrom = cols[0]

        # Where possible, convert chromosome names into
        # integers for sorting. If not possible, use
        # original names.
        try:
            chrom = int(chrom.upper().strip("CHR"))
        except Exception:
            pass

        pos = int(cols[1])
        heapq.heappush(self.heap, (chrom, pos, line))

        while not self.finishedReadingFile and len(self.heap) < 100:

            try:
                line = self.the_file.next()
                cols = line.strip().split("\t")
                chrom = cols[0]

                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
            except StopIteration:
                self.finishedReadingFile = True
                break

            heapq.heappush(self.heap, (chrom, pos, line))

        # Now take the top line
        self.chrom, self.pos, self.line = heapq.heappop(self.heap)

    def __cmp__(self, other):
        """
        Comparison function. Utilises the comparison function defined in
        the AlignedRead class.
        """
        return cmp(self.chrom, other.chrom) or cmp(self.pos, other.pos)

    def __del__(self):
        """
        Destructor
        """
        self.the_file.close()
        os.remove(self.the_file.name)

    def next(self):
        """
        Increment the iterator and yield the new value. Also, store the
        current value for use in the comparison function.
        """
        if not self.finishedReadingFile:

            try:
                line = self.the_file.next()
                #cols = line.strip().split('\t')
                cols = line.strip().split()
                chrom = cols[0]

                # Where possible, convert chromosome names into
                # integers for sorting. If not possible, use
                # original names.
                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
                heapq.heappush(self.heap, (chrom, pos, line))

            except StopIteration:
                self.finishedReadingFile = True

        if len(self.heap) != 0:
            # Now take the top line
            self.chrom, self.pos, self.line = heapq.heappop(self.heap)
        else:
            raise StopIteration


def merge_files(temp_file_names, final_file_name, is_detle_raw_file=False):
    """
    Merging output VCF/CVG files into a final file
    log.info("Merging output VCF/CVG file(s) into final file %s" %(final_file_name))
    """

    # Final output file
    if final_file_name == "-":
        output_file = sys.stdout
    else:
        output_file = Open(final_file_name, 'wb')
    the_heap = []

    # Initialise queue
    for index, file_name in enumerate(temp_file_names):
        the_file = Open(file_name, 'rb')

        for line in the_file:

            # End of this file
            if line[0] == "#":
                if index == 0:
                    output_file.write(line)
            else:
                the_file_for_queueing = FileForQueueing(the_file, line)
                heapq.heappush(the_heap, the_file_for_queueing)
                break

        # If there are no calls in the temp file, we still want to
        # remove it.
        else:
            the_file.close()

            if is_detle_raw_file:
                os.remove(file_name)

    # Merge-sort the output using a priority queue
    while len(the_heap) != 0:

        # Get file from heap in right order
        next_file = heapq.heappop(the_heap)
        output_file.write(next_file.line)

        # Put file back on heap
        try:
            next_file.next()
            heapq.heappush(the_heap, next_file)
        except StopIteration:
            continue

    # Close final output file
    if final_file_name != "-":
        output_file.close()

    #log.info("Finished merging %s file(s)"%final_file_name)
    return

