import sys
import os
import heapq
import time

import cProfile
import pstats

from basevar.log import logger

from basevar.io.BGZF.tabix import tabix_index
from basevar.io.fasta cimport FastaFile
from basevar.io.openfile import Open, FileForQueueing

def do_cprofile(filename, is_do_profiling=False, stdout=False):
    """
    Decorator for function profiling.
    """
    def wrapper(func):
        def profiled_func(*args, **kwargs):
            # Flag for do profiling or not.
            if is_do_profiling:
                profile = cProfile.Profile()
                profile.enable()
                result = func(*args, **kwargs)
                profile.disable()

                # Sort stat by internal time.
                sortby = "tottime"
                ps = pstats.Stats(profile).sort_stats(sortby)
                ps.dump_stats(filename)

                # Output profile to STDOUT
                if stdout:
                    ps.strip_dirs().sort_stats("time").print_stats()
            else:
                result = func(*args, **kwargs)

            return result
        return profiled_func

    return wrapper

cdef long int c_min(long int x, long int y):
    if x < y:
        return x
    else:
        return y

cdef long int c_max(long int x, long int y):
    if x < y:
        return y
    else:
        return x

cdef float set_minaf(int sample_size):
    """setting the resolution of MAF for basetype"""
    return min(100.0 / sample_size, 0.001)

def safe_remove(bytes fname):
    """Remove a file if it exist"""
    if not fname:
        return False

    if os.path.exists(fname):
        os.remove(fname)

    return True

def safe_makedir(bytes dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname

    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise

            num_tries += 1
            time.sleep(2)

    return dname

def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    try:
        return fname and os.path.exists(fname) and os.path.getsize(fname) > 0
    except OSError:
        return False

def get_last_modification_file(bytes dirname):
    """Find the last modification file in a directory and return it."""

    file_name_list = os.listdir(dirname)
    file_name_list.sort(key=lambda fn: os.path.getmtime(os.path.join(dirname, fn)))

    if file_name_list:

        # return the last modification file path of dirname
        return os.path.join(dirname, file_name_list[-1])
    else:
        # dirname is empty
        return ""

def vcf_header_define(ref_file_path, info=None, samples=None):
    """define header for VCF"""

    if not samples:
        samples = []

    fa = FastaFile(ref_file_path, ref_file_path + ".fai")
    fa_name = os.path.basename(fa.filename)
    contigs = ["##contig=<ID=%s,length=%d,assembly=%s>" % (c, s, fa_name) for c, s in zip(fa.refnames, fa.lengths)]
    header = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=LowQual,Description="Low quality (QUAL < 60)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AB,Number=1,Type=String,Description="Allele Base">',
        '##FORMAT=<ID=SO,Number=1,Type=String,Description="Strand orientation of the mapping base. Marked as + or -">',
        '##FORMAT=<ID=BP,Number=1,Type=String,Description="Base Probability which calculate by base quality">',

        '##INFO=<ID=CM_AF,Number=A,Type=Float,Description="An ordered, comma delimited list of allele frequencies base on LRT algorithm">',
        '##INFO=<ID=CM_CAF,Number=A,Type=Float,Description="An ordered, comma delimited list of allele frequencies just base on read count">',
        '##INFO=<ID=CM_AC,Number=A,Type=Integer,Description="An ordered, comma delimited allele depth in CMDB">',
        '##INFO=<ID=CM_DP,Number=A,Type=Integer,Description="Total Depth in CMDB">',
        '##INFO=<ID=SB_REF,Number=A,Type=Integer,Description="Read number support REF: Forward,Reverse">',
        '##INFO=<ID=SB_ALT,Number=A,Type=Integer,Description="Read number support ALT: Forward,Reverse">',
        '##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">',
        '##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Phred-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">',
        '##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">',
        '##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Phred-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">',
        '##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Phred-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">',
        '##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence Quality by Depth">',
        '\n'.join([info] + contigs if info else contigs),
        '##reference=file://{}'.format(os.path.realpath(fa.filename)),
        '\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'] + samples)
    ]

    fa.close()

    return header

def cvg_header_define(group_info):
    """define header for coverage file"""
    h = '\t'.join(['#CHROM', 'POS', 'REF', 'Depth', 'A', 'C', 'G', 'T'] +
                  ['Indels', 'FS', 'SOR', 'Strand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)'])

    header = [
        "##fileformat=CVGv1.0",
        "##Group information is the depth of A:C:G:T:Indel",
        "%s\t%s" % (h, '\t'.join(group_info)) if group_info else h
    ]

    return header

cdef list generate_regions_by_process_num(list regions, int process_num, bint convert_to_2d):
    """create regions for each process"""

    regions_for_each_process = [[] for _ in range(process_num)]

    # calculate region size for each processes
    total_regions_size = sum([e - s + 1 for _, s, e in regions])
    delta = int(total_regions_size / process_num) + 1 if total_regions_size % process_num else \
        int(total_regions_size / process_num)

    # reset
    total_regions_size = 0
    for chrid, start, end in regions:

        sub_region_size = end - start + 1
        # find index of regions_for_each_process by the total_regions_size
        i = int(total_regions_size / delta)
        total_regions_size += sub_region_size

        pre_size = 0
        if len(regions_for_each_process[i]):
            pre_size = sum([e - s + 1 for _, s, e in regions_for_each_process[i]])

        if sub_region_size + pre_size > delta:
            d = delta - pre_size
            regions_for_each_process[i].append([chrid, start, start + d - 1])

            i += 1
            for k, pos in enumerate(range(start + d - 1, end, delta)):
                s = pos + 1 if pos + 1 < end else end
                e = pos + delta if pos + delta < end else end

                regions_for_each_process[i + k].append([chrid, s, e])

        else:
            regions_for_each_process[i].append([chrid, start, end])

    cdef list region_array = []
    if convert_to_2d:
        for a1 in regions_for_each_process:
            for a2 in a1:
                region_array.append(a2)

        return region_array
    else:
        return regions_for_each_process

def fetch_next(iter_fh):
    """
    re-define the next funtion of fetching info from pysam
    prevent throunghing the 'StopIteration' exception.
    """

    if iter_fh == '':
        return ''

    try:
        line = iter_fh.next()
    except StopIteration:
        line = ''

    return line

def load_file_list(in_file):
    with Open(in_file, 'rb') as fh:
        files = [r.strip().split()[0] for r in fh if r[0] != '#']

    return files

def load_target_position(bytes referencefile, bytes posfile, bytes region_info):
    # Loading positions
    _sites = get_position_list(posfile) if posfile else {}

    fa = FastaFile(referencefile, referencefile + ".fai")
    if len(region_info):

        regions = []
        if os.path.isfile(region_info):
            regions = get_region_fromfile(region_info)

        else:
            for r in region_info.split(','):

                if ':' in r:
                    chr_id, reg = r.strip().split(':')
                    start, end = map(int, reg.split('-'))
                    regions.append([chr_id, start, end])
                else:
                    # Just chromosome id
                    ci = r.strip()
                    regions.append([ci, 1, fa.get_reference_length(ci)])

        for chrid, start, end in regions:

            if chrid not in _sites:
                _sites[chrid] = []

            _sites[chrid].append([start, end])

    # sort and merge the regions
    # [[chrid1, start1, end1], [chrid2, start2, end2], ...]
    regions = []
    for chrid, v in sorted(_sites.items(), key=lambda x: x[0]):
        for start, end in merge_region(v):
            regions.append([chrid, start, end])

    # load all the genome if no position or regions provide
    if not regions:
        sys.stderr.write('[WARNINGS] Program will load all the genome. This will take a long long time.\n')
        regions = [[ci, 1, fa.get_reference_length(ci)]
                   for ci in fa.refnames]

    fa.close()

    return regions

def get_position_list(in_site_file):
    sites = {}
    with Open(in_site_file, 'rb') as f:
        for r in f:
            tok = r.strip().split()
            if tok[0] not in sites:
                sites[tok[0]] = []

            if len(tok) < 3:
                sites[tok[0]].append([int(tok[1]), int(tok[1])])
            else:
                sites[tok[0]].append([int(tok[1]), int(tok[2])])

    return sites

def get_region_fromfile(in_region_file):
    regions = []
    with Open(in_region_file, 'rb') as f:
        for r in f:
            chr_id, reg = r.strip().split(':')
            start, end = map(int, reg.split('-'))
            regions.append([chr_id, start, end])

    return regions

def regions2dict(regions):
    """
    Convert a 2d array into a dict
    """
    reg_dict = {}
    # store the region into a dict
    for chrid, start, end in regions:

        if chrid not in reg_dict:
            reg_dict[chrid] = []

        reg_dict[chrid].append([start, end])

    return reg_dict

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
    position_region.sort(key=lambda x: x[0])

    m_region = []
    prepos, start, end = '', '', ''
    flag = False
    for s, e in position_region:
        if s > e:
            sys.stderr.write(('[ERROR]Your region start > end. It is not allow '
                              'when call Merge function\n'))
            sys.exit(1)

        if prepos == '':
            # The light is on => Get the region!
            if flag:
                m_region.append([start, end])

            start, end = s, e
            flag = True

        else:
            if prepos > s:
                sys.stderr.write('[ERROR] Array is un-sorted.\n')
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
    sorted_bc = sorted(bc.items(), key=lambda x: x[1], reverse=True)

    if len(sorted_bc):
        mi, mj = sorted_bc[-1][0], sorted_bc[0][0]
    else:
        mi, mj = 'N', 'N'

    # minor and major
    return mi, mj, sorted_bc

cdef void fast_merge_files(list temp_file_names, basestring final_file_name, bint is_del_raw_file):
    """Merge file which is already in order.
    We don't have to sort anything, just cat them together!
    """
    # Final output file
    if final_file_name == "-":
        output_file = sys.stdout
    else:
        output_file = Open(final_file_name, 'wb', isbgz=True if final_file_name.endswith(".gz") else False)

    cdef int index = 0
    cdef int total_file_num = len(temp_file_names)
    cdef basestring file_name
    cdef dict check_position_in_order = {}
    for index, file_name in enumerate(temp_file_names):

        the_file = Open(file_name, 'rb')
        for line in the_file:

            # End of this file
            if line[0] == "#":
                if index == 0:
                    # Keep the header of first file.
                    output_file.write(line.strip() + "\n")
            else:

                chrom, pos = line.strip().split()[:2]
                pos = int(pos)
                if chrom not in check_position_in_order:
                    check_position_in_order[chrom] = pos

                if pos < check_position_in_order[chrom]:
                    sys.stderr.write(("[ERROR] Previous position (%d) > the following position(%d). "
                                      "This is not allow when call Merge function. You must keep increasing "
                                      "order in of all your input files.\n"))
                    sys.exit(1)

                output_file.write(line.strip() + "\n")

        # If there are no calls in the temp file, we still want to
        # remove it.
        else:
            the_file.close()
            if is_del_raw_file:
                os.remove(file_name)

        logger.info("Fast merge process %d/%d done." % (index+1, total_file_num))

    return

def merge_files(temp_file_names, final_file_name, is_del_raw_file=False):
    """
    Merging output VCF/CVG files into a final big one
    """
    # Final output file
    if final_file_name == "-":
        output_file = sys.stdout
    else:
        output_file = Open(final_file_name, 'wb', isbgz=True if final_file_name.endswith(".gz") else False)

    the_heap = []

    # Initialise queue
    for index, file_name in enumerate(temp_file_names):
        the_file = Open(file_name, 'rb')

        for line in the_file:

            # End of this file
            if line[0] == "#":
                if index == 0:
                    output_file.write(line.strip() + "\n")
            else:
                the_file_for_queueing = FileForQueueing(the_file, line.strip() + "\n", is_del_raw_file=is_del_raw_file)
                heapq.heappush(the_heap, the_file_for_queueing)
                break

        # If there are no calls in the temp file, we still want to
        # remove it.
        else:
            the_file.close()

            if is_del_raw_file:
                os.remove(file_name)

    # Merge-sort the output using a priority queue
    while len(the_heap) != 0:

        # Get file from heap in right order
        next_file = heapq.heappop(the_heap)
        output_file.write(next_file.line.strip() + "\n")

        # Put file back on heap
        try:
            next_file.next()
            heapq.heappush(the_heap, next_file)
        except StopIteration:
            continue

    # Close final output file
    if final_file_name != "-":
        output_file.close()

    return

def merge_batch_files(temp_file_names, final_file_name, output_isbgz=False, is_del_raw_file=False, justbase=False):
    """
    Merging output batch files into a final big one.

    ``temp_file_names``: must contain the same positions but different samples per file
    """
    # Final output file
    if final_file_name == "-":
        output_file = sys.stdout
    else:
        output_file = Open(final_file_name, 'wb', isbgz=output_isbgz)

    batch_files_hd = [Open(f, 'rb') for f in temp_file_names]
    eof = False
    while not eof:

        # [CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand]
        sampleinfos = []
        for fh in batch_files_hd:

            line = fh.readline()
            if line:
                sampleinfos.append(line.strip())
            else:
                sampleinfos.append(None)
                eof = True

        # hit the end of files
        if eof:
            is_error = True if any(sampleinfos) else False
            if is_error:
                sys.stderr.write("[ERROR] %s\n[ERROR]Error happen when 'merge_batch_files', they don't have the same "
                                 "positions in above files.\n" % "\n".join(temp_file_names))
            continue

        if sampleinfos[0].startswith("#"):
            # Header line
            if sampleinfos[0].startswith("##SampleIDs="):
                sampleid_info = [s.split("=")[1] for s in sampleinfos]
                output_file.write("##SampleIDs=%s\n" % ",".join(sampleid_info))
            else:
                output_file.write("%s\n" % sampleinfos[0])

        else:
            chrid, position, ref_base = sampleinfos[0].split()[0:3]

            read_bases, read_base_quals, mapqs, read_pos_rank, strands = [], [], [], [], []
            depth = 0
            for i, line in enumerate(sampleinfos):
                # <CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand>
                col = line.split()
                if len(col) == 0:
                    sys.stderr.write("[Error] %d lines happen to be empty in batchfiles!\n" % (i + 1))
                    sys.exit(1)

                col[1], col[3] = col[1], col[3]
                if col[0] != chrid or col[1] != position or col[2] != ref_base:
                    sys.stderr.write("[Error] Error happen when 'merge_batch_files' in %d lines, chromosome "
                                     "[%s and %s] or position [%s and %s] or ref-base [%s and %s] in batchfiles "
                                     "not match with each other!\n" % (i + 1, col[0], chrid, col[1],
                                                                       position, col[2], ref_base))
                    sys.exit(1)

                depth += int(col[3])

                if not justbase:
                    mapqs.append(col[4])
                    read_bases.append(col[5])
                    read_base_quals.append(col[6])
                    read_pos_rank.append(col[7])
                    strands.append(col[8])
                else:
                    # The 5th column is readbase!
                    read_bases.append(col[4])

            # cat all the info together and create ...
            depth = str(depth)

            if not justbase:
                mapqs = ",".join(mapqs)
                read_bases = ",".join(read_bases)
                read_base_quals = ",".join(read_base_quals)
                read_pos_rank = ",".join(read_pos_rank)
                strands = ",".join(strands)

                # [CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand]
                output_file.write("%s\n" % "\t".join([chrid, position, ref_base, depth, mapqs, read_bases,
                                                      read_base_quals, read_pos_rank, strands]))
            else:
                # [CHROM POS REF Depth Readbases]
                output_file.write("%s\n" % "\t".join([chrid, position, ref_base, depth, ",".join(read_bases)]))

    else:

        for i, fh in enumerate(batch_files_hd):
            fh.close()

            if is_del_raw_file:
                os.remove(temp_file_names[i])

    # Close final output file
    if final_file_name != "-":
        output_file.close()

    return

def load_popgroup_info(samples, in_popgroup_file):
    """loading population group"""

    tmpdict = {}
    line_num = 0
    with Open(in_popgroup_file, 'rb') as f:

        # Just two columns: sample_id and group_id
        for line in f:
            line_num += 1

            try:
                sample_id, group_id = line.strip().split()[0:2]
            except ValueError:
                sys.stderr.write('[ERROR] Format error in `in_popgroup_file` it '
                                 'may not contain two columns happen in: %d '
                                 'lines in file "%s" \n' % (line_num, in_popgroup_file))
                sys.exit(1)

            tmpdict[sample_id] = group_id + '_AF'

    # group_id => [a list samples_index]
    cdef dict popgroup = {}
    cdef int i
    cdef basestring s
    # keep the sample's order
    for i, s in enumerate(samples):

        if s in tmpdict:
            if tmpdict[s] not in popgroup:
                popgroup[tmpdict[s]] = []

            # record different index of different groups
            popgroup[tmpdict[s]].append(i)

    return popgroup

def output_cvg_and_vcf(sub_cvg_files, sub_vcf_files, outcvg, outvcf=None):
    """CVG file and VCF file could use the same tabix strategy."""
    for out_final_file, sub_file_list in zip([outcvg, outvcf], [sub_cvg_files, sub_vcf_files]):

        if out_final_file:
            output_file(sub_file_list, out_final_file, del_raw_file=True)

    return

def output_file(sub_files, out_file_name, del_raw_file=False):
    if out_file_name.endswith(".gz"):
        merge_files(sub_files, out_file_name, is_del_raw_file=del_raw_file)

        # Column indices are 0-based. Note: this is different from the tabix command line
        # utility where column indices start at 1.
        tabix_index(out_file_name, force=True, seq_col=0, start_col=1, end_col=1)
    else:
        merge_files(sub_files, out_file_name, is_del_raw_file=del_raw_file)

    return
