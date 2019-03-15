"""
Package for parsing bamfile
Author: Shujia Huang
Date : 2016-07-19 14:14:21
"""
import os
import sys
import time

import pysam

from . import utils


def open_align_files(bamfiles):
    ali_files_hd = []
    for f in bamfiles:

        try:
            bf = pysam.AlignmentFile(f)

        except ValueError:
            sys.stderr.write('[ERROR] Input file: %s is not BAM/CRAM.\n' % f)
            close_align_file(ali_files_hd)
            sys.exit(1)

        ali_files_hd.append(bf)

    return ali_files_hd


def close_align_file(ali_files_hd):
    for bf in ali_files_hd:
        bf.close()
    return


def create_batchfiles_for_regions(chrid, regions, batchcount, align_files, fa, mapq, outdir,
                                  sample_ids=None, is_smart_rerun=False):
    """
    ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
    ``fa``:
        # get sequence of chrid from reference fasta
        fa = self.ref_file_hd.fetch(chrid)
    """
    # store all the batch files
    batchfiles = []
    part_num = len(align_files) / batchcount
    if part_num * batchcount < len(align_files):
        part_num += 1

    tmp_region = []
    for p in regions:
        tmp_region.extend(p)

    tmp_region = sorted(tmp_region)
    bigstart, bigend = tmp_region[0], tmp_region[-1]

    m = 0
    for i in range(0, len(align_files), batchcount):
        # Create a batch of temp files for variant discovery
        start_time = time.time()

        m += 1
        part_file_name = "BaseVar.%s.%d_%d.batch.gz" % (".".join(map(str, [chrid, bigstart, bigend])),
                                                        m, part_num)
        part_file_name = os.path.join(outdir, part_file_name)  # Join Path could fix different OS

        # store the name of batchfiles into a list.
        batchfiles.append(part_file_name)
        if is_smart_rerun and os.path.isfile(part_file_name):
            # ``part_file_name`` is exists We don't have to create it again if setting `smartrerun`
            sys.stderr.write("[INFO] %s already exists, we don't have to create it again, "
                             "when you set `smartrerun` %s\n" % (part_file_name, time.asctime()))
            continue
        else:
            sys.stderr.write("[INFO] Creating batchfile %s at %s\n" % (part_file_name, time.asctime()))

        # One batch of alignment files
        sub_align_files = align_files[i:i + batchcount]
        batch_sample_ids = None
        if sample_ids:
            batch_sample_ids = sample_ids[i:i + batchcount]

        create_single_batchfile(chrid, bigstart, bigend, regions, sub_align_files, fa, mapq, part_file_name,
                                batch_sample_ids=batch_sample_ids)

        elapsed_time = time.time() - start_time
        sys.stderr.write("[INFO] Done for batchfile %s at %s, %d seconds elapsed\n"
                         "\n" % (part_file_name, time.asctime(), elapsed_time))

    return batchfiles


def create_single_batchfile(chrid, bigstart, bigend, regions, batch_align_files, fa, mapq,
                            out_batch_file, batch_sample_ids=None):
    # One batch of alignment files
    ali_files_hd = open_align_files(batch_align_files)

    # ``iter_tokes`` is a list of iterator for each sample's input file
    iter_tokes = []
    for j, bf in enumerate(ali_files_hd):
        try:
            # 0-base
            iter_tokes.append(bf.pileup(chrid, bigstart - 1, bigend))
        except ValueError:
            sys.stderr.write("# [WARMING] Empty region %s:%d-%d in %s" %
                             (chrid, bigstart - 1, bigend, batch_align_files[j]))
            iter_tokes.append("")

    with utils.Open(out_batch_file, "wb", isbgz=True) if out_batch_file.endswith(".gz") else \
            open(out_batch_file, "w") as OUT:

        OUT.write("##fileformat=BaseVarBatchFile_v1.0\n")
        if batch_sample_ids:
            OUT.write("##SampleIDs=%s\n" % ",".join(batch_sample_ids))

        OUT.write("%s\n" % "\t".join(
            ["#CHROM", "POS", "REF", "Depth(CoveredSample)", "MappingQuality", "Readbases",
             "ReadbasesQuality", "ReadPositionRank", "Strand"]))

        # Set iteration marker: 1->iterate; 0->Do not iterate or hit the end
        sample_info = [utils.fetch_next(it) for it in iter_tokes]

        n = 0
        for start, end in regions:

            sys.stderr.write('[INFO] Fetching information in region %s for %s, '
                             'at %s\n' % (chrid + ":" + str(start) + "-" + str(end),
                                          out_batch_file,
                                          time.asctime()))

            for position in range(start, end + 1):

                if n % 100000 == 0:
                    sys.stderr.write("[INFO] loading lines %d at position %s:%d\t%s\n" %
                                     (n + 1, chrid, position, time.asctime()))
                n += 1

                ref_base = fa[position - 1]
                # ignore 'N' bases in reference
                if ref_base.upper() not in ['A', 'C', 'G', 'T']:
                    continue

                (depth, sample_bases, sample_base_quals,
                 strands, mapqs, read_pos_rank) = fetch_base_by_position(
                    position - 1,  # position for pysam is 0-base
                    sample_info,
                    iter_tokes,
                    mapq,
                    fa  # Fa sequence for indel sequence
                )

                OUT.write("%s\n" % "\t".join([
                    chrid,
                    str(position),
                    ref_base,
                    str(depth),
                    ",".join(map(str, mapqs)),
                    ",".join(sample_bases),
                    ",".join(map(str, sample_base_quals)),
                    ",".join(map(str, read_pos_rank)),
                    ",".join(strands)
                ]))

    close_align_file(ali_files_hd)
    return


def fetch_base_by_position(position, sample_info, iter_tokes, mapq_thd, fa):
    """
    """
    base_quals = []
    bases = []
    strands = []
    mapqs = []
    read_pos_rank = []
    depth = 0

    for i, sample_pos_line in enumerate(sample_info):

        bs, qs, strand, mapq, rpr, sample_info[i] = seek_position(position, sample_pos_line,
                                                                  iter_tokes[i], mapq_thd, fa)

        bases.append(bs)
        base_quals.append(qs)
        strands.append(strand)
        mapqs.append(mapq)
        read_pos_rank.append(rpr)

        if qs:
            depth += 1

    return depth, bases, base_quals, strands, mapqs, read_pos_rank


def seek_position(target_pos, sample_pos_line, sample_iter, mapq_thd, fa):
    """Get mapping info for specific position.

    `fa`: Use for scanning indels
    """
    base, strand, qual, rpr, mapq = 'N', '.', 0, 0, 0  # Init
    if sample_pos_line:

        if sample_pos_line.pos < target_pos:

            pos = sample_pos_line.pos
            while pos < target_pos:

                sample_pos_line = utils.fetch_next(sample_iter)
                if sample_pos_line:
                    pos = sample_pos_line.pos
                else:
                    # hit the end of file, break the loop.
                    break

        # sample_pos_line may hit the end of file
        if sample_pos_line and sample_pos_line.pos == target_pos:
            base, strand, qual, mapq, rpr = first_base(sample_pos_line, sample_pos_line.pos, mapq_thd, fa)

    return base, qual, strand, mapq, rpr, sample_pos_line


def scan_indel(read, target_pos, fa):
    """Just scanning indel from pysam's pileups object.

    `target_pos`: 0-base
    `fa`: for fetch sequence from reference

    The cigar string order in the array is "MIDNSHP=X" followed by a
    field for the NM tag. If the NM tag is not present, this
    field will always be 0.

        +-----+--------------+-----+
        |M    |BAM_CMATCH    |0    |
        +-----+--------------+-----+
        |I    |BAM_CINS      |1    |
        +-----+--------------+-----+
        |D    |BAM_CDEL      |2    |
        +-----+--------------+-----+
        |N    |BAM_CREF_SKIP |3    |
        +-----+--------------+-----+
        |S    |BAM_CSOFT_CLIP|4    |
        +-----+--------------+-----+
        |H    |BAM_CHARD_CLIP|5    |
        +-----+--------------+-----+
        |P    |BAM_CPAD      |6    |
        +-----+--------------+-----+
        |=    |BAM_CEQUAL    |7    |
        +-----+--------------+-----+
        |X    |BAM_CDIFF     |8    |
        +-----+--------------+-----+
        |B    |BAM_CBACK     |9    |
        +-----+--------------+-----+
        |NM   |NM tag        |10   |
        +-----+--------------+-----+

    If no cigar string is present, empty arrays will be archived.
    """
    target_indx = 0
    delta = 0
    for i, (cigar_type, cigar_len) in enumerate(read.alignment.cigar):
        # If the cigar string is : 20M2I13M
        # then alignment.cigar is: [(0, 20), (1, 2), (0, 13)]
        # and alignment.blocks looks like: [(1121815, 1121835), (1121835, 1121848)].
        # But we should find the position of Insertion, which is the next one.

        if cigar_type in [3, 4, 5, 6]:  # 'SHPN'
            delta += 1
            continue

        # mapping
        if cigar_type == 0:
            _, map_end = read.alignment.blocks[target_indx]

            # map_end is 1-base and target_pos is 0-base
            if map_end == target_pos + 1:
                target_indx += 1  # +1 Get the index of indel in alignment.cigar
                break
            else:
                # +1
                target_indx += 1

    target_indx += delta
    cigar_type, cigar_len = read.alignment.cigar[target_indx]
    if cigar_type == 1:  # Insertion

        qpos = read.query_position + 1
        indel = '+' + read.alignment.query_sequence[qpos:qpos + cigar_len]
    elif cigar_type == 2:  # Deletion

        tpos = target_pos + 1
        indel = '-' + fa[tpos:tpos + cigar_len]
    else:
        # Must just be 1 or 2
        sys.stderr.write("[ERROR] Wrong Indel CIGAR number %s %s %s %s (at) %d in %s\n" %
                         (read.alignment.cigarstring, read.alignment.cigar,
                          read.alignment.blocks, read.alignment.cigar[target_indx],
                          target_indx, read.alignment))
        sys.exit(1)

    return indel if indel else 'N'


def first_base(sample_pos_line, position, mapq_thd, fa):
    """Just get first alignement base for each sample.
    """
    base, strand, qual, rpr, mapq = 'N', '.', 0, 0, 0  # Init
    for read in sample_pos_line.pileups:

        if read.alignment.mapq < mapq_thd:
            continue

        # skip read which mapping quality less then ``mapq_thd``
        strand = '-' if read.alignment.is_reverse else '+'
        mapq = read.alignment.mapq
        if read.indel:
            base = scan_indel(read, position, fa)
            break

        elif not read.is_del and not read.is_refskip:
            rpr = read.query_position + 1
            base = read.alignment.query_sequence[read.query_position]

            # Todo: much faster than `read.alignment.qqual[read.query_position]`
            #  and I don't know why.
            qual = read.alignment.query_qualities[read.query_position]

            # Just get base from first read and skip others,
            # no matter the first one it's indel or not.
            break

    return base, strand, qual, mapq, rpr
