"""
This is a Process module for BaseType

"""
import sys
import time

from pysam import TabixFile

from .basetype import BaseType
from .algorithm import strand_bias, ref_vs_alt_ranksumtest
from .utils import fetch_next, vcf_header_define, cvg_header_define, Open


def variants_discovery(chrid, fa, batchfiles, popgroup, cmm, cvg_file_handle, vcf_file_handle):
    """Function for variants discovery.
    """
    batch_files_hd = [Open(f, 'rb') for f in batchfiles]
    n, is_empty, eof = 0, True, False
    while True:
        # Loading ...
        # [CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand]
        sampleinfos = []
        for fh in batch_files_hd:
            line = fh.readline()
            if line:
                if line.startswith("#"):
                    continue

                sampleinfos.append(line.strip().split())
            else:
                sampleinfos.append(None)
                eof = True

        # hit the end of files
        if eof:
            is_error = True if any(sampleinfos) else False
            if is_error:
                sys.stderr.write(
                    "[ERROR] %s\n[ERROR]Error happen when 'variants_discovery', they don't have the same "
                    "positions in above files.\n" % "\n".join(batchfiles))
            break

        # Empty! Just have header information.
        if not sampleinfos:
            continue

        position = int(sampleinfos[0][1])
        if n % 10000 == 0:
            sys.stderr.write("[INFO] Have been loading %d lines when hit position %s:%d\t%s\n" %
                             (n if n > 0 else 1, chrid, position, time.asctime()))
        n += 1

        ref_base = fa[position - 1]

        (depth,
         sample_bases,
         sample_base_quals,
         strands,
         mapqs,
         read_pos_rank) = _fetch_baseinfo_by_position_from_batchfiles(chrid, position, ref_base, sampleinfos)

        # ignore if coverage=0
        if depth == 0:
            continue

        # Not empty
        is_empty = False

        # Calling varaints by Basetypes and output VCF and Coverage files.
        _basetypeprocess(chrid,
                         position,
                         ref_base,
                         sample_bases,
                         sample_base_quals,
                         mapqs,
                         strands,
                         read_pos_rank,
                         popgroup,
                         cmm,
                         cvg_file_handle,
                         vcf_file_handle)

    for fh in batch_files_hd:
        fh.close()

    return is_empty


def batchfile_variants_discovery(chrid, regions, fa, batchfiles, popgroup, cmm, cvg_file_handle, vcf_file_handle):
    """Function for variants discovery.
    """
    tmp_region = []
    for p in regions:
        tmp_region.extend(p)

    tmp_region = sorted(tmp_region)
    bigstart, bigend = tmp_region[0], tmp_region[-1]

    tbfs = [TabixFile(f) for f in batchfiles]
    batch_files_hd = [f.fetch(chrid, bigstart-1, bigend) for f in tbfs]
    n, is_empty, eof = 0, True, False
    while True:
        # Loading ...
        # [CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand]
        sampleinfos = []
        for fh in batch_files_hd:
            line = fetch_next(fh)
            if line:
                sampleinfos.append(line.strip().split())
            else:
                sampleinfos.append(None)
                eof = True

        # hit the end of files
        if eof:
            is_error = True if any(sampleinfos) else False
            if is_error:
                sys.stderr.write(
                    "[ERROR] %s\n[ERROR]Error happen when 'variants_discovery', they don't have the same "
                    "positions in above files.\n" % "\n".join(batchfiles))
            break

        # Empty! Just have header information.
        if not sampleinfos:
            continue

        chrid, position = sampleinfos[0][0], int(sampleinfos[0][1])
        if n % 10000 == 0:
            sys.stderr.write("[INFO] Have been loading %d lines when hit position %s:%d\t%s\n" %
                             (n if n > 0 else 1, chrid, position, time.asctime()))
        n += 1

        ref_base = fa[position - 1]

        (depth,
         sample_bases,
         sample_base_quals,
         strands,
         mapqs,
         read_pos_rank) = _fetch_baseinfo_by_position_from_batchfiles(chrid, position, ref_base, sampleinfos)

        # ignore if coverage=0
        if depth == 0:
            continue

        # Not empty
        is_empty = False

        # Calling varaints by Basetypes and output VCF and Coverage files.
        _basetypeprocess(chrid,
                         position,
                         ref_base,
                         sample_bases,
                         sample_base_quals,
                         mapqs,
                         strands,
                         read_pos_rank,
                         popgroup,
                         cmm,
                         cvg_file_handle,
                         vcf_file_handle)

    for fh in tbfs:
        fh.close()

    return is_empty


def _fetch_baseinfo_by_position_from_batchfiles(chrid, position, ref_base, infolines):
    sample_bases, sample_base_quals, mapqs, read_pos_rank, strands = [], [], [], [], []
    depth = 0
    for i, col in enumerate(infolines):
        # <CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand>
        if len(col) == 0:
            sys.stderr.write("[Error] %d lines happen to be empty in batchfiles!\n" % (i + 1))
            sys.exit(1)

        col[1], col[3] = map(int, [col[1], col[3]])
        if col[0] != chrid or col[1] != position or col[2] != ref_base:
            sys.stderr.write("[Error] %d lines, chromosome [%s and %s] or position [%d and %d] "
                             "or ref-base [%s and %s] in batchfiles not match with each other!\n" %
                             (i + 1, col[0], chrid, col[1], position, col[2], ref_base))
            sys.exit(1)

        depth += col[3]
        mapqs.append(col[4])
        sample_bases.append(col[5].upper())
        sample_base_quals.append(col[6])
        read_pos_rank.append(col[7])
        strands.append(col[8])

    # cat all the info together and create ...
    mapqs = map(int, ",".join(mapqs).split(","))
    sample_bases = ",".join(sample_bases).split(",")
    sample_base_quals = map(int, ",".join(sample_base_quals).split(","))
    read_pos_rank = map(int, ",".join(read_pos_rank).split(","))
    strands = ",".join(strands).split(",")

    return depth, sample_bases, sample_base_quals, strands, mapqs, read_pos_rank


def _basetypeprocess(chrid, position, ref_base, bases, base_quals, mapqs, strands, read_pos_rank,
                     popgroup, cmm, cvg_file_handle, vcf_file_handle):

    _out_cvg_file(chrid, position, ref_base, bases, strands, popgroup, cvg_file_handle, cmm=cmm)

    # Call variant if ``vcf_file_handle`` is not None
    if vcf_file_handle:

        bt = BaseType(ref_base.upper(), bases, base_quals, cmm=cmm)
        bt.lrt()

        if len(bt.alt_bases()) > 0:

            popgroup_bt = {}
            for group, index in popgroup.items():
                group_sample_bases, group_sample_base_quals = [], []
                for i in index:
                    group_sample_bases.append(bases[i])
                    group_sample_base_quals.append(base_quals[i])

                group_bt = BaseType(ref_base.upper(), group_sample_bases,
                                    group_sample_base_quals, cmm=cmm)

                basecombination = [ref_base.upper()] + bt.alt_bases()
                group_bt.lrt(basecombination)

                popgroup_bt[group] = group_bt

            _out_vcf_line(chrid,
                          position,
                          ref_base,
                          bases,
                          mapqs,
                          read_pos_rank,
                          base_quals,
                          strands,
                          bt,
                          popgroup_bt,
                          vcf_file_handle,
                          cmm=cmm)
    return


def _base_depth_and_indel(bases, cmm=None):
    # coverage info for each position
    base_depth = {b: 0 for b in cmm.BASE}
    indel_depth = {}

    for b in bases:

        if b == "N":
            continue

        if b in base_depth:
            # ignore all bases('*') which not match ``cmm.BASE``
            base_depth[b] += 1
        else:
            # Indel
            indel_depth[b] = indel_depth.get(b, 0) + 1

    indel_string = ','.join(
        [k + '|' + str(v) for k, v in indel_depth.items()]
    ) if indel_depth else "."

    return [base_depth, indel_string]


def output_header(fa_file_name, sample_ids, pop_group_sample_dict, out_cvg_handle, out_vcf_handle=None):

    info, group = [], []
    if pop_group_sample_dict:
        for g in pop_group_sample_dict.keys():
            g_id = g.split('_AF')[0]  # ignore '_AF'
            group.append(g_id)
            info.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="Allele frequency in the %s '
                        'populations calculated base on LRT, in the range (0,1)">' % (g_id, g_id))

    if out_vcf_handle:
        vcf_header = vcf_header_define(fa_file_name, info="\n".join(info), samples=sample_ids)
        out_vcf_handle.write("%s\n" % "\n".join(vcf_header))

    out_cvg_handle.write('%s\n' % "\n".join(cvg_header_define(group)))

    return


def _out_cvg_file(chrid, position, ref_base, bases, strands, popgroup, out_file_handle, cmm=None):
    """output coverage information into `out_file_handle`"""

    # coverage info for each position
    base_depth, indel_string = _base_depth_and_indel(bases, cmm=cmm)

    # base depth and indels for each subgroup
    group_cvg = {}
    for group, index in popgroup.items():

        group_sample_bases = []
        for i in index:
            group_sample_bases.append(bases[i])

        bd, ind = _base_depth_and_indel(group_sample_bases, cmm=cmm)
        group_cvg[group] = [bd, ind]

    fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = 0, -1, 0, 0, 0, 0
    if sum(base_depth.values()) > 0:
        base_sorted = sorted(base_depth.items(),
                             key=lambda x: x[1],
                             reverse=True)

        b1, b2 = base_sorted[0][0], base_sorted[1][0]
        fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
            ref_base.upper(),
            [b1 if b1 != ref_base.upper() else b2],
            bases,
            strands
        )

    if sum(base_depth.values()):

        group_info = []
        if group_cvg:
            for k in popgroup.keys():
                depth, indel = group_cvg[k]

                indel = [indel] if indel != "." else []
                s = ':'.join(map(str, [depth[b] for b in cmm.BASE]) + indel)
                group_info.append(s)

        out_file_handle.write(
            '\t'.join(
                [chrid, str(position), ref_base, str(sum(base_depth.values()))] +
                [str(base_depth[b]) for b in cmm.BASE] +
                [indel_string] +
                [str(fs), str(sor), ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev]))] + group_info
            ) + '\n'
        )

    return


def _out_vcf_line(chrid, position, ref_base, bases, mapqs, read_pos_rank, sample_base_qual,
                  strands, bt, pop_group_bt, out_file_handle, cmm=None):
    """output vcf lines into `out_file_handle`"""

    alt_gt = {b: './' + str(k + 1) for k, b in enumerate(bt.alt_bases())}
    samples = []

    for k, b in enumerate(bases):

        # For sample FORMAT
        if b != 'N' and b[0] not in ['-', '+']:
            # For the base which not in bt.alt_bases()
            if b not in alt_gt:
                alt_gt[b] = './.'

            gt = '0/.' if b == ref_base.upper() else alt_gt[b]

            samples.append(gt + ':' + b + ':' + strands[k] + ':' +
                           str(round(bt.qual_pvalue[k], 6)))
        else:
            samples.append('./.')  # 'N' base or indel

    # Rank Sum Test for mapping qualities of REF versus ALT reads
    mq_rank_sum = ref_vs_alt_ranksumtest(ref_base.upper(), bt.alt_bases(),
                                         zip(bases, mapqs))

    # Rank Sum Test for variant appear position among read of REF versus ALT
    read_pos_rank_sum = ref_vs_alt_ranksumtest(ref_base.upper(), bt.alt_bases(),
                                               zip(bases, read_pos_rank))

    # Rank Sum Test for base quality of REF versus ALT
    base_q_rank_sum = ref_vs_alt_ranksumtest(ref_base.upper(), bt.alt_bases(),
                                             zip(bases, sample_base_qual))

    # Variant call confidence normalized by depth of sample reads
    # supporting a variant.
    ad_sum = sum([bt.depth[b] for b in bt.alt_bases()])
    qd = round(float(bt.var_qual() / ad_sum), 3)

    # Strand bias by fisher exact test and Strand bias estimated by the
    # Symmetric Odds Ratio test
    fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
        ref_base.upper(), bt.alt_bases(), bases, strands)

    # base=>[CAF, allele depth], CAF = Allele frequency by read count
    caf = {b: ['%f' % round(bt.depth[b] / float(bt.total_depth), 6),
               bt.depth[b]] for b in bt.alt_bases()}

    info = {'CM_DP': str(int(bt.total_depth)),
            'CM_AC': ','.join(map(str, [caf[b][1] for b in bt.alt_bases()])),
            'CM_AF': ','.join(map(str, [bt.af_by_lrt[b] for b in bt.alt_bases()])),
            'CM_CAF': ','.join(map(str, [caf[b][0] for b in bt.alt_bases()])),
            'MQRankSum': str(mq_rank_sum),
            'ReadPosRankSum': str(read_pos_rank_sum),
            'BaseQRankSum': str(base_q_rank_sum),
            'QD': str(qd),
            'SOR': str(sor),
            'FS': str(fs),
            'SB_REF': str(ref_fwd) + ',' + str(ref_rev),
            'SB_ALT': str(alt_fwd) + ',' + str(alt_rev)}

    if pop_group_bt:

        for group, g_bt in pop_group_bt.items():
            af = ','.join(map(str, [g_bt.af_by_lrt[b] if b in g_bt.af_by_lrt else 0
                                    for b in bt.alt_bases()]))
            info[group] = af

    out_file_handle.write('\t'.join([chrid, str(position), '.', ref_base,
                                     ','.join(bt.alt_bases()), str(bt.var_qual()),
                                     '.' if bt.var_qual() > cmm.QUAL_THRESHOLD else 'LowQual',
                                     ';'.join([k + '=' + v for k, v in sorted(
                                         info.items(), key=lambda x: x[0])]),
                                     'GT:AB:SO:BP'] + samples) + '\n')
    return
