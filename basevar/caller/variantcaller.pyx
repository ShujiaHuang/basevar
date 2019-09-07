# cython: profile=True
"""This is a Process module for BaseType
"""
import sys

from basevar.log import logger

from basevar.io.openfile import Open
from basevar.utils import CommonParameter, vcf_header_define, cvg_header_define
from basevar.caller.algorithm cimport strand_bias
from basevar.caller.algorithm cimport ref_vs_alt_ranksumtest

from basevar.caller.basetype cimport BaseType
from basevar.caller.batch cimport BatchInfo


cdef bint variants_discovery(bytes chrid, list batchfiles, dict popgroup, float min_af,
                             int batch_count, cvg_file_handle, vcf_file_handle):
    """Function for variants discovery
    """
    cdef list sampleinfos = []
    cdef list batch_files_hd = [Open(f, 'rb') for f in batchfiles]
    cdef bint is_empty = True
    cdef bint eof = False
    cdef bint is_error = False

    cdef int total_sample_num = batch_count * len(batchfiles)
    cdef BatchInfo batchinfo = BatchInfo(chrid, total_sample_num)

    cdef int n = 0
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
                logger.warning(
                    "%s\n[ERROR]Error happen when 'variants_discovery', they don't have the same "
                    "positions in above files." % "\n".join(batchfiles))
            break

        # Empty!! Just get header information.
        if not sampleinfos:
            continue

        # reset position and ref_base
        batchinfo.position = int(sampleinfos[0][1])
        batchinfo.ref_base = sampleinfos[0][2]

        if n % 10000 == 0:
            logger.info("Have been loading %d lines when hit position %s:%d" %
                        (n if n > 0 else 1, chrid, batchinfo.position))
        n += 1

        # data in ``batchinfo`` will been updated automatically in this funcion
        _fetch_baseinfo_by_position_from_batchfiles(sampleinfos, batch_count, batchinfo)

        # ignore if coverage=0
        if batchinfo.depth == 0:
            continue

        # Not empty
        is_empty = False

        # Calling varaints position one by one and output files.
        _basetypeprocess(batchinfo,
                         popgroup,
                         min_af,
                         cvg_file_handle,
                         vcf_file_handle)

    for fh in batch_files_hd:
        fh.close()

    return is_empty


cdef void _fetch_baseinfo_by_position_from_batchfiles(list infolines, int batch_count, BatchInfo batchinfo):

    # reset depth
    batchinfo.depth = 0

    cdef char *c_t4
    cdef char *c_t5
    cdef char *c_t6
    cdef char *c_t7
    cdef char *c_t8

    cdef int n = 0, index = 0
    for i, col in enumerate(infolines):
        # <CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand>
        if len(col) == 0:
            logger.error(" %d lines happen to be empty in batchfiles!" % (i + 1))
            sys.exit(1)

        col[1], col[3] = map(int, [col[1], col[3]])
        if col[0] != batchinfo.chrid or col[1] != batchinfo.position or col[2] != batchinfo.ref_base:
            logger.error("%d lines, chromosome [%s and %s] or position [%d and %d] "
                         "or ref-base [%s and %s] in batchfiles not match with each other!\n" %
                         (i + 1, col[0], batchinfo.chrid, col[1], batchinfo.position, col[2], batchinfo.ref_base))
            sys.exit(1)

        batchinfo.depth += col[3]
        if col[3] > 0:

            c_t4, c_t5, c_t6, c_t7, c_t8 = col[4:9]
            for n in range(batch_count):

                # if catch segmentation fault then the problem would probably be here!
                batchinfo.mapqs[index] = atoi(strsep(&c_t4, ","))
                batchinfo.sample_bases[index] = strsep(&c_t5, ",")  # must all be all upper charater in batchfile!
                batchinfo.sample_base_quals[index] = atoi(strsep(&c_t6, ","))
                batchinfo.read_pos_rank[index] = atoi(strsep(&c_t7, ","))
                batchinfo.strands[index] = strsep(&c_t8, ",")[0] # It's char not string

                # move to the next
                index += 1
        else:
            for n in range(batch_count):

                batchinfo.mapqs[index] = 0
                batchinfo.sample_bases[index] = "N"
                batchinfo.sample_base_quals[index] = 0
                batchinfo.read_pos_rank[index] = 0
                batchinfo.strands[index] = "."

                # move to the next
                index += 1

    return


cdef void _basetypeprocess(BatchInfo batchinfo, dict popgroup, float min_af, cvg_file_handle, vcf_file_handle):

    _out_cvg_file(batchinfo, popgroup, cvg_file_handle)

    cdef dict popgroup_bt = {}
    cdef bint is_variant = True
    cdef BaseType bt, group_bt

    cdef list bases = []
    cdef list strands = []
    cdef list base_quals = []
    cdef int i = 0
    if vcf_file_handle:

        # change later
        for i in range(batchinfo.size):
            bases.append(batchinfo.sample_bases[i])
            base_quals.append(batchinfo.sample_base_quals[i])
            strands.append(chr(batchinfo.strands[i]))

        bt = BaseType(batchinfo.ref_base.upper(), bases, base_quals, min_af)
        is_variant = bt.lrt(None)  # do not need to set specific_base_combination

        if is_variant:

            popgroup_bt = {}
            for group, index in popgroup.items():
                group_sample_bases, group_sample_base_quals = [], []
                for i in index:
                    group_sample_bases.append(bases[i])
                    group_sample_base_quals.append(base_quals[i])

                group_bt = BaseType(batchinfo.ref_base.upper(), group_sample_bases, group_sample_base_quals, min_af)

                group_bt.lrt([batchinfo.ref_base.upper()] + bt.alt_bases)
                popgroup_bt[group] = group_bt


            mapqs, read_pos_rank = [], []
            for i in range(batchinfo.size):
                mapqs.append(batchinfo.mapqs[i])
                read_pos_rank.append(batchinfo.read_pos_rank[i])

            _out_vcf_line(batchinfo, bt, popgroup_bt, vcf_file_handle)
    return


cdef list _base_depth_and_indel(char **bases, int size):
    # coverage info for each position
    cdef dict base_depth = {b: 0 for b in CommonParameter.BASE}
    cdef dict indel_depth = {}

    cdef int i = 0
    # for b in bases:
    for i in range(size):

        # The size of `base[i]` would be 1 except indel!
        if bases[i][0] == 'N':
            continue

        if bases[i] in base_depth:
            # ignore all bases('*') which not match ``cmm.BASE``
            base_depth[bases[i]] += 1
        else:
            # Indel
            indel_depth[bases[i]] = indel_depth.get(bases[i], 0) + 1

    cdef bytes indels = bytes(','.join(
        [k + '|' + str(v) for k, v in indel_depth.items()]
    ) if indel_depth else ".")

    return [base_depth, indels]

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

cdef void _out_cvg_file(BatchInfo batchinfo, dict popgroup, out_file_handle):
    """output coverage information into `out_file_handle`"""

    # coverage info for each position
    cdef dict base_depth
    cdef bytes indels
    base_depth, indels = _base_depth_and_indel(batchinfo.sample_bases, batchinfo.size)

    # base depth and indels for each subgroup
    cdef char **group_sample_bases

    cdef dict group_cvg = {}
    cdef bytes group
    cdef list index
    cdef int index_size

    cdef int i = 0
    cdef dict sub_bd
    cdef bytes sub_inds
    for group, index in popgroup.items():

        index_size = len(index)
        group_sample_bases = <char**>(calloc(index_size, sizeof(char*)))
        assert group_sample_bases != NULL, "Could not allocate memory for ``group_sample_bases`` in _out_cvg_file."

        for i in range(index_size):
            group_sample_bases[i] = batchinfo.sample_bases[index[i]]

        sub_bd, sub_inds = _base_depth_and_indel(group_sample_bases, index_size)
        group_cvg[group] = [sub_bd, sub_inds]

        free(group_sample_bases)

    cdef double fs, sor
    cdef int ref_fwd, ref_rev, alt_fwd, alt_rev
    fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = 0, -1, 0, 0, 0, 0

    cdef bytes ref_base = batchinfo.ref_base
    cdef bytes b1, b2
    if sum(base_depth.values()) > 0:
        base_sorted = sorted(base_depth.items(), key=lambda x: x[1], reverse=True)
        b1, b2 = base_sorted[0][0], base_sorted[1][0]
        fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
            ref_base.upper(),  # reference
            [b1 if b1 != ref_base.upper() else b2],  # alt-allele
            batchinfo.sample_bases,
            batchinfo.strands,
            batchinfo.size
        )

    cdef list group_info
    if sum(base_depth.values()):

        group_info = []
        if group_cvg:
            for k in group_cvg.keys():
                depth, indel = group_cvg[k]
                indel = [indel] if indel != "." else []
                s = ':'.join(map(str, [depth[b] for b in CommonParameter.BASE]) + indel)
                group_info.append(s)

        out_file_handle.write(
            '\t'.join(
                [batchinfo.chrid, str(batchinfo.position), ref_base, str(sum(base_depth.values()))] +
                [str(base_depth[b]) for b in CommonParameter.BASE] +
                [indels] +
                [str("%.3f" % fs), str("%.3f" % sor), ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev]))] +
                group_info
            ) + '\n'
        )

    return

cdef void _out_vcf_line(BatchInfo batchinfo, BaseType bt, dict pop_group_bt, out_file_handle):
    """output vcf lines into `out_file_handle`"""

    cdef dict alt_gt = {b: './' + str(k + 1) for k, b in enumerate(bt.alt_bases)}
    cdef list samples = []
    cdef int k
    cdef char *b
    # for k, b in enumerate(bases):
    for k in range(batchinfo.size):

        b = batchinfo.sample_bases[k]
        # For sample FORMAT
        if b[0] not in ['N', '-', '+']:
            # For the base which not in bt.alt_bases()
            if b not in alt_gt:
                alt_gt[b] = './.'

            gt = '0/.' if b == batchinfo.ref_base.upper() else alt_gt[b]

            samples.append(gt + ':' + b + ':' + chr(batchinfo.strands[k]) + ':' +
                           str(round(bt.qual_pvalue[k], 6)))
        else:
            samples.append('./.')  # 'N' base or indel

    # Rank Sum Test for mapping qualities of REF versus ALT reads
    mq_rank_sum = ref_vs_alt_ranksumtest(batchinfo.ref_base.upper(), bt.alt_bases, batchinfo.sample_bases,
                                         batchinfo.mapqs, batchinfo.size)

    # Rank Sum Test for variant appear position among read of REF versus ALT
    read_pos_rank_sum = ref_vs_alt_ranksumtest(batchinfo.ref_base.upper(), bt.alt_bases, batchinfo.sample_bases,
                                               batchinfo.read_pos_rank, batchinfo.size)

    # Rank Sum Test for base quality of REF versus ALT
    base_q_rank_sum = ref_vs_alt_ranksumtest(batchinfo.ref_base.upper(), bt.alt_bases, batchinfo.sample_bases,
                                             batchinfo.sample_base_quals, batchinfo.size)

    # Variant call confidence normalized by depth of sample reads
    # supporting a variant.
    ad_sum = sum([bt.depth[b] for b in bt.alt_bases])
    qd = round(float(bt.var_qual/ad_sum), 3)

    # Strand bias by fisher exact test and Strand bias estimated by the
    # Symmetric Odds Ratio test
    fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
        batchinfo.ref_base.upper(),
        bt.alt_bases,
        batchinfo.sample_bases,
        batchinfo.strands,
        batchinfo.size
    )

    # base=>[CAF, allele depth], CAF = Allele frequency by read count
    caf = {b: ['%f' % round(bt.depth[b] / float(bt.total_depth), 6),
               bt.depth[b]] for b in bt.alt_bases}

    info = {'CM_DP': str(int(bt.total_depth)),
            'CM_AC': ','.join(map(str, [caf[b][1] for b in bt.alt_bases])),
            'CM_AF': ','.join(map(str, [bt.af_by_lrt[b] for b in bt.alt_bases])),
            'CM_CAF': ','.join(map(str, [caf[b][0] for b in bt.alt_bases])),
            'MQRankSum': str("%.3f" % mq_rank_sum) if mq_rank_sum != -1 else 'nan',
            'ReadPosRankSum': str("%.3f" % read_pos_rank_sum) if read_pos_rank_sum != -1 else 'nan',
            'BaseQRankSum': str("%.3f" % base_q_rank_sum) if base_q_rank_sum != -1 else 'nan',
            'QD': str(qd),
            'SOR': str("%.3f" % sor),
            'FS': str("%.3f" % fs),
            'SB_REF': str(ref_fwd) + ',' + str(ref_rev),
            'SB_ALT': str(alt_fwd) + ',' + str(alt_rev)}

    cdef BaseType g_bt
    cdef bytes group

    if pop_group_bt:
        for group, g_bt in pop_group_bt.items():
            af = ','.join(map(str, [g_bt.af_by_lrt[bb] if bb in g_bt.af_by_lrt else 0
                                    for bb in bt.alt_bases]))
            info[group] = af

    out_file_handle.write('\t'.join([batchinfo.chrid, str(batchinfo.position), '.', batchinfo.ref_base,
                                     ','.join(bt.alt_bases), str(bt.var_qual),
                                     '.' if bt.var_qual > CommonParameter.QUAL_THRESHOLD else 'LowQual',
                                     ';'.join([kk + '=' + vv for kk, vv in sorted(
                                         info.items(), key=lambda x: x[0])]),
                                     'GT:AB:SO:BP'] + samples) + '\n')
    return
