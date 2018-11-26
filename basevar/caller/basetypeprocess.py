"""
This is a Process module for BaseType by mpileup

"""
from .basetype import BaseType
from .algorithm import strand_bias, ref_vs_alt_ranksumtest


def basetypeprocess(chrid, position, ref_base, sample_bases, base_quals, mapqs, strands,
                    indels, read_pos_rank, popgroup, cmm, cvg_file_handle, vcf_file_handle):

    _out_cvg_file(chrid, position, ref_base, sample_bases, strands, indels,
                  popgroup, cvg_file_handle, cmm=cmm)

    # Call variant if ``vcf_file_handle`` is not None
    if vcf_file_handle:

        bt = BaseType(ref_base.upper(), sample_bases, base_quals, cmm=cmm)
        bt.lrt()

        if len(bt.alt_bases()) > 0:

            popgroup_bt = {}
            for group, index in popgroup.items():
                group_sample_bases, group_sample_base_quals = [], []
                for i in index:
                    group_sample_bases.append(sample_bases[i])
                    group_sample_base_quals.append(base_quals[i])

                group_bt = BaseType(ref_base.upper(), group_sample_bases,
                                    group_sample_base_quals, cmm=cmm)

                basecombination = [ref_base.upper()] + bt.alt_bases()
                group_bt.lrt(basecombination)

                popgroup_bt[group] = group_bt

            _out_vcf_line(chrid,
                          position,
                          ref_base,
                          sample_bases,
                          mapqs,
                          read_pos_rank,
                          base_quals,
                          strands,
                          bt,
                          popgroup_bt,
                          vcf_file_handle,
                          cmm=cmm)
    return


def _base_depth_and_indel(bases, indels, cmm=None):
    # coverage info for each position
    base_depth = {b: 0 for b in cmm.BASE}

    for b in bases:

        # ignore all bases('*') which not match ``cmm.BASE``
        if b in base_depth:
            base_depth[b] += 1

    # deal with indels
    indel_dict = {}
    for ind in indels:

        if len(ind) == 0:
            # non indels
            continue

        indel_dict[ind] = indel_dict.get(ind, 0) + 1

    indel_string = ','.join(
        [k + '|' + str(v) for k, v in indel_dict.items()]) if indel_dict else '.'

    return [base_depth, indel_string]


def _out_cvg_file(chrid, position, ref_base, sample_bases, strands,
                 indels, popgroup, out_file_handle, cmm=None):
    """output coverage information into `out_file_handle`"""

    # coverage info for each position
    base_depth, indel_string = _base_depth_and_indel(sample_bases, indels, cmm=cmm)

    # base depth and indels for each subgroup
    group_cvg = {}
    for group, index in popgroup.items():

        group_sample_bases, group_sample_indels = [], []
        for i in index:
            group_sample_bases.append(sample_bases[i])
            group_sample_indels.append(indels[i])

        bd, ind = _base_depth_and_indel(group_sample_bases, group_sample_indels, cmm=cmm)
        group_cvg[group] = [bd, ind]

    fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = 0, -1, 0, 0, 0, 0
    if sample_bases:
        base_sorted = sorted(base_depth.items(),
                             key=lambda x: x[1],
                             reverse=True)

        b1, b2 = base_sorted[0][0], base_sorted[1][0]
        fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
            ref_base.upper(),
            [b1 if b1 != ref_base.upper() else b2],
            sample_bases,
            strands
        )

    if sum(base_depth.values()):

        group_info = []
        if group_cvg:
            for k in popgroup.keys():
                depth, indel_str = group_cvg[k]
                s = ':'.join([str(depth[b]) for b in cmm.BASE] + [indel_str])
                group_info.append(s)

        out_file_handle.write('\t'.join(
            [chrid, str(position), ref_base, str(sum(base_depth.values()))] +
            [str(base_depth[b]) for b in cmm.BASE] + [indel_string] +
            [str(fs), str(sor), ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev]))] +
            group_info) + '\n')

    return


def _out_vcf_line(chrid, position, ref_base, sample_base, mapqs, read_pos_rank, sample_base_qual,
                  strands, bt, pop_group_bt, out_file_handle, cmm=None):
    """output vcf lines into `out_file_handle`"""

    alt_gt = {b: './'+str(k+1) for k, b in enumerate(bt.alt_bases())}
    samples = []

    for k, b in enumerate(sample_base):

        # For sample FORMAT
        if b != 'N':
            # For the base which not in bt.alt_bases()
            if b not in alt_gt:
                alt_gt[b] = './.'

            gt = '0/.' if b == ref_base.upper() else alt_gt[b]

            samples.append(gt + ':' + b + ':' + strands[k] + ':' +
                           str(round(bt.qual_pvalue[k], 6)))
        else:
            samples.append('./.')  # 'N' base

    # Rank Sum Test for mapping qualities of REF versus ALT reads
    mq_rank_sum = ref_vs_alt_ranksumtest(ref_base.upper(), bt.alt_bases(),
                                         zip(sample_base, mapqs))

    # Rank Sum Test for variant appear position among read of REF versus ALT
    read_pos_rank_sum = ref_vs_alt_ranksumtest(ref_base.upper(), bt.alt_bases(),
                                               zip(sample_base, read_pos_rank))

    # Rank Sum Test for base quality of REF versus ALT
    base_q_rank_sum = ref_vs_alt_ranksumtest(ref_base.upper(), bt.alt_bases(),
                                             zip(sample_base, sample_base_qual))

    # Variant call confidence normalized by depth of sample reads
    # supporting a variant.
    ad_sum = sum([bt.depth[b] for b in bt.alt_bases()])
    qd = round(float(bt.var_qual() / ad_sum), 3)

    # Strand bias by fisher exact test and Strand bias estimated by the
    # Symmetric Odds Ratio test
    fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
        ref_base.upper(), bt.alt_bases(), sample_base, strands)

    # base=>[CAF, allele depth], CAF = Allele frequency by read count
    caf = {b: ['%f' % round(bt.depth[b]/float(bt.total_depth), 6),
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
            'SB_REF': str(ref_fwd)+','+str(ref_rev),
            'SB_ALT': str(alt_fwd)+','+str(alt_rev)}

    if pop_group_bt:

        for group, g_bt in pop_group_bt.items():
            af = ','.join(map(str, [g_bt.af_by_lrt[b] if b in g_bt.af_by_lrt else 0
                                    for b in bt.alt_bases()]))
            info[group] = af

    out_file_handle.write('\t'.join([chrid, str(position), '.', ref_base,
                          ','.join(bt.alt_bases()), str(bt.var_qual()),
                          '.' if bt.var_qual() > cmm.QUAL_THRESHOLD else 'LowQual',
                          ';'.join([k+'='+v for k, v in sorted(
                              info.items(), key=lambda x:x[0])]),
                              'GT:AB:SO:BP'] + samples) + '\n')
    return
