#!/usr/bin/env python

#coding:utf-8
from odps.udf import annotate
from odps.udf import BaseUDTF
from odps.distcache import get_cache_archive

from algorithm import strand_bias, ref_vs_alt_ranksumtest
from basetype import CommonParameter
from basetype import BaseType

# for import scipy
def include_package_path(res_name):
    import os, sys
    sys.path.append('work/scipy.zip')

@annotate("string,string,string,string,string,string,string->string")
class BaseVar(BaseUDTF):

    def __init__(self, cmm=CommonParameter()):
        self.cmm = cmm
        if not cmm.debug:
            include_package_path('scipy.zip')
        res_file = get_cache_file('target_sample.txt')
        name_idx = dict()
        i = 0
        for line in res_file:
            name_idx[line.strip()] = i
            i += 1
        res_file.close()
        popgroup_file = get_cache_file('popgroup.txt')
        self.popgroup = {}
        for line in popgroup_file:
            sample_id, group_id = line.strip().split()[0:2]
            if group_id not in self.popgroup:
                self.popgroup[group_id] = []
            self.popgroup[group_id].append(name_idx[sample_id])

    def process(self, mode, chrid, pos, base_ref, part0, part1, part2):
        #tokens = '\t'.join([part0, part1, part2]).strip().split('\t')
        tokens = part0.strip().split('\t')
        tokens.extend(part1.strip().split('\t'))
        tokens.extend(part2.strip().split('\t'))

        sample_count = len(tokens) / 6

        bases = ['N'] * sample_count
        quals = [0] * sample_count
        strands = ['.'] * sample_count
        mapqs = [0] * sample_count
        read_pos_ranks = [0] * sample_count
        indels = []

        for i in xrange(0, sample_count):
            # TODO b != '0' q != '*'
            read_base       = tokens[i * 6]
            read_quality    = int(tokens[i * 6 + 1])
            mapping_quality = int(tokens[i * 6 + 2])
            read_pos_rank   = int(tokens[i * 6 + 3])
            indel           = tokens[i * 6 + 4]
            strand          = tokens[i * 6 + 5]

            bases[i] = read_base
            quals[i] = read_quality
            mapqs[i] = mapping_quality
            read_pos_ranks[i] = read_pos_rank
            if indel:
                indels.append(indel)
            strands[i] = strand

        if mode == 'coverage':
            cvg_line = self._out_cvg_line(chrid, pos, base_ref, bases, strands, indels)
            if cvg_line:
                self.forward(cvg_line)
        elif mode == 'vcf':
            bt = BaseType(base_ref, bases, quals, cmm=self.cmm)
            bt.lrt()
            if len(bt.alt_bases()) > 0:
                popgroup_bt = {}
                for group, index in self.popgroup.items():
                    group_sample_bases = [bases[i] for i in index]
                    group_sample_base_quals = [quals[i] for i in index]
                    group_bt = BaseType(base_ref.upper(), group_sample_bases, group_sample_base_quals, cmm=self.cmm)
                    basecombination = [base_ref.upper()] + bt.alt_bases()
                    group_bt.lrt(basecombination)
                    popgroup_bt[group] = group_bt
                self.forward(self._out_vcf_line(chrid,
                                                pos,
                                                base_ref,
                                                bases,
                                                mapqs,
                                                read_pos_ranks,
                                                quals,
                                                strands,
                                                bt,
                                                popgroup_bt)
                )
        else:
            raise Exception('unknown mode %s' % mode)

    def _out_cvg_line(self, chrid, position, ref_base, sample_base,
                      strands, indels):
        # coverage info for each position

        has_valid_sample = False
        base_depth = {b: 0 for b in self.cmm.BASE}
        for k, b in enumerate(sample_base):
            if b in base_depth:
                has_valid_sample = True
                base_depth[b] += 1
        if not has_valid_sample:
            return ''
        # deal with indels
        indel_dict = {}
        for ind in indels:
            indel_dict[ind] = indel_dict.get(ind, 0) + 1

        indel_string = ','.join(
            [k + ':' + str(v) for k, v in indel_dict.items()]) if indel_dict else '.'

        fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = 0, -1, 0, 0, 0, 0
        if sample_base:
            base_sorted = sorted(base_depth.items(),
                                 key=lambda x: x[1],
                                 reverse=True)

            b1, b2 = base_sorted[0][0], base_sorted[1][0]
            fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
                ref_base,
                [b1 if b1 != ref_base else b2],
                sample_base,
                strands
            )

        return '\t'.join(
            [chrid, str(position), ref_base, str(sum(base_depth.values()))] +
            [str(base_depth[b]) for b in self.cmm.BASE] + [indel_string]) + '\t' + str(fs) + '\t' + str(sor) + '\t' + ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev]))

    def _out_vcf_line(self, chrid, position, ref_base, sample_base, mapqs,
                      read_pos_ranks, sample_base_qual, strands, bt):

        alt_gt = {b:'./'+str(k+1) for k,b in enumerate(bt.alt_bases())}
        samples = []

        for k, b in enumerate(sample_base):
            # For sample FORMAT
            if b != 'N':
                # For the base which not in bt.alt_bases()
                if b not in alt_gt:
                    alt_gt[b] = './.'

                gt = '0/.' if b == ref_base else alt_gt[b]

                samples.append(gt+':'+b+':'+strands[k]+':'+
                               str(round(bt.qual_pvalue[k], 6)))
            else:
                samples.append('./.') ## 'N' base

        # Rank Sum Test for mapping qualities of REF versus ALT reads
        mq_rank_sum = ref_vs_alt_ranksumtest(ref_base, bt.alt_bases(),
                                             zip(sample_base, mapqs))

        # Rank Sum Test for variant appear position among read of REF versus ALT
        read_pos_rank_sum = ref_vs_alt_ranksumtest(ref_base, bt.alt_bases(),
                                                   zip(sample_base, read_pos_ranks))

        # Rank Sum Test for base quality of REF versus ALT
        base_q_rank_sum = ref_vs_alt_ranksumtest(ref_base, bt.alt_bases(),
                                                 zip(sample_base, sample_base_qual))

        # Variant call confidence normalized by depth of sample reads
        # supporting a variant.
        ad_sum = sum([bt.depth[b] for b in bt.alt_bases()])
        qd = round(float(bt.var_qual() / ad_sum), 3)

        # Strand bias by fisher exact test
        # Normally you remove any SNP with FS > 60.0 and an indel with FS > 200.0
        fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
            ref_base, bt.alt_bases(), sample_base, strands)

        # base=>[AF, allele depth]
        af = {b:['%f' % round(bt.depth[b]/float(bt.total_depth), 6),
                 bt.depth[b]] for b in bt.alt_bases()}

        info = {'CM_DP': str(int(bt.total_depth)),
                'CM_AC': ','.join(map(str, [af[b][1] for b in bt.alt_bases()])),
                'CM_AF': ','.join(map(str, [af[b][0] for b in bt.alt_bases()])),
                'CM_EAF': ','.join(map(str, [bt.eaf[b] for b in bt.alt_bases()])),
                'MQRankSum': str(mq_rank_sum),
                'ReadPosRankSum': str(read_pos_rank_sum),
                'BaseQRankSum': str(base_q_rank_sum),
                'QD': str(qd),
                'SOR': str(sor),
                'FS': str(fs),
                'SB_REF': str(ref_fwd)+','+str(ref_rev),
                'SB_ALT': str(alt_fwd)+','+str(alt_rev)}

        return '\t'.join([chrid, str(position), '.', ref_base,
                         ','.join(bt.alt_bases()), str(bt.var_qual()),
                         '.' if bt.var_qual() > self.cmm.QUAL_THRESHOLD else 'LowQual',
                         ';'.join([k+'='+v for k, v in sorted(
                            info.items(), key=lambda x:x[0])]),
                            'GT:AB:SO:BP'] + samples)

# for local test
if __name__ == '__main__':
    import sys
    sys.path.append('.')
    if len(sys.argv) < 3:
        print 'usage: %s vcf|coverage input_wide_table_file' % sys.argv[0]
        sys.exit(1)

    mode = sys.argv[1]
    cmm = CommonParameter()
    cmm.debug = True
    basevar = BaseVar(cmm)
    with open(sys.argv[2]) as f:
        for l in f:
            token = l.split(',', 3)
            token[-1] = token[-1].rstrip('\n')
            basevar.process(mode, token[0], token[1], token[2], token[3])
