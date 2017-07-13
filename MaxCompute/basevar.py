#coding:utf-8
from odps.udf import annotate
from odps.udf import BaseUDTF
from odps.distcache import get_cache_archive

from mpileup import first_base
from algorithm import strand_bias
from basetype import CommonParameter
from basetype import BaseType

# for import scipy
def include_package_path(res_name):
    import os, sys
    archive_files = get_cache_archive(res_name)
    dir_names = sorted([os.path.dirname(os.path.normpath(f.name)) for f in archive_files
                       if '.dist_info' not in f.name], key=lambda v: len(v))
    sys.path.append(os.path.dirname(dir_names[0]))

@annotate("string,string,string,string,string->string")
class BaseVar(BaseUDTF):

    def __init__(self):
        self.cmm = CommonParameter()
        include_package_path('scipy.zip')

    def process(self, mode, chrid, pos, base_ref, one):
        print '%s processing %s %s %s' % (mode, chrid, pos, base_ref)
        bases = []
        quals = []
        strands = []
        indels = []

        for b, q, s in zip(*[iter(one.split('\t'))]*3):
            if b != '0' and q != '*':
            # TODO hardcode is_scan_indel=False
               strand, base, qual, indel = first_base(
                   base_ref, q, s, is_scan_indel=self.cmm.scan_indel)
               bases.append(base)
               quals.append(ord(qual) - 33)
               strands.append(strand)
               indels.extend(indel)

        if mode == 'coverage':
            self.forward(self._out_cvg_line(chrid, pos, base_ref, bases, strands, indels))
        elif mode == 'vcf':
            bt = BaseType(base_ref.upper(), bases, quals, cmm=self.cmm)
            bt.lrt()
            if len(bt.alt_bases()) > 0:
                self.forward(self._out_vcf_line(chrid, pos, base_ref,
                                                bases, strands, bt))
        else:
            raise Exception('unknown mode %s' % mode)

    def _out_cvg_line(self, chrid, position, ref_base, sample_base,
                      strands, indels):

        # TODO
        self.total_subsamcol = None

        # coverage info for each position
        base_depth = {b: 0 for b in self.cmm.BASE}
        for k, b in enumerate(sample_base):

            if self.total_subsamcol and k not in self.total_subsamcol:
                # set un-selected bases to be 'N' which will be filted later
                sample_base[k] = 'N'
                continue

            # ignore all bases('*') which not match ``cmm.BASE``
            if b in base_depth:
                base_depth[b] += 1

        # deal with indels
        indel_dict = {}
        for ind in indels:
            indel_dict[ind] = indel_dict.get(ind, 0) + 1

        indel_string = ','.join([k + ':' + str(v)
                                 for k, v in indel_dict.items()]) if indel_dict else '.'

        fs, ref_fwd, ref_rev, alt_fwd, alt_rev = 0, 0, 0, 0, 0
        if sample_base:
            base_sorted = sorted(base_depth.items(),
                                 lambda x, y: cmp(x[1], y[1]),
                                 reverse=True)

            b1, b2 = base_sorted[0][0], base_sorted[1][0]
            fs, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
                ref_base, [b1 if b1 != ref_base.upper() else b2],
                sample_base, strands)

        return '\t'.join([chrid, str(position), ref_base, str(sum(base_depth.values()))] + [str(base_depth[b]) for b in self.BASE] + [indel_string]) + '\t' + str(fs) + '\t' + ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev]))

    def _out_vcf_line(self, chrid, position, ref_base, sample_base,
                      strands, bt):
        #
        alt_gt = {b:'./'+str(k+1) for k,b in enumerate(bt.alt_bases())}
        samples = []

        for k, b in enumerate(sample_base):

            # For sample FORMAT
            if b != 'N':
                # For the base which not in bt.alt_bases()
                if b not in alt_gt: alt_gt[b] = './.'
                gt = '0/.' if b==ref_base.upper() else alt_gt[b]

                samples.append(gt+':'+b+':'+strands[k]+':'+
                               str(round(bt.qual_pvalue[k], 6)))
            else:
                samples.append('./.') ## 'N' base

        # Strand bias by fisher exact test
        # Normally you remove any SNP with FS > 60.0 and an indel with FS > 200.0
        fs, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
            ref_base, bt.alt_bases(), sample_base, strands)

        # base=>[AF, allele depth]
        af = {b:['%f' % round(bt.depth[b]/float(bt.total_depth), 6),
                 bt.depth[b]] for b in bt.alt_bases()}

        info = {'CM_DP': str(int(bt.total_depth)),
                'CM_AC': ','.join(map(str, [af[b][1] for b in bt.alt_bases()])),
                'CM_AF': ','.join(map(str, [af[b][0] for b in bt.alt_bases()])),
                'CM_EAF': ','.join(map(str, [bt.eaf[b] for b in bt.alt_bases()])),
                'FS': str(fs),
                'SB_REF': str(ref_fwd)+','+str(ref_rev),
                'SB_ALT': str(alt_fwd)+','+str(alt_rev)}

        return '\t'.join([chrid, str(position), '.', ref_base,
                         ','.join(bt.alt_bases()), str(bt.var_qual()),
                         '.' if bt.var_qual() > self.cmm.QUAL_THRESHOLD else 'LowQual',
                         ';'.join([k+'='+v for k, v in sorted(
                            info.items(), key=lambda x:x[0])]),
                            'GT:AB:SO:BP'] + samples)
