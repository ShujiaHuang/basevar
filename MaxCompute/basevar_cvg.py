#coding:utf-8
from odps.udf import annotate
from odps.distcache import get_cache_archive

from mpileup import first_base
from algorithm import strand_bias

def include_package_path(res_name):
    import os, sys
    archive_files = get_cache_archive(res_name)
    dir_names = sorted([os.path.dirname(os.path.normpath(f.name)) for f in archive_files
                       if '.dist_info' not in f.name], key=lambda v: len(v))
    sys.path.append(os.path.dirname(dir_names[0]))

@annotate("string,string,string,string->string")
class BaseVarCoverage(object):

    def __init__(self):
        include_package_path('scipy.zip')

    def evaluate(self, chrid, pos, base_ref, one):
        print 'processing %s %s %s' % (chrid, pos, base_ref)
        bases = []
        quals = []
        strands = []
        indels = []

        for b, q, s in zip(*[iter(one.split('\t'))]*3):
            if b != '0' and q != '*':
            # TODO hardcode is_scan_indel=False
               strand, base, qual, indel = first_base(
                   base_ref, q, s, is_scan_indel=False)
               bases.append(base)
               quals.append(ord(qual) - 33)
               strands.append(strand)
               indels.extend(indel)

        return self._out_cvg_line(chrid, pos, base_ref, bases, strands, indels)


    def _out_cvg_line(self, chrid, position, ref_base, sample_base,
                      strands, indels):

        # TODO
        self.BASE = ['A', 'C', 'G', 'T']
        self.total_subsamcol = None

        # coverage info for each position
        base_depth = {b: 0 for b in self.BASE}
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
