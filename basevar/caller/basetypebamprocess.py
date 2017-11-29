"""
This is a Process module for BaseType by BAM/CRAM

"""
import sys
import time
import multiprocessing

import pysam

from . import utils
from .basetype import BaseType
from . import bam
from .algorithm import strand_bias, ref_vs_alt_ranksumtest


class BaseVarSingleProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """
    def __init__(self, ref_file, aligne_files, out_vcf_file, out_cvg_file,
                 regions, samples, cmm=None):
        """
        Store input file, options and output file name.

        Parameters:
        ===========

            samples: list like
                    A list of sample id
        """
        self.ref_file_hd = pysam.FastaFile(ref_file)
        self.aligne_files = aligne_files
        self.out_vcf_file = out_vcf_file
        self.out_cvg_file = out_cvg_file
        self.samples = samples
        self.cmm = cmm
        self.regions = {}

        # store the region into a dict
        for chrid, start, end in regions:

            if chrid not in self.regions:
                self.regions[chrid] = []

            self.regions[chrid].append([start, end])

        # Cache a batch of aligne file handle
        self.ali_files_hd = []
        for f in self.aligne_files:

            if f.endswith('.bam'):
                bf = pysam.AlignmentFile(f, 'rb')
            elif f.endswith('.cram'):
                bf = pysam.AlignmentFile(f, 'rc')
            else:
                sys.stderr.write('[ERROR] Input file: %s is not BAM nor CRAM.\n' % f)
                self._close_aligne_file()
                sys.exit(1)

            self.ali_files_hd.append(bf)

    def _close_aligne_file(self):
        self.ref_file_hd.close()
        for bf in self.ali_files_hd:
            bf.close()

    def run(self):
        """
        Run the process of calling variant and output
        """
        vcf_header = utils.vcf_header_define()
        with open(self.out_vcf_file, 'w') as VCF, open(self.out_cvg_file, 'w') as CVG:

            CVG.write('\t'.join(['#CHROM', 'POS', 'REF', 'Depth'] + self.cmm.BASE +
                                ['Indel', 'FS', 'SOR', 'Strand_Coverage(REF_FWD,'
                                                       'REF_REV,ALT_FWD,ALT_REV)\n']))

            # set header
            VCF.write('\n'.join(vcf_header) + '\n')
            VCF.write('\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
                                 'INFO\tFORMAT'] + self.samples) + '\n')

            for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):
                # ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
                # ``iter_tokes`` is a list of iterator for each sample's input file
                tmp_region = []
                for p in regions:  # covert to 1d-array
                    tmp_region.extend(p)

                tmp_region = sorted(tmp_region)

                start, end = tmp_region[0], tmp_region[-1]
                iter_tokes = []
                sample_info = []

                for i, bf in enumerate(self.ali_files_hd):
                    try:
                        # 0-base
                        iter_tokes.append(bf.pileup(chrid, start-1, end))
                    except ValueError:
                        if self.cmm.debug:
                            print >> sys.stderr, ("# [WARMING] Empty region",
                                                  chrid, start-1, end,
                                                  self.aligne_files[i])
                        iter_tokes.append('')

                # get all the reference sequence of chrid
                fa = self.ref_file_hd.fetch(chrid)

                # Set iteration marker: 1->iterate; 0->do not
                # iterate or hit the end
                go_iter = [1] * len(iter_tokes)
                n = 0
                for start, end in regions:
                    for position in xrange(start, end + 1):

                        if n % 100000 == 0:
                            sys.stderr.write("[INFO] loading lines %d at position %s:%d\t%s\n" %
                                             (n, chrid, position, time.asctime()))
                        n += 1

                        sample_info = [utils.fetch_next(iter_tokes[i]) if g else sample_info[i]
                                       for i, g in enumerate(go_iter)]

                        # sample_base, sample_base_qual, strands, mapqs and
                        # read_pos_rank are listed the same orde with each other.
                        (sample_base, sample_base_qual, strands, mapqs, read_pos_rank,
                         indels) = bam.fetch_base_by_position(
                            position - 1,  # postion in pysam is 0-base
                            sample_info,
                            go_iter,
                            iter_tokes,
                            fa,  # Fa sequence for indel sequence
                            is_scan_indel=True
                        )

                        ref_base = fa[position-1]
                        # ignore positions if coverage=0 or ref base is 'N' base
                        if not sample_base or ref_base in ['N', 'n']:
                            continue

                        self._out_cvg_file(chrid, position, ref_base, sample_base,
                                           strands, indels, CVG)

                        bt = BaseType(ref_base.upper(), sample_base,
                                      sample_base_qual, cmm=self.cmm)
                        bt.lrt()
                        if len(bt.alt_bases()) > 0:
                            self._out_vcf_line(chrid,
                                               position,
                                               ref_base,
                                               sample_base,
                                               mapqs,
                                               read_pos_rank,
                                               sample_base_qual,
                                               strands,
                                               bt,
                                               VCF)

        self._close_aligne_file()

    def _out_cvg_file(self, chrid, position, ref_base, sample_base,
                      strands, indels, out_file_handle):
        # coverage info for each position

        base_depth = {b: 0 for b in self.cmm.BASE}
        for k, b in enumerate(sample_base):

            # ignore all bases('*') which not match ``cmm.BASE``
            if b in base_depth:
                base_depth[b] += 1

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
                ref_base.upper(),
                [b1 if b1 != ref_base.upper() else b2],
                sample_base,
                strands
            )

        if sum(base_depth.values()):
            out_file_handle.write('\t'.join(
                [chrid, str(position), ref_base, str(sum(base_depth.values()))] +
                [str(base_depth[b]) for b in self.cmm.BASE] + [indel_string]) +
                      '\t' + str(fs) + '\t' + str(sor) + '\t' +
                      ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev])) + '\n')

        return

    def _out_vcf_line(self, chrid, position, ref_base, sample_base, mapqs,
                      read_pos_rank, sample_base_qual, strands, bt, out_file_handle):

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

        # base=>[AF, allele depth]
        af = {b: ['%f' % round(bt.depth[b] / float(bt.total_depth), 6),
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

        out_file_handle.write('\t'.join([chrid, str(position), '.', ref_base,
                              ','.join(bt.alt_bases()), str(bt.var_qual()),
                              '.' if bt.var_qual() > self.cmm.QUAL_THRESHOLD else 'LowQual',
                              ';'.join([k+'='+v for k, v in sorted(
                                  info.items(), key=lambda x:x[0])]),
                                  'GT:AB:SO:BP'] + samples) + '\n')
        return self


###############################################################################
class BaseVarMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.
    """
    def __init__(self, ref_in_file, aligne_files, out_vcf_file, out_cvg_file,
                 regions, cmm=None):
        """
        Constructor.

        regions: 2d-array like, required
                It's region info , format like: [[chrid, start, end], ...]
        """
        multiprocessing.Process.__init__(self)

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        samples_id = self._load_sample_id(aligne_files)

        self.single_process = BaseVarSingleProcess(ref_in_file,
                                                   aligne_files,
                                                   out_vcf_file,
                                                   out_cvg_file,
                                                   regions,
                                                   samples_id,
                                                   cmm=cmm)

    def _load_sample_id(self, aligne_files):
        """loading sample id in bam/cram files from RG tag"""

        sample_id = []
        for al in aligne_files:
            bf = pysam.AlignmentFile(al)

            if 'RG' not in bf.header:
                sys.stderr.write('[ERROR] Bam file format error: missiong '
                                 '@RG in the header.\n')
                bf.close()
                sys.exit(1)

            sample_id.append(bf.header['RG'][0]['SM'])
            bf.close()

        return sample_id

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()
