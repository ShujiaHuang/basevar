"""
This is a process module of BaseType for fusion file
"""
import sys
import time
import multiprocessing

import pysam

from . import utils
from .basetype import BaseType
from .algorithm import strand_bias, ref_vs_alt_ranksumtest


class BaseTypeFusionSingleProcess(object):
    """
    simple class to repesent a single process.
    """
    def __init__(self, in_ref_file, in_fusion_files, in_popgroup_file,
                 out_vcf_file, out_cvg_file, regions, cmm=None):
        """
        Store input file, options and output file name.

        Parameters:
            ``in_fusion_files``: Array like
                fusion file is 0-base coordinate system
        ===========

            samples: list like
                A list of sample id
        """
        self.ref_file_hd = pysam.FastaFile(in_ref_file)
        self.in_fusion_files = in_fusion_files
        self.out_vcf_file = out_vcf_file
        self.out_cvg_file = out_cvg_file
        self.cmm = cmm
        self.regions = {}

        # store the region into a dict
        for chrid, start, end in regions:

            if chrid not in self.regions:
                self.regions[chrid] = []

            self.regions[chrid].append([start, end])

        # Cache a batch of file handle which index by tabix
        self.tb_files = [pysam.TabixFile(f) for f in self.in_fusion_files]

        # loading sample'id from the header of fusion files
        self.samples = []
        for i, tf in enumerate(self.tb_files):

            try:
                # get sample ID: '##RG\tSM:SAMPLE_ID'
                header = [h for h in tf.header if h.startswith('##RG\tSM:')][0]
                self.samples.append(header.split(':')[-1])
            except IndexError:

                sys.stderr.write('[ERROR] File header has no sample tag mark '
                                 'by "SM:", Please check %s!' % self.in_fusion_files[i])
                self._close_file()
                sys.exit(1)

        # loading population group
        # group_id => [a list samples_index]
        self.popgroup = {}
        if in_popgroup_file and len(in_popgroup_file):

            tmpdict = {}
            with open(in_popgroup_file) as f:
                # Just two columns: sample_id and group_id
                for line in f:
                    sample_id, group_id = line.strip().split()[0:2]
                    tmpdict[sample_id] = group_id + '_AF'

            for i, s in enumerate(self.samples):

                if s in tmpdict:
                    if tmpdict[s] not in self.popgroup:
                        self.popgroup[tmpdict[s]] = []

                    # record different index of different groups
                    self.popgroup[tmpdict[s]].append(i)

    def _close_file(self):

        # close the reference file
        self.ref_file_hd.close()

        # close all of the fusion files
        for tb in self.tb_files:
            tb.close()

        return

    def run(self):
        """
        Run the process of calling variant and output
        """
        vcf_header = utils.vcf_header_define()
        if self.popgroup:
            for g in self.popgroup.keys():
                g_id = g.split('_')[0]  # ignore '_AF'
                vcf_header.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="Allele '
                                  'frequency in the %s populations calculated base on LRT, in '
                                  'the range (0,1)">' % (g_id, g_id))

        with open(self.out_vcf_file, 'w') as VCF, open(self.out_cvg_file, 'w') as CVG:

            CVG.write('\t'.join(['#CHROM', 'POS', 'REF', 'Depth'] + self.cmm.BASE +
                                ['Indel', 'FS', 'SOR', 'Strand_Coverage(REF_FWD,'
                                 'REF_REV,ALT_FWD,ALT_REV)\n']))

            # set header
            VCF.write('\n'.join(vcf_header) + '\n')
            VCF.write('\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
                                 'INFO\tFORMAT'] + self.samples) + '\n')

            for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):
                # ``regions`` is a 2-D array : [[start1,end1], [start2,end2], ...]
                # ``iter_tokes`` is a list of iterator for each sample's input file
                tmp_region = []
                for p in regions:  # covert to 1d-array
                    tmp_region.extend(p)

                tmp_region = sorted(tmp_region)
                start, end = tmp_region[0], tmp_region[-1]
                iter_tokes = []

                # All the fusion file is 0-base system
                for i, bt in enumerate(self.tb_files):
                    try:
                        iter_tokes.append(bt.fetch(chrid, start-1, end))
                    except ValueError:
                        if self.cmm.debug:
                            sys.stderr.write("# [WARMING] Empty region %s:%d-%d in %s" %
                                             (chrid, start-1, end, self.in_fusion_files[i]))
                        iter_tokes.append('')

                # get sequence of chrid from reference fasta
                fa = self.ref_file_hd.fetch(chrid)

                # Init the sample information
                n = 0
                sample_info = [utils.fetch_next(it) for it in iter_tokes]
                for start, end in regions:

                    sys.stderr.write('[INFO] Fetching info from %d samples in region %s'
                                     ' at %s\n' % (len(iter_tokes),
                                                   chrid + ":" + str(start) + "-" + str(end),
                                                   time.asctime())
                                     )

                    for position in xrange(start, end + 1):

                        if n % 100000 == 0:
                            sys.stderr.write("[INFO] loading lines %d at position %s:%d\t%s\n" %
                                             (n+1, chrid, position, time.asctime()))

                        n += 1
                        # sample_base, sample_base_qual, strands, mapqs and
                        # read_pos_rank are listed the same orde with each other.
                        (sample_bases, sample_base_quals, strands, mapqs, read_pos_rank,
                         indels) = self.fetch_base_by_position(
                            fa,  # Fa sequence for indel sequence
                            position - 1,  # postion for fusion file is 0-base
                            sample_info,
                            iter_tokes
                        )

                        ref_base = fa[position-1]

                        # ignore positions if coverage=0 or ref base is 'N' base
                        if (not sample_bases) or (ref_base in ['N', 'n']):
                            continue

                        self._out_cvg_file(chrid, position, ref_base, sample_bases,
                                           strands, indels, CVG)

                        # These two lines just for debug.
                        # VCF.write("\t".join([chrid, str(position), ".", ref_base, ".\t.\t.\t.\t."]+[":".join(map(str, [a,b,c])) for a,b,c in zip(sample_bases,sample_base_quals,strands)]))
                        # sys.exit(1)

                        if ref_base.upper() not in ['A', 'C', 'G', 'T', 'N']:
                            continue

                        bt = BaseType(ref_base.upper(), sample_bases,
                                      sample_base_quals, cmm=self.cmm)
                        bt.lrt()

                        if len(bt.alt_bases()) > 0:

                            popgroup_bt = {}
                            for group, index in self.popgroup.items():
                                group_sample_bases = [sample_bases[i] for i in index]
                                group_sample_base_quals = [sample_base_quals[i] for i in index]

                                group_bt = BaseType(ref_base.upper(), group_sample_bases,
                                                    group_sample_base_quals, cmm=self.cmm)
                                basecombination = [ref_base.upper()] + bt.alt_bases()
                                group_bt.lrt(basecombination)

                                popgroup_bt[group] = group_bt

                            self._out_vcf_line(chrid,
                                               position,
                                               ref_base,
                                               sample_bases,
                                               mapqs,
                                               read_pos_rank,
                                               sample_base_quals,
                                               strands,
                                               bt,
                                               popgroup_bt,
                                               VCF)

        self._close_file()

    def fetch_base_by_position(self, fa, position, sample_info, iter_tokes):
        """
        sample_info example: 0-base
            '#CHROM  START  END  TYPE  MAPQ  SO  Read_POS   BASE_QUAL'
            'chr3	1002728	1002776	<NON_REF>	60	+	16	G'
        """

        base_quals = []
        bases = []
        strands = []
        mapqs = []
        read_pos_rank = []
        indels = []

        for i, sample_pos_line in enumerate(sample_info):
            # 'chr3	1002728	1002776	<NON_REF>	60	+	16	G'

            bs, qs, strand, mapq, rpr, indel, sample_info[i] = (
                self.seek_position(fa, position, sample_pos_line, iter_tokes[i])
            )

            bases.append(bs)
            base_quals.append(qs)
            strands.append(strand)
            mapqs.append(mapq)
            read_pos_rank.append(rpr)

            if indel:
                indels.append(indel.upper())

        return bases, base_quals, strands, mapqs, read_pos_rank, indels

    def seek_position(self, fa, target_pos, sample_line, sample_iter):
        """
        Get mapping info for specific position.
            `fa`: Just for scan indel
        """
        base, strand, indel, qual, rpr, mapq = 'N', '.', '', '!', 0, 0  # Init
        if sample_line:
            # 'chr3	1002728	1002776	<NON_REF>	60	+	16	G'
            tok = sample_line.strip().split()
            start = int(tok[1])
            end = int(tok[2])

            # target_pos must be in the [``start``, ``end``), and not allow == ``end``
            # cause this it's the 0-base coordinate system
            if target_pos >= end:

                while target_pos >= end:

                    sample_line = utils.fetch_next(sample_iter)
                    if sample_line:

                        tok = sample_line.strip().split()
                        start = int(tok[1])
                        end = int(tok[2])
                    else:
                        # The end of file. Break the loop.
                        break

            # In case of sample_line may hit the end of file
            if sample_line and target_pos >= start:

                mapq = int(tok[4])
                strand = tok[5]
                rpr = target_pos - start + int(tok[6])
                qual = tok[7]

                base = tok[3]
                if base == '<NON_REF>':
                    # reference base.
                    # target_pos is 0-base
                    base = fa[target_pos].upper()
                elif base == '.':
                    # deletion
                    indel = '-' + fa[start:end]

                elif len(base) > 1:
                    # insertion
                    indel = '+' + base

        qual = ord(qual) - 33
        return base, qual, strand, mapq, rpr, indel, sample_line

    def _out_cvg_file(self, chrid, position, ref_base, sample_bases,
                      strands, indels, out_file_handle):
        # coverage info for each position

        base_depth = {b: 0 for b in self.cmm.BASE}

        for k, b in enumerate(sample_bases):

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

            out_file_handle.write('\t'.join(
                [chrid, str(position), ref_base, str(sum(base_depth.values()))] +
                [str(base_depth[b]) for b in self.cmm.BASE] + [indel_string]) +
                      '\t' + str(fs) + '\t' + str(sor) + '\t' +
                      ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev])) + '\n')

        return

    def _out_vcf_line(self, chrid, position, ref_base, sample_base, mapqs, read_pos_rank,
                      sample_base_qual, strands, bt, pop_group_bt, out_file_handle):

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
                              '.' if bt.var_qual() > self.cmm.QUAL_THRESHOLD else 'LowQual',
                              ';'.join([k+'='+v for k, v in sorted(
                                  info.items(), key=lambda x:x[0])]),
                                  'GT:AB:SO:BP'] + samples) + '\n')
        return self


###############################################################################
class BaseVarFusionMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.
    """
    def __init__(self, ref_in_file, aligne_files, pop_group_file,
                 out_vcf_file, out_cvg_file, regions, cmm=None):
        """
        Constructor.

        regions: 2d-array like, required
                It's region info , format like: [[chrid, start, end], ...]
        """
        multiprocessing.Process.__init__(self)

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.single_process = BaseTypeFusionSingleProcess(ref_in_file,
                                                          aligne_files,
                                                          pop_group_file,
                                                          out_vcf_file,
                                                          out_cvg_file,
                                                          regions,
                                                          cmm=cmm)

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()

