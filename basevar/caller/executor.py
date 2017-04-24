"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``BaseVar.py`` is "Sauron", and this module could just be called by it.
"""

import sys
import os
import argparse
import heapq
import pysam

from . import utils
from . import mpileup
from .variantcaller import BaseType
from .algorithm import strand_bias


class FileForQueueing(object):
    """
    """
    def __init__(self, the_file, line):
        """
        Store the file, and initialise the current value
        """
        self.the_file = the_file
        self.finishedReadingFile = False
        self.heap = []

        line = line
        cols = line.strip().split("\t")
        chrom = cols[0]

        # Where possible, convert chromosome names into
        # integers for sorting. If not possible, use
        # original names.
        try:
            chrom = int(chrom.upper().strip("CHR"))
        except Exception:
            pass

        pos = int(cols[1])
        heapq.heappush(self.heap, (chrom, pos, line))

        while not self.finishedReadingFile and len(self.heap) < 100:

            try:
                line = self.the_file.next()
                cols = line.strip().split("\t")
                chrom = cols[0]

                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
            except StopIteration:
                self.finishedReadingFile = True
                break

            heapq.heappush(self.heap, (chrom, pos, line))

        # Now take the top line
        self.chrom, self.pos, self.line = heapq.heappop(self.heap)

    def __cmp__(self, other):
        """
        Comparison function. Utilises the comparison function defined in
        the AlignedRead class.
        """
        return cmp(self.chrom, other.chrom) or cmp(self.pos, other.pos)

    def __del__(self):
        """
        Destructor
        """
        self.the_file.close()
        os.remove(self.the_file.name)

    def next(self):
        """
        Increment the iterator and yield the new value. Also, store the
        current value for use in the comparison function.
        """
        if not self.finishedReadingFile:

            try:
                line = self.the_file.next()
                #cols = line.strip().split('\t')
                cols = line.strip().split()
                chrom = cols[0]

                # Where possible, convert chromosome names into
                # integers for sorting. If not possible, use
                # original names.
                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
                heapq.heappush(self.heap, (chrom, pos, line))

            except StopIteration:
                self.finishedReadingFile = True

        if len(self.heap) != 0:
            # Now take the top line
            self.chrom, self.pos, self.line = heapq.heappop(self.heap)
        else:
            raise StopIteration


def mergeVCFFiles(temp_file_names, final_file_name, log):
    """
    """
    log.info("Merging output VCF file(s) into final file %s" %(final_file_name))

    # Final output file
    if final_file_name == "-":
        outputVCF = sys.stdout
    else:
        outputVCF = utils.Open(final_file_name, 'wb')
    the_heap = []

    # Initialise queue
    for index, file_name in enumerate(temp_file_names):
        the_file = utils.Open(file_name, 'rb')

        for line in the_file:

            # End of this file
            if line[0] == "#":
                if index == 0:
                    outputVCF.write(line)
            else:
                the_file_for_queueing = FileForQueueing(the_file, line)
                heapq.heappush(the_heap, the_file_for_queueing)
                break

        # If there are no calls in the temp file, we still want to
        # remove it.
        else:
            the_file.close()
            os.remove(file_name)

    # Merge-sort the output using a priority queue
    while len(the_heap) != 0:

        # Get file from heap in right order
        next_file = heapq.heappop(the_heap)
        outputVCF.write(next_file.line)

        # Put file back on heap
        try:
            next_file.next()
            heapq.heappush(the_heap, next_file)
        except StopIteration:
            continue

    # Close final output file
    if final_file_name != "-":
        outputVCF.close()

    log.info("Finished merging VCF file(s)")


class Runner(object):

    def __init__(self, cmm=utils.CommonParameter()):
        """init function
        """
        self.cmm = cmm

    def _common_init(self, optp):
        """Common init function for getting positions and region.

        :param optp:
        :return:
        """
        # common parameters
        optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                          help='skip unlisted positions (chr pos)', default='')
        optp.add_argument('-R', '--regions',
                          metavar='chr:start-end', dest='regions',
                          help='skip positions not in (chr:start-end)', default='')
        optp.add_argument('-l', '--mpileup-list', dest='infilelist', metavar='FILE',
                          help='The input mpileup file list.', default='')
        optp.add_argument('-s', '--sample-list', dest='samplelistfile',
                          metavar='FILE', help='The sample list.')
        optp.add_argument('-S', '--subsample-list', dest='subsample', metavar='FILE',
                          help='Skip samples not in subsample-list, one sample per row.')

        opt = optp.parse_args()

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] At least one mpileup to input\n')

        if len(opt.samplelistfile) == 0:
            optp.error('[ERROR] Must input the sample\'s ID list file by (-s)')

        if len(opt.positions) == 0 and len(opt.regions) == 0:
            optp.error('[ERROR] The list of position (-L or -R) is required.\n')

        # Loading positions
        _sites = utils.get_list_position(opt.positions) if opt.positions else {}
        if len(opt.regions):
            chrid, reg = opt.regions.strip().split(':')
            reg = map(int, reg.split('-'))
            if chrid not in _sites:
                _sites[chrid] = []

            _sites[chrid].append([reg[0], reg[1]])

        # merge and sorted the regions
        self.sites = {k:utils.merge_region(v) for k, v in _sites.items()}

        # Load all the mpileup files
        self.files = [f for f in sys.argv if '.mpileup.gz' in f]
        if opt.infilelist:
            self.files.extend(utils.load_file_list(opt.infilelist))

        # open mpileup files which store a batch of sample info by tabix
        self.tb_files = [pysam.TabixFile(f) for f in self.files]

        # load all the samples, 2D array
        self.sample_id = []
        with open(opt.samplelistfile) as I:

            samplefiles = []
            for r in I:
                if r.startswith('#'): continue
                samplefiles.append(r.strip().split()[0])

            for f in samplefiles:
                with open(f) as I:
                    self.sample_id.append([s.strip().split()[0] for s in I])

        self.total_sample = []
        for s in self.sample_id:
            self.total_sample.extend(s)

        # loading subsample if provide
        self.total_subsamcol = []
        if opt.subsample:
            subsample = []
            with open(opt.subsample) as I:
                for r in I:
                    if r.startswith('#'): continue
                    subsample.append(r.strip().split()[0])

            subsample = set(subsample)
            for i, s in enumerate(self.total_sample):
                if s in subsample:
                    self.total_subsamcol.append(i) # get index in total_sample

            self.total_subsamcol = set(self.total_subsamcol)

        return opt

    # def _fetch_base_by_position(self, position, sample_info, go_iter, iter_tokes,
    #                             is_scan_indel=False):
    #
    #     sample_base_qual = []
    #     sample_base = []
    #     strands = []
    #     indels = []
    #     ref_base = ''
    #     for i, tb_sample_line in enumerate(sample_info):
    #
    #         sample_info[i], ref_base_t, bs, qs, strand, go_iter[i], indel = (
    #             mpileup.seek_position(position, tb_sample_line,
    #                                   len(self.sample_id[i]),
    #                                   iter_tokes[i],
    #                                   is_scan_indel=is_scan_indel)
    #         )
    #
    #         sample_base.extend(bs)
    #         strands.extend(strand)
    #         sample_base_qual.extend([ord(q) - 33 for q in qs])
    #         indels.extend(indel)
    #
    #         if not ref_base:
    #             ref_base = ref_base_t
    #
    #     return ref_base, sample_base, sample_base_qual, strands, indels

    def basetype(self):

        optp = argparse.ArgumentParser()
        optp.add_argument('basetype')
        optp.add_argument('-m', '--min_af', dest='min_af', type='float',
                          metavar='MINAF', default=0.001,
                          help='The effective base frequence threshold. [0.001]')
        optp.add_argument('-o', '--outprefix', dest='outprefix',
                          metavar='FILE', default='out',
                          help='The prefix of output files. [out]')

        self.opt = self._common_init(optp)
        self.cmm.MINAF = self.opt.min_af

        # output file name
        self.out_vcf_file = self.opt.outprefix + '.vcf'
        self.out_cvg_file = self.opt.outprefix + '.cvg.tsv' # position coverage

        vcf_header = utils.vcf_header_define()
        with open(self.out_vcf_file, 'w') as VCF, open(self.out_cvg_file, 'w') as CVG:

            CVG.write('\t'.join(['#CHROM', 'POS', 'REF', 'Depth'] + self.cmm.BASE +
                                ['Indel', 'FS', 'Strand_cvg']) + '\n')

            VCF.write('\n'.join(vcf_header) + '\n')
            VCF.write('\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
                                 'INFO\tFORMAT'] + self.total_sample) + '\n')

            for chrid, regions in sorted(self.sites.items(), key=lambda x:x[0]):
                # ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
                # fetch the position data from each mpileup files
                # `iter_tokes` is a list of iterator for each sample's mpileup
                tmp_region = []
                for p in regions: tmp_region.extend(p) # 1d-array
                tmp_region = sorted(tmp_region)

                start, end = tmp_region[0], tmp_region[-1]
                iter_tokes = []
                sample_info = []
                for i, tb in enumerate(self.tb_files):
                    try:
                        iter_tokes.append(tb.fetch(chrid, start-1, end))
                    except ValueError:
                        if self.cmm.debug:
                            print >> sys.stderr, ("# [WARMING] Empty region",
                                                  chrid, start-1, end,
                                                  self.files[i])
                        iter_tokes.append('')

                # Set iteration marker: 1->iterate; 0->donot
                # iterate or hit the end
                go_iter = [1] * len(iter_tokes)
                for start, end in regions:
                    for position in xrange(start, end+1):

                        sample_info = [mpileup.fetch_next(iter_tokes[i])
                                       if g else sample_info[i]
                                       for i, g in enumerate(go_iter)]

                        ref_base, sample_base, sample_base_qual, strands, indels = (
                            mpileup.fetch_base_by_position(position, self.sample_id,
                                                           sample_info, go_iter,
                                                           iter_tokes,
                                                           is_scan_indel=True)
                        )

                        if self.total_subsamcol:
                            for k, b in enumerate(sample_base):
                                if k not in self.total_subsamcol:
                                    # set un-selected bases to be 'N' which
                                    # will be filted later
                                    sample_base[k] = 'N'

                        bt = BaseType(ref_base.upper(), sample_base, sample_base_qual,
                                      cmm=self.cmm)
                        bt.lrt()

                        # ACGT count and mark the refbase
                        if not ref_base:
                            # mark '*' if coverage is 0
                            ref_base = '*'

                        self._out_cvg_file(chrid, position, ref_base,
                                           sample_base, strands, indels, CVG)

                        if len(bt.alt_bases()) > 0:
                            self._out_vcf_line(chrid, position, ref_base,
                                               sample_base, strands, bt, VCF)

        self._close_tabix()

        return

    def coverage(self):

        optp = argparse.ArgumentParser()
        optp.add_argument('converage')
        optp.add_argument('-o', '--outfile', dest='outfile',
                          metavar='FILE', default='out.tsv',
                          help='The name of output file. [out.tsv]')
        self.opt = self._common_init(optp)

        # output file name, position coverage
        self.out_cvg_file = self.opt.outfile

        with open(self.out_cvg_file, 'w') as CVG:

            CVG.write('\t'.join(['#CHROM','POS', 'REF', 'Depth'] +
                                self.cmm.BASE + ['Indel', 'FS', 'Strand_cvg'])+ '\n')

            for chrid, regions in sorted(self.sites.items(), key=lambda x: x[0]):
                # ``regions`` is a 2-D array: [[start1,end1], [start2,end2], ...]
                # fetch the position data from each mpileup files
                # `iter_tokes` is a list of iterator for each sample's mpileup
                tmp_region = []
                for p in regions: tmp_region.extend(p)
                tmp_region = sorted(tmp_region)
                start, end = tmp_region[0], tmp_region[-1]

                iter_tokes = []
                sample_info = []

                for i, tb in enumerate(self.tb_files):
                    try:
                        iter_tokes.append(tb.fetch(chrid, start-1, end))
                    except ValueError:
                        if self.cmm.debug:
                            print >> sys.stderr, ("# [WARMING] Empty region",
                                                  chrid, start-1, end,
                                                  self.files[i])
                        iter_tokes.append('')

                # Set iteration marker: 1->iterate; 0->hit the end
                go_iter = [1] * len(iter_tokes)
                for start, end in regions:
                    for position in xrange(start, end+1):

                        sample_info = [mpileup.fetch_next(iter_tokes[i])
                                       if g else sample_info[i]
                                       for i, g in enumerate(go_iter)]

                        ref_base, sample_base, sample_base_qual, strands, indels = (
                            mpileup.fetch_base_by_position(position, self.sample_id,
                                                           sample_info, go_iter,
                                                           iter_tokes,
                                                           is_scan_indel=True)
                        )

                        # ACGT count and mark refbase
                        if not ref_base:
                            # mark '*' if the coverage is 0
                            ref_base = '*'

                        self._out_cvg_file(chrid, position, ref_base,
                                           sample_base, strands, indels, CVG)

        self._close_tabix()

        return

    def _out_cvg_file(self, chrid, position, ref_base, sample_base,
                      strands, indels, out_file_handle):
        # coverage info for each position

        base_depth = {b: 0 for b in self.cmm.BASE}
        for k, b in enumerate(sample_base):

            if self.total_subsamcol and k not in self.total_subsamcol:
                # set un-selected bases to be 'N' which will be filted later
                sample_base[k] = 'N'
                continue

            # ignore all bases which not match ``cmm.BASE``
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

        out_file_handle.write('\t'.join(
            [chrid, str(position), ref_base, str(sum(base_depth.values()))] +
            [str(base_depth[b]) for b in self.cmm.BASE] + [indel_string]) +
                  '\t' + str(fs) + '\t' +
                  ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev])) + '\n')

        return

    def _out_vcf_line(self, chrid, position, ref_base, sample_base,
                      strands, bt, out_file_handle):
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

        out_file_handle.write('\t'.join([chrid, str(position), '.', ref_base,
                         ','.join(bt.alt_bases()), str(bt.var_qual()),
                         '.' if bt.var_qual() > self.cmm.QUAL_THRESHOLD else 'LowQual',
                         ';'.join([k+'='+v for k, v in sorted(
                            info.items(), key=lambda x:x[0])]),
                            'GT:AB:SO:BP'] + samples) + '\n')
        return

    def _close_tabix(self):
        for tb in self.tb_files:
            tb.close()



