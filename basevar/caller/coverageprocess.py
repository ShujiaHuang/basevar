"""
This is a Process module for calculatinig coverage

"""
import sys
import time
import multiprocessing

import pysam

from . import utils
from . import bam


def fetch_base_by_position(position, sample_info, go_iter, iter_tokes, fa,
                           is_scan_indel=True):
    bases = []
    strands = []
    indels = []

    for i, sample_pos_line in enumerate(sample_info):
        bs, strand, indel, sample_info[i], go_iter[i] = seek_position(
            position, sample_pos_line, iter_tokes[i], fa, is_scan_indel=is_scan_indel)

        bases.append(bs)
        strands.append(strand)

        if indel:
            indels.append(indel)

    return bases, strands, indels


def seek_position(target_pos, sample_pos_line, sample_iter, fa,
                  is_scan_indel=False):
    """Get mapping info for specific position.

    `fa`: Just for scan indel
    """
    base, strand, indel = [], [], []  # Init
    go_iter_mark = 0  # 1->iterate; 0->donot iterate or hit the end
    if sample_pos_line:

        if sample_pos_line.pos < target_pos:

            pos = sample_pos_line.pos
            while pos < target_pos:

                sample_pos_line = utils.fetch_next(sample_iter)
                if sample_pos_line:
                    pos = sample_pos_line.pos
                else:
                    # The end of file. Break the loop.
                    break

        # In case sample_pos_line may hit the end of file
        if sample_pos_line and sample_pos_line.pos == target_pos:

            go_iter_mark = 1  # keep iterate
            base, strand, indel = all_base(
                sample_pos_line,
                is_scan_indel={'yes': is_scan_indel,
                               'pos': sample_pos_line.pos,
                               'fa': fa}
            )

        else:
            # sample_pos_line.pos > target_pos
            go_iter_mark = 0

    return base, strand, indel, sample_pos_line, go_iter_mark


def all_base(sample_pos_line, is_scan_indel=None):
    """Just get the all base for each sample.
    """
    # set the default value for is_scan_indel
    if not is_scan_indel:
        is_scan_indel = {'yes': False}

    base, strand, indel = [], [], []
    # skip read which mapping quality less then 30
    for read in [al for al in sample_pos_line.pileups if al.alignment.mapq >= 10]:

        if is_scan_indel['yes'] and read.indel:
            indel.append(bam.scan_indel(read,
                                        is_scan_indel['pos'],
                                        is_scan_indel['fa']))

        if not read.is_del and not read.is_refskip:
            # skip the base which base_quality < 20
            if read.alignment.query_qualities[read.query_position] < 20:
                continue

            base.append(read.alignment.query_sequence[read.query_position])
            strand.append('-' if read.alignment.is_reverse else '+')

    return base, strand, indel


class CvgSingleProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """

    def __init__(self, ref_file, aligne_files, out_cvg_file,
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
        with open(self.out_cvg_file, 'w') as CVG:

            CVG.write('\t'.join(['#CHROM', 'POS', 'REF', 'Depth'] + self.cmm.BASE +
                                ['Indel'] + self.samples) + '\n')

            # set header
            for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):
                # ``regions`` is a 2-D array : [[start1, end1], [start2, end2], ...]
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
                        iter_tokes.append(bf.pileup(chrid, start - 1, end))
                    except ValueError:
                        iter_tokes.append('')

                # get all the reference sequence of chrid
                fa = self.ref_file_hd.fetch(chrid)

                # Set iteration marker: 1->iterate; 0->do not
                # iterate or hit the end
                go_iter = [1] * len(iter_tokes)
                n = 0
                for start, end in regions:
                    for position in range(start, end + 1):

                        if n % 100000 == 0:
                            sys.stderr.write("[INFO] loading lines %d at position %s:%d\t%s\n" %
                                             (n, chrid, position, time.asctime()))

                        n += 1
                        sample_info = [utils.fetch_next(iter_tokes[i]) if g else sample_info[i]
                                       for i, g in enumerate(go_iter)]

                        # sample_base, sample_base_qual, strands, mapqs and
                        # read_pos_rank are listed the same orde with each other.
                        sample_base, strands, indels = fetch_base_by_position(
                            position - 1,  # postion in pysam is 0-base
                            sample_info,
                            go_iter,
                            iter_tokes,
                            fa,  # Fa sequence for indel sequence
                            is_scan_indel=True
                        )

                        ref_base = fa[position - 1]
                        # ignore positions if coverage=0 or ref base is 'N' base
                        if not sample_base or ref_base in ['N', 'n']:
                            continue

                        self._out_cvg_file(chrid, position, ref_base, sample_base, indels, CVG)

        self._close_aligne_file()

    def _out_cvg_file(self, chrid, position, ref_base, sample_base, indels, out_file_handle):

        """
        :param chrid:
        :param position:
        :param ref_base:
        :param sample_base: 2d array,
        :param indels: 2d array
        :return:
        """

        # coverage info for each position
        base_depth = {b: 0 for b in self.cmm.BASE}
        for k, bs in enumerate(sample_base):

            # ignore all bases('*') which not match ``cmm.BASE``
            for b in bs:
                if b in base_depth:
                    base_depth[b] += 1

        # deal with indels
        indel_dict = {}
        for sample_ind in indels:
            for ind in sample_ind:
                indel_dict[ind] = indel_dict.get(ind, 0) + 1

        indel_string = ','.join(
            [k + ':' + str(v) for k, v in indel_dict.items()]) if indel_dict else '.'

        if sum(base_depth.values()):
            out_file_handle.write('\t'.join(
                [chrid, str(position), ref_base, str(sum(base_depth.values()))] +
                [str(base_depth[b]) for b in self.cmm.BASE] + [indel_string] +
                map(str, [len(s) for s in sample_base])
            ) + '\n')

        return


###############################################################################
class CvgMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.
    """

    def __init__(self, ref_in_file, aligne_files, out_cvg_file,
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

        self.single_process = CvgSingleProcess(ref_in_file,
                                               aligne_files,
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
