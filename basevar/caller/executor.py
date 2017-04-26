"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``BaseVar.py`` is "Sauron", and this module could just be called by it.
"""
from __future__ import division

import sys
import os
import argparse
import heapq
import time

from . import utils
from .variantcaller import BaseVarMultiProcess


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


def merge_files(temp_file_names, final_file_name):
    """
    Merging output VCF/CVG files into a final file
    log.info("Merging output VCF/CVG file(s) into final file %s" %(final_file_name))
    """

    # Final output file
    if final_file_name == "-":
        output_file = sys.stdout
    else:
        output_file = utils.Open(final_file_name, 'wb')
    the_heap = []

    # Initialise queue
    for index, file_name in enumerate(temp_file_names):
        the_file = utils.Open(file_name, 'rb')

        for line in the_file:

            # End of this file
            if line[0] == "#":
                if index == 0:
                    output_file.write(line)
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
        output_file.write(next_file.line)

        # Put file back on heap
        try:
            next_file.next()
            heapq.heappush(the_heap, next_file)
        except StopIteration:
            continue

    # Close final output file
    if final_file_name != "-":
        output_file.close()

    #log.info("Finished merging %s file(s)"%final_file_name)
    return


class Runner(object):

    def __init__(self, cmm=utils.CommonParameter()):
        """init function
        """
        optp = argparse.ArgumentParser()
        optp.add_argument('basetype')
        optp.add_argument('-o', '--outprefix', dest='outprefix', metavar='FILE', default='out',
                          help='The prefix of output files. [out]')
        optp.add_argument('-l', '--mpileup-list', dest='infilelist', metavar='FILE',
                          help='The input mpileup file list.', default='')
        optp.add_argument('-L', '--positions', metavar='FILE', dest='positions',
                          help='skip unlisted positions (chr pos)', default='')
        optp.add_argument('-R', '--regions', metavar='chr:start-end', dest='regions',
                          help='skip positions not in (chr:start-end)', default='')
        optp.add_argument('-s', '--sample-list', dest='samplelistfile', metavar='FILE',
                          help='The sample list.')
        optp.add_argument('-S', '--subsample-list', dest='subsample', metavar='FILE',
                          help='Skip samples not in subsample-list, one sample per row.')
        optp.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int,
                          help='Number of processer to use. [1]', default=1)
        optp.add_argument('-m', '--min_af', dest='min_af', type=float, metavar='MINAF', default=0.001,
                          help='By setting min AF to skip uneffective caller positions to accelerate program speed. [0.001]')

        opt = optp.parse_args()
        self.opt = opt

        self.cmm = cmm
        # reset threshold of init min allele frequence by read depth
        self.cmm.MINAF = self.opt.min_af

        if len(sys.argv) == 2 and len(opt.infilelist) == 0:
            optp.error('[ERROR] At least input one mpileup file\n')

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
        # [[chrid1, start1, end1], [chrid2, start2, end2], ...]
        self.regions = []
        for chrid, v in sorted(_sites.items(), key=lambda x: x[0]):
            for start, end in utils.merge_region(v):
                self.regions.append([chrid, start, end])

        # Load all the mpileup files
        self.mpileupfiles = [f for f in sys.argv if '.mpileup.gz' in f]
        if opt.infilelist:
            self.mpileupfiles.extend(utils.load_file_list(opt.infilelist))

    def basetype(self):
        """
        Run variant caller
        """

        # Always create process manager even if nCPU==1, so that we can
        # listen for signals from main thread
        processes = []
        regions_for_each_process = [[] for _ in range(self.opt.nCPU)]
        for i, region in enumerate(self.regions):
            regions_for_each_process[i % self.opt.nCPU].append(region)

        out_vcf_names = set()
        out_cvg_names = set()

        for i in range(self.opt.nCPU):
            sub_vcf_file = self.opt.outprefix + '_temp_%s'%i + '.vcf'
            sub_cvg_file = self.opt.outprefix + '_temp_%s'%i + '.cvg.tsv'

            out_vcf_names.add(sub_vcf_file)
            out_cvg_names.add(sub_cvg_file)
            processes.append(BaseVarMultiProcess(self.mpileupfiles,
                                                 sub_vcf_file,
                                                 sub_cvg_file,
                                                 regions_for_each_process[i],
                                                 self.opt,
                                                 cmm=self.cmm))

        for p in processes:
            p.start()

        # listen for signal while any process is alive
        while True in [p.is_alive() for p in processes]:
            try:
                time.sleep(1)

            except KeyboardInterrupt:
                print 'KeyboardInterrupt detected, terminating all processes...'
                for p in processes:
                    p.terminate()

                sys.exit(1)

        # make sure all process are finished
        for p in processes:
            p.join()

        # Final output file name
        out_vcf_file = self.opt.outprefix + '.vcf'
        out_cvg_file = self.opt.outprefix + '.cvg.tsv' # position coverage

        merge_files(out_vcf_names, out_vcf_file)
        merge_files(out_cvg_names, out_cvg_file)

        return



