"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``BaseVar.py`` is "Sauron", and 'executor.py' module could just be called by it.
"""
from __future__ import division

import os
import sys
import time

from pysam import AlignmentFile, TabixFile, tabix_index

from . import utils
from . import CallerProcess, process_runner

from .basetypebam import BaseVarProcess
from .batchgenerator import BatchProcess
from .basetypebatch import BaseVarBatchProcess
from .coverageprocess import CvgSingleProcess
from .vqsr import vqsr


def _generate_regions_for_each_process(regions, process_num=1):
    """create regions for each process"""

    regions_for_each_process = [[] for _ in range(process_num)]

    # calculate region size for each processes
    total_regions_size = sum([e - s + 1 for _, s, e in regions])
    delta = int(total_regions_size / process_num) + 1 if total_regions_size % process_num else \
        int(total_regions_size / process_num)

    # reset
    total_regions_size = 0
    for chrid, start, end in regions:

        sub_region_size = end - start + 1
        # find index of regions_for_each_process by the total_regions_size
        i = int(total_regions_size / delta)
        total_regions_size += sub_region_size

        pre_size = 0
        if len(regions_for_each_process[i]):
            pre_size = sum([e - s + 1 for _, s, e in regions_for_each_process[i]])

        if sub_region_size + pre_size > delta:
            d = delta - pre_size
            regions_for_each_process[i].append([chrid, start, start + d - 1])

            i += 1
            for k, pos in enumerate(range(start + d - 1, end, delta)):
                s = pos + 1 if pos + 1 < end else end
                e = pos + delta if pos + delta < end else end

                regions_for_each_process[i + k].append([chrid, s, e])

        else:
            regions_for_each_process[i].append([chrid, start, end])

    return regions_for_each_process


def _load_sample_id_from_bam(bamfiles, filename_has_samplename=True):
    """loading sample id in BAM/CRMA files from RG tag"""

    sys.stderr.write('[INFO] Start loading all samples\' id from alignment files\n')
    if filename_has_samplename:
        sys.stderr.write('[INFO] loading samples\' id from filename because you '
                         'set "--filename-has-samplename"\n')
    sample_id = []
    for i, al in enumerate(bamfiles):

        if i % 1000 == 0:
            sys.stderr.write("[INFO] loading %d/%d alignment files ... %s\n" %
                             (i + 1, len(bamfiles), time.asctime()))

        if filename_has_samplename:
            filename = os.path.basename(al)

            # sample id should be the first element separate by ".",
            # e.g: "CL100045504_L02_61.sorted.rmdup.realign.BQSR.bam", "CL100045504_L02_61" is sample id.
            sample_id.append(filename.split(".")[0])

        else:

            # This situation will take a very long time to get sampleID from BAM header.
            bf = AlignmentFile(al)
            if 'RG' not in bf.header:
                sys.stderr.write('[ERROR] Bam file format error: missing @RG in the header.\n')
                bf.close()
                sys.exit(1)

            sample_id.append(bf.header['RG'][0]['SM'])
            bf.close()

    sys.stderr.write('[INFO] Finish loading all %d samples\' ID\n\n' % len(sample_id))
    return sample_id


def _output_cvg_and_vcf(sub_cvg_files, sub_vcf_files, outcvg, outvcf=None):
    """CVG file and VCF file could use the same tabix strategy."""
    for out_final_file, sub_file_list in zip([outcvg, outvcf], [sub_cvg_files, sub_vcf_files]):

        if out_final_file:
            _output_file(sub_file_list, out_final_file, del_raw_file=True)

    return


def _output_file(sub_files, out_file_name, del_raw_file=False):
    if out_file_name.endswith(".gz"):
        utils.merge_files(sub_files, out_file_name, output_isbgz=True, is_del_raw_file=del_raw_file)

        # Column indices are 0-based. Note: this is different from the tabix command line
        # utility where column indices start at 1.
        tabix_index(out_file_name, force=True, seq_col=0, start_col=1, end_col=1)
    else:
        utils.merge_files(sub_files, out_file_name, is_del_raw_file=del_raw_file)

    return


class BaseTypeRunner(object):

    def __init__(self, args):
        """init function
        """
        # setting parameters
        self.nCPU = args.nCPU
        self.mapq = args.mapq
        self.reference_file = args.referencefile
        self.pop_group_file = args.pop_group_file
        self.batchcount = args.batchcount

        self.outvcf = args.outvcf if args.outvcf else None
        self.outcvg = args.outcvg if args.outcvg else None
        self.outbatchfile = args.outbatchfile if args.outbatchfile else None

        self.smartrerun = True if args.smartrerun else False

        # Loading positions or load all the genome regions
        regions = utils.load_target_position(self.reference_file, args.positions, args.regions)
        self.regions_for_each_process = _generate_regions_for_each_process(regions, process_num=self.nCPU)

        self.alignfiles = args.input
        if args.infilelist:
            self.alignfiles += utils.load_file_list(args.infilelist)

        # setting the resolution of MAF
        self.min_af = utils.set_minaf(len(self.alignfiles)) if (args.min_af is None) else args.min_af

        sys.stderr.write('[INFO] Finish loading arguments and we have %d BAM/CRAM files for '
                         'variants calling. %s\n' % (len(self.alignfiles), time.asctime()))

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.sample_id = _load_sample_id_from_bam(self.alignfiles, filename_has_samplename=args.filename_has_samplename)

    def batch_generator(self):
        """Create batchfile for the input align files"""
        sys.stderr.write('[INFO] Start create batch file ... %s\n' % time.asctime())

        out_batch_names = []
        processes = []
        # Always create process manager even if nCPU==1, so that we can
        # listen signals from main thread
        for i in range(self.nCPU):
            sub_batch_file = self.outbatchfile + '_temp_%s' % i
            out_batch_names.append(sub_batch_file)

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s]\n' % (i + 1, self.nCPU, sub_batch_file))

            processes.append(CallerProcess(BatchProcess,
                                           self.reference_file,
                                           self.alignfiles,
                                           self.regions_for_each_process[i],
                                           self.sample_id,
                                           mapq=self.mapq,
                                           batchcount=self.batchcount,
                                           out_batch_file=sub_batch_file,
                                           rerun=self.smartrerun))

        process_runner(processes)

        # Final output
        _output_file(out_batch_names, self.outbatchfile, del_raw_file=True)

        return processes

    def basevar_caller(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call variants by BaseType ... %s\n' % time.asctime())

        out_vcf_names = []
        out_cvg_names = []

        processes = []
        # Always create process manager even if nCPU==1, so that we can
        # listen signals from main thread

        for i in range(self.nCPU):
            sub_cvg_file = self.outcvg + '_temp_%s' % i
            out_cvg_names.append(sub_cvg_file)

            if self.outvcf:
                sub_vcf_file = self.outvcf + '_temp_%s' % i
                out_vcf_names.append(sub_vcf_file)
            else:
                sub_vcf_file = None

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s, %s]\n' % (i + 1, self.nCPU, sub_vcf_file,
                                             sub_cvg_file))

            processes.append(CallerProcess(BaseVarProcess,
                                           self.reference_file,
                                           self.alignfiles,
                                           self.pop_group_file,
                                           self.regions_for_each_process[i],
                                           self.sample_id,
                                           min_af=self.min_af,
                                           mapq=self.mapq,
                                           batchcount=self.batchcount,
                                           out_cvg_file=sub_cvg_file,
                                           out_vcf_file=sub_vcf_file,
                                           rerun=self.smartrerun))

        process_runner(processes)

        # Final output
        _output_cvg_and_vcf(out_cvg_names, out_vcf_names, self.outcvg, outvcf=self.outvcf)

        return processes


class BaseTypeBatchRunner(object):

    def __init__(self, args):
        """init function
        """
        # setting parameters

        self.reference_file = args.referencefile
        self.outcvg = args.outcvg
        self.outvcf = args.outvcf if args.outvcf else None

        self.nCPU = args.nCPU
        self.pop_group_file = args.pop_group_file

        # Loading positions or load all the genome regions
        regions = utils.load_target_position(self.reference_file, args.positions, args.regions)
        self.regions_for_each_process = _generate_regions_for_each_process(regions, process_num=self.nCPU)

        self.batch_files = args.input
        if args.infilelist:
            self.batch_files += utils.load_file_list(args.infilelist)

        # loading all the sample id from batch_files
        self.sample_id = self._load_sample_id_from_batchheader()

        # setting the resolution of MAF
        self.min_af = utils.set_minaf(len(self.sample_id)) if (args.min_af is None) else args.min_af
        sys.stderr.write('[INFO] Finish loading arguments and we have %d BAM/CRAM files for '
                         'variants calling. %s\n' % (len(self.batch_files), time.asctime()))

    def _load_sample_id_from_batchheader(self):
        """loading sample id from header of batchfile"""

        sys.stderr.write('[INFO] Start loading all samples\' id from header of batch files\n')
        sample_id = []
        for i, al in enumerate(self.batch_files):

            if i % 1000 == 0:
                sys.stderr.write("[INFO] loading %d/%d batch files ... %s\n" %
                                 (i + 1, len(self.batch_files), time.asctime()))

            bf = TabixFile(al)
            get_sample_id = False
            for h in bf.header:
                if h.startswith("##SampleIDs"):
                    get_sample_id = True
                    for sample in h.split("=")[-1].split(","):
                        sample_id.append(sample)

                    break
            bf.close()

            if not get_sample_id:
                sys.stderr.write('[ERROR] batch file format error: missing ##SampleIDs in %s header.\n' % al)
                sys.exit(1)

        sys.stderr.write('[INFO] Finish load all %d samples\' ID '
                         'from file header\n\n' % len(sample_id))
        return sample_id

    def basevar_caller(self):
        """Run variant caller"""
        sys.stderr.write('[INFO] Start call variants by BaseType ... %s\n' % time.asctime())

        out_vcf_names = []
        out_cvg_names = []

        processes = []
        # Always create process manager even if nCPU==1, so that we can
        # listen signals from main thread
        for i in range(self.nCPU):

            sub_cvg_file = self.outcvg + '_temp_%s' % i
            out_cvg_names.append(sub_cvg_file)

            if self.outvcf:
                sub_vcf_file = self.outvcf + '_temp_%s' % i
                out_vcf_names.append(sub_vcf_file)
            else:
                sub_vcf_file = None

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s, %s]\n' % (i + 1, self.nCPU, sub_vcf_file,
                                             sub_cvg_file))

            processes.append(CallerProcess(BaseVarBatchProcess,
                                           self.reference_file,
                                           self.batch_files,
                                           self.pop_group_file,
                                           self.regions_for_each_process[i],
                                           self.sample_id,
                                           min_af=self.min_af,
                                           out_cvg_file=sub_cvg_file,
                                           out_vcf_file=sub_vcf_file))

        process_runner(processes)

        # Final output
        _output_cvg_and_vcf(out_cvg_names, out_vcf_names, self.outcvg, outvcf=self.outvcf)

        return processes


class VQSRRuner(object):
    """Runner for VQSR"""

    def __init__(self, args):
        """Init function"""
        self.args = args

    def run(self):
        vqsr.main(self.args)
        return


class CoverageRunner(object):

    def __init__(self, args, cmm=utils.CommonParameter()):
        """init function
        """
        self.referencefile = args.referencefile
        self.nCPU = args.nCPU
        self.outputfile = args.outputfile
        self.cmm = cmm

        # Loading positions if not provid we'll load all the genome
        self.regions = utils.load_target_position(self.referencefile, args.positions, args.regions)

        # Get all the input alignement files
        self.alignefiles = utils.load_file_list(args.infilelist)
        sys.stderr.write('[INFO] Finish loading parameters and input file list %s\n' % time.asctime())

    def run(self):
        """
        Run variant caller
        """
        sys.stderr.write('[INFO] Start call varaintis by BaseType ... %s\n' % time.asctime())

        # Always create process manager even if nCPU==1, so that we can
        # listen for signals from main thread
        regions_for_each_process = [[] for _ in range(self.nCPU)]
        if len(self.regions) < self.nCPU:
            # We cut the region evenly to fit nCPU if regions < nCPU
            for chrid, start, end in self.regions:
                delta = int((end - start + 1) / self.nCPU)
                if delta == 0:
                    delta = 1

                for i, pos in enumerate(range(start - 1, end, delta)):
                    s = pos + 1 if pos + 1 < end else end
                    e = pos + delta if pos + delta < end else end

                    regions_for_each_process[i % self.nCPU].append([chrid, s, e])

        else:
            for i, region in enumerate(self.regions):
                regions_for_each_process[i % self.nCPU].append(region)

        out_cvg_names = set()
        processes = []
        for i in range(self.nCPU):
            sub_cvg_file = self.outputfile + '_temp_%s' % i

            out_cvg_names.add(sub_cvg_file)
            processes.append(CallerProcess(CvgSingleProcess,
                                           self.referencefile,
                                           self.alignefiles,
                                           sub_cvg_file,
                                           regions_for_each_process[i],
                                           cmm=self.cmm))

        process_runner(processes)

        # Final output file name
        utils.merge_files(out_cvg_names, self.outputfile, is_del_raw_file=True)

        return processes


class MergeRunner(object):
    """Runner for merging files"""

    def __init__(self, args):
        """init function"""

        self.inputfiles = args.input
        if args.infilelist:
            self.inputfiles += utils.load_file_list(args.infilelist)

        self.outputfile = args.outputfile

    def run(self):
        _output_file(self.inputfiles, self.outputfile)
        return


class NearbyIndelRunner(object):
    """Add Nearby Indel density and type information for each variants of VCF"""

    def __init__(self, args):
        """init function"""
        args.nearby_dis_around_indel = int(args.nearby_dis_around_indel)
        self.in_vcf_file = args.in_vcf_file
        self.in_cvg_file = args.in_cvg_file
        self.nearby_dis_around_indel = args.nearby_dis_around_indel

        sys.stderr.write('[INFO] basevar NearbyIndel'
                         '\n\t-i %s'
                         '\n\t-c %s'
                         '\n\t-d %d' % (args.in_vcf_file,
                                        args.in_cvg_file,
                                        args.nearby_dis_around_indel)
                         )

    def run(self):
        from .other import NearbyIndel

        nbi = NearbyIndel(self.in_vcf_file, self.in_cvg_file, nearby_distance=self.nearby_dis_around_indel)
        nbi.run()

        return self


class PopulationMatrixRunner(object):
    """Create Population matrix"""
    def __init__(self, args):
        # setting parameters
        self.alignfiles = args.input
        self.input_postion_file = args.positions
        self.output = args.output_file

        self.nCPU = args.nCPU
        self.mapq = args.mapq
        self.reference_file = args.referencefile
        self.batchcount = args.batchcount
        self.smartrerun = True if args.smartrerun else False

        # Loading positions or load all the genome regions
        sites, regions = self.load_position()
        self.sites = sites  # A dict
        self.regions_for_each_process = _generate_regions_for_each_process(regions, process_num=self.nCPU)

        if args.infilelist:
            self.alignfiles += utils.load_file_list(args.infilelist)

        sys.stderr.write('[INFO] Finish loading arguments and we have %d BAM/CRAM files for '
                         'variants calling. %s\n' % (len(self.alignfiles), time.asctime()))

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.sample_id = _load_sample_id_from_bam(self.alignfiles, filename_has_samplename=args.filename_has_samplename)

    def load_position(self):
        sites, regiondict = {}, {}
        with utils.Open(self.input_postion_file, 'rb') as f:
            for r in f:
                # chr1	15777	G	A
                tok = r.strip().split()
                k = tok[0] + ':' + tok[1]
                sites[k] = [tok[2].upper(), tok[3].upper()]

                if tok[0] not in regiondict:
                    regiondict[tok[0]] = []

                regiondict[tok[0]].append([int(tok[1]), int(tok[1])])

        # merge and sort the regions make sure all the positions are in order
        regions = []
        for chrid, v in sorted(regiondict.items(), key=lambda x: x[0]):
            for start, end in utils.merge_region(v):
                regions.append([chrid, start, end])

        return sites, regions

    def create_matrix(self):
        """Create matrix from the input align files"""
        sys.stderr.write('[INFO] Start create population matrix ... %s\n' % time.asctime())

        out_basebatch_names = []
        processes = []
        # Always create process manager even if nCPU==1, so that we can
        # listen signals from main thread
        for i in range(self.nCPU):
            sub_basebatch_file = self.output + '_temp_%s' % i
            out_basebatch_names.append(sub_basebatch_file)

            sys.stderr.write('[INFO] Process %d/%d output to temporary files:'
                             '[%s]\n' % (i + 1, self.nCPU, sub_basebatch_file))
            processes.append(CallerProcess(BatchProcess,
                                           self.reference_file,
                                           self.alignfiles,
                                           self.regions_for_each_process[i],
                                           self.sample_id,
                                           mapq=self.mapq,
                                           batchcount=self.batchcount,
                                           out_batch_file=sub_basebatch_file,
                                           rerun=self.smartrerun,
                                           justbase=True))

        process_runner(processes)

        # Final output: regions_for_each_process keep the same position order with the input file 'args.positions'
        with utils.Open(self.output, 'wb', isbgz=True if self.output.endswith(".gz") else False) as OUT:
            for index, sub_file_name in enumerate(out_basebatch_names):
                sample_id = []
                with utils.Open(sub_file_name, 'rb') as IN:
                    for line in IN:
                        col = line.strip().split()
                        if line[0] == "#":
                            if index == 0:
                                if col[0].startswith("##SampleIDs="):
                                    sample_id = col[0].split("=")[1].split(",")

                                if col[0].startswith("#CHROM"):
                                    col.pop()
                                    OUT.write("#CHROM\tPOS\tREF\tALT\tDepth\tALT_Freq\t%s\n" % "\t".join(sample_id))
                        else:
                            # #CHROM  POS  REF  Depth(CoveredSample)  Readbases
                            bases = col[4].split(",")

                            k = col[0] + ":" + col[1]
                            (ref_base, target_alt) = self.sites[k] if k in self.sites else ('', '')
                            if ref_base and ref_base != col[2].upper():
                                sys.stderr.write(
                                    "[Error] Error happen when final output in create_matrix(), "
                                    "reference base (%s != %s) in %s, in %s and %s!\n" % (
                                        col[2], ref_base, k, sub_file_name, self.input_postion_file)
                                )
                                sys.exit(1)

                            depth = float(col[3])
                            target_alt_num = 0
                            for i in range(len(bases)):
                                if bases[i] != "N":
                                    if bases[i] == target_alt:
                                        target_alt_num += 1
                                        bases[i] = '1'
                                    elif bases[i] == ref_base:
                                        bases[i] = '0'
                                    else:
                                        bases[i] = '.'
                                else:
                                    bases[i] = '.'

                            OUT.write("%s\t%.6f\t%s\n" % ("\t".join([col[0], col[1], col[2], target_alt, col[3]]),
                                                          target_alt_num/depth if depth > 0 else 0,
                                                          "\t".join(bases)))

                # remove the tmp file
                os.remove(sub_file_name)

        if self.output.endswith(".gz"):
            tabix_index(self.output, force=True, seq_col=0, start_col=1, end_col=1)
