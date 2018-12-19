"""
This is a Process module for BaseType by BAM/CRAM

"""
import os
import sys
import time
import multiprocessing

import pysam

from . import utils
from . import bam
from .basetypeprocess import basetypeprocess


class BaseVarSingleProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """
    def __init__(self, ref_file, aligne_files, in_popgroup_file, regions, samples, batchcount=1000,
                 out_vcf_file=None, out_cvg_file=None, rerun=False, cmm=None):
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
        self.batchcount = batchcount
        self.samples = samples
        self.smartrerun = rerun
        self.cmm = cmm
        self.regions = {}

        # store the region into a dict
        for chrid, start, end in regions:

            if chrid not in self.regions:
                self.regions[chrid] = []

            self.regions[chrid].append([start, end])

        # loading population group
        # group_id => [a list samples_index]
        self.popgroup = {}
        if in_popgroup_file and len(in_popgroup_file):
            self.popgroup = utils.load_popgroup_info(self.samples, in_popgroup_file)

        # create batch file for variant discovery

    def _open_aligne_files(self, bamfiles):

        ali_files_hd = []
        for f in bamfiles:

            try:
                bf = pysam.AlignmentFile(f)

            except ValueError:
                sys.stderr.write('[ERROR] Input file: %s is not BAM nor CRAM.\n' % f)
                self._close_aligne_file(ali_files_hd)
                sys.exit(1)

            ali_files_hd.append(bf)

        return ali_files_hd

    def _close_aligne_file(self, ali_files_hd):

        for bf in ali_files_hd:
            bf.close()

        return

    def create_batch_file_for_region(self, chrid, regions, batchcount):
        """
        ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
        """
        # Create a batch of temp files for variant discovery
        # ``self.out_vcf_file`` may be None if just calculating position depth
        batchfile_dir = utils.safe_makedir(os.path.split(os.path.realpath(self.out_cvg_file))[0] + "/batchfiles")

        # store all the batch files
        batchfiles = []
        part_num = len(self.aligne_files) / batchcount
        if part_num * batchcount < len(self.aligne_files):
            part_num += 1

        tmp_region = []
        for p in regions:
            tmp_region.extend(p)

        tmp_region = sorted(tmp_region)
        bigstart, bigend = tmp_region[0], tmp_region[-1]

        m = 0
        for i in range(0, len(self.aligne_files), batchcount):

            m += 1
            part_file_name = "%s/BaseVar.%s.%d_%d.batch" % (batchfile_dir,
                                                            ".".join(map(str, [chrid, bigstart, bigend])),
                                                            m,
                                                            part_num)
            # store the name of batchfile into a list.
            batchfiles.append(part_file_name)

            if self.smartrerun and os.path.isfile(part_file_name):
                # ``part_file_name`` exists and we do not have create it again if setting `smartrerun`
                sys.stderr.write("[INFO] %s has been created we don't need to create again, "
                                 "when you set `smartrerun` %s\n" % (part_file_name, time.asctime()))
                continue
            else:
                sys.stderr.write("[INFO] Creating batchfile %s at %s\n" % (part_file_name, time.asctime()))

            # One batch of alignment files
            sub_alignfiles = self.aligne_files[i:i+batchcount]
            ali_files_hd = self._open_aligne_files(sub_alignfiles)

            # ``iter_tokes`` is a list of iterator for each sample's input file
            iter_tokes = []
            for i, bf in enumerate(ali_files_hd):
                try:
                    # 0-base
                    iter_tokes.append(bf.pileup(chrid, bigstart-1, bigend))
                except ValueError:
                    if self.cmm.debug:
                        sys.stderr.write("# [WARMING] Empty region %s:%d-%d in %s" %
                                         (chrid, bigstart-1, bigend, self.aligne_files[i]))
                    iter_tokes.append("")

            with open(part_file_name, "w") as OUT:
                OUT.write("%s\n" % "\t".join(["#CHROM", "POS", "REF", "MappingQuality", "Readbases", "ReadbasesQuality",
                                               "ReadPositionRank", "Strand"]))

                # Set iteration marker: 1->iterate; 0->Do not iterate or hit the end
                sample_info = [utils.fetch_next(it) for it in iter_tokes]

                # get sequence of chrid from reference fasta
                fa = self.ref_file_hd.fetch(chrid)

                n = 0
                for start, end in regions:

                    sys.stderr.write('[INFO] Fetching info from %d samples in region %s'
                                     ' at %s\n' % (len(self.aligne_files),
                                                   chrid + ":" + str(start) + "-" + str(end),
                                                   time.asctime()))

                    for position in range(start, end + 1):

                        if n % 100000 == 0:
                            sys.stderr.write("[INFO] loading lines %d at position %s:%d\t%s\n" %
                                             (n + 1, chrid, position, time.asctime()))
                        n += 1

                        ref_base = fa[position - 1]
                        # ref base is 'N' base, very important
                        if ref_base.upper() not in ['A', 'C', 'G', 'T']:
                            continue

                        sample_bases, sample_base_quals, strands, mapqs, read_pos_rank = bam.fetch_base_by_position(
                            position - 1,  # postion for pysam is 0-base
                            sample_info,
                            iter_tokes,
                            fa  # Fa sequence for indel sequence
                        )

                        OUT.write("%s\n" % "\t".join([
                            chrid,
                            str(position),
                            ref_base,
                            ",".join(map(str, mapqs)),
                            ",".join(sample_bases),
                            ",".join(map(str, sample_base_quals)),
                            ",".join(map(str, read_pos_rank)),
                            ",".join(strands)
                        ]))

            self._close_aligne_file(ali_files_hd)
            sys.stderr.write("[INFO] Done for batchfile %s at %s\n\n" % (part_file_name, time.asctime()))

        return batchfiles

    def fetch_baseinfo_by_position(self, chrid, position, ref_base, infolines):

        sample_bases, sample_base_quals, mapqs, read_pos_rank, strands = [], [], [], [], []
        for i, line in enumerate(infolines):
            # <#CHROM  POS REF MappingQuality  Readbases  ReadbasesQuality  ReadPositionRank, Strand>
            if len(line) == 0:
                sys.stderr.write("[Error] %d lines happen to be empty in batchfiles!\n" % (i+1))
                sys.exit(1)

            if line.startswith("#"):
                continue

            col = line.strip().split()
            col[1] = int(col[1])
            if col[0] != chrid or col[1] != position or col[2] != ref_base:
                sys.stderr.write("[Error] %d lines, chromosome [%s and %s] or position [%d and %d] "
                                 "or ref-base [%s and %s] in batchfiles not match with each other!\n" %
                                 (i+1, col[0], chrid, col[1], position, col[2], ref_base))
                sys.exit(1)

            mapqs.append(col[3])
            sample_bases.append(col[4].upper())
            sample_base_quals.append(col[5])
            read_pos_rank.append(col[6])
            strands.append(col[7])

        # cat all the info together and create ...
        mapqs = map(int, ",".join(mapqs).split(","))
        sample_bases = ",".join(sample_bases).split(",")
        sample_base_quals = map(int, ",".join(sample_base_quals).split(","))
        read_pos_rank = map(int, ",".join(read_pos_rank).split(","))
        strands = ",".join(strands).split(",")

        return sample_bases, sample_base_quals, strands, mapqs, read_pos_rank

    def run(self):
        """
        Run the process of calling variant and output files.
        """
        vcf_header = utils.vcf_header_define()
        group = []  # Just for the header of CVG file
        if self.popgroup:
            for g in self.popgroup.keys():
                g_id = g.split('_')[0]  # ignore '_AF'
                group.append(g_id)
                vcf_header.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="Allele frequency in the %s '
                                  'populations calculated base on LRT, in the range (0,1)">' % (g_id, g_id))

        CVG = open(self.out_cvg_file, 'w')
        CVG.write('##fileformat=CVGv1.0\n')
        CVG.write('##Group information is the depth of A:C:G:T:Indel\n')
        CVG.write('\t'.join(['#CHROM', 'POS', 'REF', 'Depth'] + self.cmm.BASE +
                            ['Indel', 'FS', 'SOR', 'Strand_Coverage(REF_FWD,'
                             'REF_REV,ALT_FWD,ALT_REV)\t%s\n' % '\t'.join(group)]))

        VCF = open(self.out_vcf_file, 'w') if self.out_vcf_file else None
        if VCF:  # set header if VCF is not None
            VCF.write('\n'.join(vcf_header) + '\n')
            VCF.write('\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
                                 'INFO\tFORMAT'] + self.samples) + '\n')

        n = 0
        for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):

            batchfiles = self.create_batch_file_for_region(chrid, regions, self.batchcount)
            batch_files_hd = [open(f) for f in batchfiles]
            _ = [utils.fetch_next(f) for f in batch_files_hd]  # ignore the header

            # get sequence of chrid from reference fasta
            fa = self.ref_file_hd.fetch(chrid)
            for start, end in regions:
                for position in range(start, end + 1):

                    if n % 100000 == 0:
                        sys.stderr.write("[INFO] Have been loading %d lines when hit position %s:%d\t%s\n" %
                                         (n+1, chrid, position, time.asctime()))
                    n += 1

                    ref_base = fa[position - 1]
                    # ignore while ref base is 'N' base, very important
                    if ref_base.upper() not in ['A', 'C', 'G', 'T']:
                        continue

                    # The order of position in all the batchfiles will be the same and match exactly
                    # with `(start, end in regions)`
                    info = [utils.fetch_next(it) for it in batch_files_hd]

                    (sample_bases,
                     sample_base_quals,
                     strands,
                     mapqs,
                     read_pos_rank) = self.fetch_baseinfo_by_position(chrid, position, ref_base, info)

                    # # ignore positions if coverage=0
                    if sum(read_pos_rank) == 0:
                        continue

                    # Calling varaints by Basetypes and output VCF and Coverage files.
                    basetypeprocess(chrid,
                                    position,
                                    ref_base,
                                    sample_bases,
                                    sample_base_quals,
                                    mapqs,
                                    strands,
                                    read_pos_rank,
                                    self.popgroup,
                                    self.cmm,
                                    CVG,
                                    VCF)

            for fh in batch_files_hd:
                fh.close()

        CVG.close()
        if VCF: VCF.close()
        self.ref_file_hd.close()


###############################################################################
class BaseVarMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.
    """
    def __init__(self, ref_in_file, aligne_files, pop_group_file, regions, samples_id, batchcount=1000,
                 out_vcf_file=None, out_cvg_file=None, rerun=False, cmm=None):
        """
        Constructor.

        regions: 2d-array like, required
                It's region info , format like: [[chrid, start, end], ...]
        """
        multiprocessing.Process.__init__(self)

        # loading all the sample id from aligne_files
        # ``samples_id`` has the same size and order as ``aligne_files``
        self.single_process = BaseVarSingleProcess(ref_in_file,
                                                   aligne_files,
                                                   pop_group_file,
                                                   regions,
                                                   samples_id,
                                                   batchcount=batchcount,
                                                   out_cvg_file=out_cvg_file,
                                                   out_vcf_file=out_vcf_file,
                                                   rerun=rerun,
                                                   cmm=cmm)

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()
