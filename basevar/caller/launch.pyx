# cython: profile=True
"""
This module will contain all the executor steps of BaseVar.

We have many important modules in BaseVar while this one is
the lord to rule them all, in a word, it's "The Ring".

``runner.py`` is "Sauron", and 'launch.pyx' module could just be called by it.
"""
import os
import sys

from libc.time cimport time, time_t, localtime, asctime
from libc.stdio cimport fprintf, stderr, stdout
from libc.stdlib cimport calloc, free
from libc.stdlib cimport exit, EXIT_FAILURE

##############################################################################
from basevar import utils
from basevar.utils cimport generate_regions_by_process_num, fast_merge_files, load_popgroup_info
from basevar.io.openfile import Open

from basevar.datatype.strarray cimport strarray_init_by_pylist, strarray_destroy, strarray_append, \
    convert_strarray_to_list

from basevar.caller.do import CallerProcess, process_runner
from basevar.caller.basetypeprocess cimport BaseVarProcess

from basevar.io.BGZF.tabix import tabix_index
from basevar.io.bam cimport get_sample_names
from basevar.caller.vqsr import vqsr
from basevar.caller.other import NearbyIndel

cdef class BaseTypeRunner(object):
    def __cinit__(self, args):
        """init function
        """
        # setting parameters
        cdef unsigned i = 0
        strarray_init_by_pylist(&self.align_files, args.input)

        if args.infilelist:
            with Open(args.infilelist, "rb") as fh:
                for line in fh:
                    if line[0] != '#':
                        strarray_append(&self.align_files, line.strip().split()[0])

        self.nCPU = args.nCPU
        self.reference_file = args.referencefile
        self.outcvg = args.outcvg
        self.outvcf = args.outvcf if args.outvcf else None
        self._set_cmd_options(args)

        # setting the resolution of MAF
        self.options.min_af = utils.set_minaf(self.align_files.size) if (args.min_af is None) else args.min_af
        fprintf(stdout, "\n[INFO] Finish loading arguments and we have %d BAM/CRAM files for "
                        "variants calling.\n", self.align_files.size)

        # Loading positions and sorted(and merge) if not been provided we'll load all the genome
        sorted_regions = utils.load_target_position(self.reference_file, args.positions, args.regions)
        self.regions_for_each_process = generate_regions_by_process_num(
            sorted_regions, process_num=self.nCPU, convert_to_2d=False)

        # ``samples_id`` has the same size and order as ``aligne_files``
        self.sample_ids = get_sample_names(self.align_files, True if args.filename_has_samplename else False)

        if self.options.batch_count > self.sample_ids.size:
            fprintf(stdout, "[WARNING] --batch-count (%d) is bigger than the number of alignment files (%ld). Reset "
                            "batch count to be %ld\n", self.options.batch_count, self.sample_ids.size, self.sample_ids.size)
            self.options.batch_count = self.sample_ids.size

        # loading population group
        # dict: group_id => [a list samples_index]
        if args.pop_group_file and len(args.pop_group_file):
            self.pop_group = load_popgroup_info(&self.sample_ids, args.pop_group_file)

    def __dealloc__(self):
        """Clean up memory.
        """
        strarray_destroy(&self.align_files)
        strarray_destroy(&self.sample_ids)

        return

    cdef void _set_cmd_options(self, args):
        # set option into C-Type struct

        self.options.mapq = args.mapq
        self.options.min_base_qual = args.min_base_qual
        self.options.r_len = args.r_len
        self.options.batch_count = args.batch_count
        self.options.smartrerun = args.smartrerun
        self.options.max_reads = args.max_reads
        self.options.is_compress_read = args.is_compress_read
        self.options.qual_bin_size = args.qual_bin_size
        self.options.trim_overlapping = args.trim_overlapping
        self.options.trim_soft_clipped = args.trim_soft_clipped
        self.options.filter_duplicates = args.filter_duplicates
        self.options.filter_reads_with_unmapped_mates = args.filter_reads_with_unmapped_mates
        self.options.filter_reads_with_distant_mates = args.filter_reads_with_distant_mates
        self.options.filter_read_pairs_with_small_inserts = args.filter_read_pairs_with_small_inserts
        self.options.verbosity = args.verbosity
        self.options.min_af = args.min_af

        return

    cdef bint basevar_caller(self):
        """
        Run variant caller
        """
        cdef time_t timer
        time(&timer)
        fprintf(stdout, '[INFO] Start calling variants by BaseType. %s\n', asctime(localtime(&timer)))

        cdef list out_vcf_names = []
        cdef list out_cvg_names = []
        cdef list successful_marker_files = []
        cdef list processes = []
        cdef int i = 0

        # I have to convert StringArray to be py list, because multiprocess in
        # CallerProcess could just accept Python Object.
        cdef list py_sample_ids_list = convert_strarray_to_list(&self.sample_ids)
        cdef list py_align_files_list = convert_strarray_to_list(&self.align_files)

        # Always create process manager even if nCPU==1, so that we can
        # listen signals from main thread
        for i in range(self.nCPU):
            sub_cvg_file = self.outcvg + '.temp_%d_%d' % (i+1, self.nCPU)
            out_cvg_names.append(sub_cvg_file)
            successful_marker_files.append(sub_cvg_file + ".PROCESS.AND_VCF_DONE_SUCCESSFULLY")

            if self.outvcf:
                sub_vcf_file = self.outvcf + '.temp_%d_%d' % (i + 1, self.nCPU)
                out_vcf_names.append(sub_vcf_file)
            else:
                sub_vcf_file = None

            print("[INFO] Process %d/%d output to temporary files:[%s, %s]." %
                  (i+1, self.nCPU, sub_vcf_file, sub_cvg_file))

            tmp_dir, name = os.path.split(os.path.realpath(sub_cvg_file))
            cache_dir = tmp_dir + "/Batchfiles.%s.WillBeDeletedWhenJobsFinish" % name

            if self.options.smartrerun and os.path.isfile(sub_cvg_file) and (not os.path.exists(cache_dir)):
                # if `cache_dir` is not exist and `sub_cvg_file` is exists means
                # `sub_cvg_file and sub_vcf_file` has been finished successfully.
                continue

            cache_dir = utils.safe_makedir(cache_dir)
            processes.append(CallerProcess(BaseVarProcess,
                                           py_sample_ids_list,
                                           py_align_files_list,
                                           self.reference_file,
                                           self.regions_for_each_process[i],
                                           self.pop_group,
                                           self.options,
                                           cache_dir,
                                           out_cvg_file=sub_cvg_file,
                                           out_vcf_file=sub_vcf_file))

        process_runner(processes)
        cdef bint all_process_success = True
        for p in processes:
            if p.exitcode != 0:
                all_process_success = False

        # double check for processes!
        cdef int *fail_process_num = <int *>calloc(self.nCPU, sizeof(int))
        for i, m in enumerate(successful_marker_files):
            if not os.path.exists(m):
                all_process_success = False
                fail_process_num[i] = i+1

            os.remove(m)  # remove any way

        # Final output if all the processes are ending successful!
        if all_process_success:
            utils.output_cvg_and_vcf(out_cvg_names, out_vcf_names, self.outcvg, outvcf=self.outvcf)
            fprintf(stdout, "[INFO] All the processes are done successfully.\n")
            free(fail_process_num)
        else:
            for i in range(self.nCPU):
                if fail_process_num[i] != 0:
                    fprintf(stderr, "[ERROR] The program is fail in [%d] processes. Abort!\n", fail_process_num[i])
            free(fail_process_num)
            exit(EXIT_FAILURE)

        return all_process_success

    # cdef bint basevar_caller_singleprocess(self):
    #     """
    #     Run variant caller --------- Just for Testting, when we done, please delete this function!!!!!!
    #     """
    #     cdef time_t timer
    #     time(&timer)
    #     fprintf(stdout, '[INFO] Start call variants by BaseType ... %s\n', asctime(localtime(&timer)))
    #
    #     out_vcf_names = []
    #     out_cvg_names = []
    #
    #     sub_cvg_file = self.outcvg + '_temp'
    #     out_cvg_names.append(sub_cvg_file)
    #
    #     if self.outvcf:
    #         sub_vcf_file = self.outvcf + '_temp'
    #         out_vcf_names.append(sub_vcf_file)
    #     else:
    #         sub_vcf_file = None
    #
    #     tmp_dir, name = os.path.split(os.path.realpath(sub_cvg_file))
    #     cache_dir = tmp_dir + "/Batchfiles.%s.WillBeDeletedWhenJobsFinish" % name
    #
    #     if self.options.smartrerun and os.path.isfile(sub_cvg_file) and (not os.path.exists(cache_dir)):
    #         # if `cache_dir` is not exist and `sub_cvg_file` is exists means
    #         # `sub_cvg_file and sub_vcf_file` has been finish successfully.
    #         return True
    #
    #     # cache_dir = utils.safe_makedir(cache_dir)
    #     # bp = BaseVarProcess(self.sample_id,
    #     #                     self.alignfiles,
    #     #                     self.reference_file,
    #     #                     self.regions_for_each_process[0],
    #     #                     out_cvg_file=sub_cvg_file,
    #     #                     out_vcf_file=sub_vcf_file,
    #     #                     cache_dir=cache_dir,
    #     #                     options=self.options)
    #     #
    #     # bp.run()
    #     #
    #     # # Final output
    #     # utils.output_cvg_and_vcf(out_cvg_names, out_vcf_names, self.outcvg, outvcf=self.outvcf)
    #     return True


class VQSRRunner(object):
    """Runner for VQSR"""
    def __init__(self, args):
        """Init function"""
        self.opt = args
        return

    def run(self):
        vqsr.run_VQSR(self.opt)
        return


class ApplyVQSRRunner(object):
    """Apply VQSR"""
    def __init__(self, args):
        self.opt = args
        return

    def run(self):
        vqsr.apply_VQSR(self.opt)


class MergeRunner(object):
    """Runner for merging files"""

    def __init__(self, args):
        """init function"""

        self.inputfiles = args.input
        if args.infilelist:
            self.inputfiles += utils.load_file_list(args.infilelist)

        self.outputfile = args.outputfile

    def run(self):
        # utils.output_file(self.inputfiles, self.outputfile)

        # A new merge methods!
        if self.outputfile.endswith(".gz"):
            fast_merge_files(self.inputfiles, self.outputfile, False)

            # Column indices are 0-based. Note: this is different from the tabix command line
            # utility where column indices start at 1.
            tabix_index(self.outputfile, force=True, seq_col=0, start_col=1, end_col=1)
        else:
            fast_merge_files(self.inputfiles, self.outputfile, False)

        return


class NearbyIndelRunner(object):
    """Add Nearby Indel density and type information for each variants of VCF"""

    def __init__(self, args):
        """init function"""
        args.nearby_dis_around_indel = int(args.nearby_dis_around_indel)
        self.in_vcf_file = args.in_vcf_file
        self.in_cvg_file = args.in_cvg_file
        self.output_file = args.outputfile
        self.nearby_dis_around_indel = args.nearby_dis_around_indel

        sys.stderr.write('[INFO] basevar NearbyIndel'
                         '\n\t-I %s'
                         '\n\t-C %s'
                         '\n\t-D %d'
                         '\n\t-O %s\n' % (args.in_vcf_file,
                                          args.in_cvg_file,
                                          args.nearby_dis_around_indel,
                                          args.outputfile))

    def run(self):
        nbi = NearbyIndel(self.in_vcf_file, self.in_cvg_file, self.output_file,
                          nearby_distance=self.nearby_dis_around_indel)
        nbi.run()
        return
