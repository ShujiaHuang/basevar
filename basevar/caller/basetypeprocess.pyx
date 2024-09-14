# cython: profile=True
"""
This is a Process module for BaseType by BAM/CRAM

"""
import os

from libc.stdio cimport fprintf, sprintf, stdout, stderr, remove
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.time cimport clock_t, clock, CLOCKS_PER_SEC

from basevar import utils
from basevar.utils cimport c_min, c_max
from basevar.utils cimport BaseTypeCmdOptions
from basevar.io.fasta import FastaFile

from basevar.datatype.strarray cimport StringArray, strarray_init, strarray_init_by_pylist, \
    strarray_append, strarray_destroy
from basevar.datatype.genomeregion cimport genome_region_init_by_pylist, genome_region_array_destroy

from basevar.caller.variantcaller import output_header
from basevar.caller.variantcaller cimport variants_discovery
from basevar.caller.batchcaller cimport generate_batchfile

cdef bint REMOVE_BATCH_FILE = True

cdef class BaseVarProcess:
    """
    simple class to repesent a single BaseVar process.
    """
    def __cinit__(self, list samples, list align_files, ref_file, regions, pop_group, BaseTypeCmdOptions options,
                  bytes cache_dir, out_cvg_file, out_vcf_file=None):
        """Constructor.

        Store input file, options and output file name.

        Parameters:
        ===========
            samples: list like
                A list of sample id

            regions: 2d-array like, required
                    It's region info , format like: [[chrid, start, end], ...]
        """
        strarray_init_by_pylist(&self.samples, samples)
        strarray_init_by_pylist(&self.align_files, align_files)

        self.fa_file_hd = FastaFile(ref_file, ref_file + ".fai")

        self.out_vcf_file = out_vcf_file
        self.out_cvg_file = out_cvg_file

        self.options = options
        self.cache_dir = cache_dir  # char*

        genome_region_init_by_pylist(&self.regions, regions)
        self.pop_group = pop_group

    def __dealloc__(self):
        """Clean up memory.
        """
        strarray_destroy(&self.samples)
        strarray_destroy(&self.align_files)
        genome_region_array_destroy(&self.regions)
        return

    def run(self):
        self.run_variant_discovery_by_batchfiles()
        return

    cdef void run_variant_discovery_by_batchfiles(self):

        VCF = open(self.out_vcf_file, "w") if self.out_vcf_file else None
        CVG = open(self.out_cvg_file, "w")

        _py_samples = []
        cdef size_t k = 0
        for k in range(self.samples.size):
            _py_samples.append(str(self.samples.array[k]))

        output_header(self.fa_file_hd.filename, _py_samples, self.pop_group, CVG, out_vcf_handle=VCF)

        if self.options.smartrerun:
            utils.safe_remove(utils.get_last_modification_file(self.cache_dir))

        cdef unsigned int batch_count = self.options.batch_count
        cdef unsigned int part_num = int(self.align_files.size/batch_count)
        if part_num * batch_count < self.align_files.size:
            part_num += 1

        cdef StringArray sub_align_files, batch_sample_ids
        cdef StringArray total_batch_files, batch_files
        strarray_init(&total_batch_files, part_num)

        cdef bint is_empty = True
        cdef size_t i = 0, j = 0, m = 0
        cdef char part_file_name[1024]

        cdef clock_t start_time, tot_start_time
        for k in range(self.regions.size):
            tot_start_time = clock()

            # set cache for fa sequence, this could make the program much faster
            # And remember that ``fa_file_hd`` is 0-base system, (start and end will be converted to
            # be 0-base in get_sequence function)
            self.fa_file_hd.set_cache_sequence(
                self.regions.array[k].chrom,
                c_max(0, self.regions.array[k].start - 5 * self.options.r_len),
                c_min(self.regions.array[k].end+self.options.r_len + 5 * self.options.r_len,
                      self.fa_file_hd.get_reference_length(self.regions.array[k].chrom)-1)
            )

            m = 0
            i = 0
            strarray_init(&batch_files, part_num)
            for i in range(0, self.align_files.size, batch_count):
                start_time = clock()

                m += 1
                # Join Path could fit different OS
                sprintf(part_file_name,"%s/basevar.%s.%ld.%ld.%d_%d.batch.gz", self.cache_dir,
                        self.regions.array[k].chrom, self.regions.array[k].start,
                        self.regions.array[k].end, m, part_num)

                strarray_append(&batch_files, part_file_name)

                # collect together will be convenient when we want to clear up these temporary files.
                strarray_append(&total_batch_files, part_file_name)

                if self.options.smartrerun and os.path.isfile(str(part_file_name)):
                    fprintf(stdout, "[INFO] %s already exists, we don't have to create it again, "
                           "when you set `smartrerun`\n", part_file_name)
                    continue
                else:
                    fprintf(stdout, "[INFO] Creating batch file %s\n", part_file_name)

                # One batch of alignment files, the size and order are the same with ``sub_align_files`` and
                # ``batch_sample_ids``
                strarray_init(&sub_align_files, batch_count)

                # could this loop be speeding up by prange in cython parallel?
                # https://cython.readthedocs.io/en/latest/src/userguide/parallelism.html
                # http://nealhughes.net/parallelcomp2/
                for j in range(i, i+batch_count):
                    strarray_append(&sub_align_files, self.align_files.array[i])

                if self.samples.size:
                    strarray_init(&batch_sample_ids, batch_count)
                    for j in range(i, i+batch_count):
                        strarray_append(&batch_sample_ids, self.samples.array[i])

                generate_batchfile(self.regions.array[k],   # 1-base in GenomeRegion
                                   &sub_align_files,
                                   &batch_sample_ids,
                                   self.fa_file_hd,
                                   part_file_name,
                                   self.options)

                strarray_destroy(&sub_align_files)
                strarray_destroy(&batch_sample_ids)
                fprintf(stdout, "[INFO] Done for batchfile %s , %.1f seconds elapsed.\n",
                        part_file_name, <double>(clock() - start_time)/CLOCKS_PER_SEC)

            # create batchfiles done
            fprintf(stdout, "[INFO] Created all batch files in %s:%ld-%ld for %ld samples are done,"
                            "%.1f seconds elapsed.\n", self.regions.array[k].chrom, self.regions.array[k].start,
                    self.regions.array[k].end, self.samples.size, <double>(clock() - tot_start_time)/CLOCKS_PER_SEC)

            start_time = clock()
            fprintf(stdout, "\n************************ Variants calling process *************************\n\n")
            try:
                _is_empty = variants_discovery(self.regions.array[k].chrom, &batch_files, self.pop_group,
                                               self.options.min_af, CVG, VCF)
            except Exception, e:
                fprintf(stderr, "[ERROR] Variants discovery in region %s:%ld-%ld\n",
                        self.regions.array[k].chrom, self.regions.array[k].start, self.regions.array[k].end+1)
                exit(EXIT_FAILURE)

            if not _is_empty:
                is_empty = False

            strarray_destroy(&batch_files)  # clear
            fprintf(stdout, "[INFO] Running variants_discovery in %s:%s-%s done, %.1f seconds elapsed.\n",
                    self.regions.array[k].chrom, self.regions.array[k].start, self.regions.array[k].end,
                    <double>(clock()-start_time)/CLOCKS_PER_SEC)

        CVG.close()
        if VCF:
            VCF.close()

        self.fa_file_hd.close()

        if is_empty:
            print("\n***************************************************************************\n"
                  "[WARNING] No reads are satisfy with the mapping quality (>=%d) in all of your\n"
                  "input files. We get nothing in %s \n\n" % (int(self.options.mapq), self.out_cvg_file))
            if VCF:
                print("and %s " % self.out_vcf_file)

        if REMOVE_BATCH_FILE:

            for i in total_batch_files.size:
                remove(total_batch_files.array[i])

            try:
                os.removedirs(self.cache_dir)
            except OSError:
                print("[WARNING] Directory not empty: %s, please delete it by yourself\n" % self.cache_dir)

        strarray_destroy(&total_batch_files)

        # double check
        with open(self.out_cvg_file + ".PROCESS.AND_VCF_DONE_SUCCESSFULLY", "w") as OUT:
            OUT.write("The process done.\n")

        return
