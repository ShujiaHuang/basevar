# cython: profile=True
"""This is a Process module for BaseType
"""
from libc.stdio cimport fprintf, stdout, stderr
from libc.stdlib cimport exit, EXIT_FAILURE, atoi, calloc, free
from libc.string cimport strsep
# from libc.time cimport clock_t, clock, CLOCKS_PER_SEC

from basevar.utils import vcf_header_define, cvg_header_define
from basevar.datatype.strarray cimport StringArray
from basevar.io.openfile import Open
# from basevar.io.fasta cimport FastaFile
# from basevar.io.bam cimport load_data_from_bamfile
# from basevar.io.htslibWrapper cimport Samfile

from basevar.caller.algorithm cimport strand_bias
from basevar.caller.algorithm cimport ref_vs_alt_ranksumtest

from basevar.caller.basetype cimport BaseType
from basevar.caller.batch cimport BatchGenerator, BatchInfo, PositionBatchCigarArray

cdef int INITIAL_CIGAR_ARRAY_SIZE = 10000
cdef int QUAL_THRESHOLD = 60
cdef list BASE = ['A', 'C', 'G', 'T']

def output_header(fa_file_name, sample_ids, pop_group_sample_dict, out_cvg_handle, out_vcf_handle=None):
    info, group = [], []
    if pop_group_sample_dict:
        for g in pop_group_sample_dict.keys():
            g_id = g.split('_AF')[0]  # ignore '_AF'
            group.append(g_id)
            info.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="Allele frequency in the %s '
                        'populations calculated base on LRT, in the range (0,1)">' % (g_id, g_id))

    if out_vcf_handle:
        vcf_header = vcf_header_define(fa_file_name, info="\n".join(info), samples=sample_ids)
        out_vcf_handle.write("%s\n" % "\n".join(vcf_header))

    out_cvg_handle.write('%s\n' % "\n".join(cvg_header_define(group)))

    return

cdef bint variants_discovery(const char *chrid, const StringArray *batchfiles, dict pop_group, float min_af,
                             cvg_file_handle, vcf_file_handle):
    """Function for variants discovery.
    """
    cdef list sampleinfos = []
    cdef list batch_files_hd = []

    cdef unsigned k = 0
    for k in range(batchfiles.size):
        batch_files_hd.append(Open(batchfiles.array[k], 'rb'))

    cdef bint is_empty = True
    cdef bint eof = False
    cdef bint is_error = False

    cdef int *batch_count = <int*> (calloc(batchfiles.size, sizeof(int)))

    cdef int total_sample_num = 0
    cdef BatchInfo batchinfo

    cdef int n = 0, i = 0
    cdef bint is_first = True
    while True:
        # [CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand]
        sampleinfos = []
        for i, fh in enumerate(batch_files_hd):
            line = fh.readline()
            if line:
                if line.startswith("#"):
                    if line.startswith("##SampleIDs="):
                        batch_count[i] = len(line.strip().split("=")[-1].split(","))
                        total_sample_num += batch_count[i]

                    continue

                elif is_first:
                    is_first = False
                    # initial here
                    batchinfo = BatchInfo(chrid, size=total_sample_num)

                sampleinfos.append(line.strip().split())
            else:
                sampleinfos.append(None)
                eof = True

        # hit the end of files
        if eof:
            is_error = True if any(sampleinfos) else False
            if is_error:
                for i in range(batchfiles.size):
                    fprintf(stdout, "%s\n", batchfiles.array[i])

                fprintf(stdout, "[ERROR]Error happen when 'variants_discovery', they don't "
                                "have the same positions in above files.\n")
            break

        # Empty!! Just get header information.
        if not sampleinfos:
            continue

        # reset position and ref_base
        batchinfo.position = int(sampleinfos[0][1])
        batchinfo.ref_base = sampleinfos[0][2]

        n += 1
        if n % 10000 == 0 or n == 1:
            fprintf(stdout, "[INFO] Have been loaded %d lines when hit position %s:%d\n",
                    n, chrid, batchinfo.position)

        # data in ``batchinfo`` will been updated automatically in this funcion
        _fetch_baseinfo_by_position_from_batchfiles(sampleinfos, batch_count, batchinfo)

        # ignore if coverage=0
        if batchinfo.depth == 0:
            continue

        # Not empty
        is_empty = False

        # Calling varaints position one by one and output files.
        _basetypeprocess(batchinfo, pop_group, min_af, cvg_file_handle, vcf_file_handle)

    for fh in batch_files_hd:
        fh.close()

    free(batch_count)
    fprintf(stdout, "[INFO] Have totally loaded %d lines in thes batch cluster. And the last hit poisition is %s:%d\n",
                    n, chrid, batchinfo.position)

    return is_empty


cdef void _fetch_baseinfo_by_position_from_batchfiles(list infolines, int *batch_count, BatchInfo batchinfo):

    # reset depth
    batchinfo.depth = 0

    cdef char *c_t4
    cdef char *c_t5
    cdef char *c_t6
    cdef char *c_t7
    cdef char *c_t8

    cdef int n = 0, index = 0, i = 0
    for i, col in enumerate(infolines):
        # <CHROM POS REF Depth MappingQuality Readbases ReadbasesQuality ReadPositionRank Strand>
        if len(col) == 0:
            fprintf(stdout, " %d lines happen to be empty in batchfiles!", i + 1)
            exit(EXIT_FAILURE)

        col[1], col[3] = map(int, [col[1], col[3]])
        if col[0] != batchinfo.chrid or col[1] != batchinfo.position or col[2] != batchinfo.ref_base:
            print ("%d lines, chromosome [%s and %s] or position [%d and %s] or ref-base [%s and %s] in "
                   "batchfiles not match with each other!\n" %
                   (i + 1, col[0], batchinfo.chrid, col[1], batchinfo.position, col[2], batchinfo.ref_base))
            exit(EXIT_FAILURE)

        batchinfo.depth += col[3]
        if col[3] > 0:

            c_t4, c_t5, c_t6, c_t7, c_t8 = col[4:9]
            for n in range(batch_count[i]):

                # if catch segmentation fault then the problem would probably be here!
                batchinfo.mapqs[index] = atoi(strsep(&c_t4, ","))
                batchinfo.sample_bases[index] = strsep(&c_t5, ",")  # must all be all upper charater in batchfile!
                batchinfo.sample_base_quals[index] = atoi(strsep(&c_t6, ","))
                batchinfo.read_pos_rank[index] = atoi(strsep(&c_t7, ","))
                batchinfo.strands[index] = strsep(&c_t8, ",")[0] # It's char not string

                # move to the next
                index += 1
        else:
            for n in range(batch_count[i]):

                batchinfo.mapqs[index] = 0
                batchinfo.sample_bases[index] = "N"
                batchinfo.sample_base_quals[index] = 0
                batchinfo.read_pos_rank[index] = 0
                batchinfo.strands[index] = "."

                # move to the next
                index += 1
    return


#####################################################################################################################
cdef void push_data_into_position_cigar_array(list regions_batch_cigar, list batch_generator_array, int sample_size):

    cdef PositionBatchCigarArray position_batch_cigar_array
    cdef BatchGenerator batch_generator
    cdef BatchInfo batch_info

    cdef int region_size = len(batch_generator_array)
    cdef int position_number
    cdef int i, j

    for i in range(region_size):
        batch_generator = batch_generator_array[i]
        position_number = len(batch_generator.batch_heap)

        for j in range(position_number):
            batch_info = batch_generator.batch_heap[j]
            if sample_size != -1:
                batch_info.set_size(sample_size)

            position_batch_cigar_array = regions_batch_cigar[i][j]
            position_batch_cigar_array.append(batch_info)
            batch_info.set_empty()

    return

# cdef bint variant_discovery_in_regions(FastaFile fa,
#                                        # StringArray *align_files,
#                                        # StringArray *samples,
#                                        # GenomeRegionArray *regions,
#                                        dict popgroup,
#                                        basestring out_cvg_file_name,
#                                        basestring out_vcf_file_name,
#                                        object options):
#
#     """
#     [ERROR] This function is imcomplete! ====================================================================

#     ``regions`` is a 2-D array, 1-base system
#         [[chr1, start1, end1], [chr1, start2, end2], ...]
#
#     ``samples``: The sample id of align_files
#     ``fa``:
#         # get sequence of chrom_name from reference fasta
#         fa = self.ref_file_hd.fetch(chrid)
#     """
#     cdef bytes chrom
#     cdef unsigned long start, end
#
#     cdef unsigned long region_size = len(regions)
#     cdef unsigned long sample_size = len(samples)
#
#     ### initial ###
#     cdef list regions_batch_cigar = []
#     cdef list positions_batch_cigar = []
#     cdef list batch_generators = []
#
#     cdef unsigned long _pos
#     for chrom, start, end in regions:
#
#         batch_generators.append(BatchGenerator(chrom, start, end, fa, options.batch_count, options))
#         positions_batch_cigar = []
#
#         for _pos in range(start, end+1):
#             # Position in positions_batch_cigar must be the same as which in `BatchGenerator.batch_heap`
#             positions_batch_cigar.append(PositionBatchCigarArray(
#                 chrom, _pos, fa.get_character(chrom, _pos-1), min(sample_size, INITIAL_CIGAR_ARRAY_SIZE))
#             )
#
#         # The size of ``regions_batch_cigar`` will be the same as ``batch_generators``
#         regions_batch_cigar.append(positions_batch_cigar)
#     ### initial done ###
#
#     fprintf(stdout, "[INFO] Done for allocating memory to ``PositionBatchCigarArray`` and ``BatchGenerator`` array.")
#
#     cdef Samfile reader
#
#     cdef int i = 0, k = 0, n = 0
#     cdef int buffer_sample_index = 0
#     cdef int size_of_batch_heap
#
#     cdef clock_t start_time
#     start_time = clock()
#     for i in range(sample_size):
#
#         reader = Samfile(align_files[i])  # Match samples[i]
#         reader.open("r", True)
#         for k ,(chrom, start, end) in enumerate(regions):
#
#             try:
#                 # load the whole mapping reads to ``be_generator`` in [chrom_name, start, end]
#                 load_data_from_bamfile(reader, samples[i], chrom, start, end, batch_generators[k],
#                                        buffer_sample_index, options)
#
#             except Exception, e:
#                 print ("Exception in region %s:%s-%s. Error: %s" % (chrom, start, end, e))
#                 exit(EXIT_FAILURE)
#
#         reader.close()
#
#         if buffer_sample_index + 1 == options.batch_count:
#             # Compress a batch data into ``PositionBatchCigarArray`` will rest depth to be 0
#             push_data_into_position_cigar_array(regions_batch_cigar, batch_generators, -1)
#             buffer_sample_index = 0
#         else:
#             buffer_sample_index += 1
#
#         # output some logger
#         n += 1
#         if n % 1000 == 0:
#             fprintf(stdout, "[INFO]Finish Loading %d bamfiles in all the regions, %.1f seconds "
#                         "elapsed in total.", n, <double>(clock() - start_time)/CLOCKS_PER_SEC)
#
#     fprintf(stdout, "[INFO] Finish Loading all %d bamfiles in all the regions, %d seconds "
#                 "elapsed in total.", n, <double>(clock() - start_time)/CLOCKS_PER_SEC)
#
#     if buffer_sample_index > 0:
#         push_data_into_position_cigar_array(regions_batch_cigar, batch_generators, buffer_sample_index)
#
#     if out_vcf_file_name:
#         VCF = Open(out_vcf_file_name, "wb", isbgz=True) if out_vcf_file_name.endswith(".gz") else \
#             open(out_vcf_file_name, "w")
#     else:
#         VCF = None
#
#     CVG = Open(out_cvg_file_name, "wb", isbgz=True) if out_cvg_file_name.endswith(".gz") else \
#         open(out_cvg_file_name, "w")
#
#     output_header(fa.filename, samples, popgroup, CVG, out_vcf_handle=VCF)
#     cdef bint is_empty = _variants_discovery(regions_batch_cigar, popgroup, options.min_af, CVG, VCF)
#
#     CVG.close()
#     if VCF:
#         VCF.close()
#
#     return is_empty


cdef bint _variants_discovery(list regions_batch_cigar, dict popgroup, float min_af, CVG, VCF):
    """Function for variants discovery.
    
    Parameter:
        ``start``: 1-base system
        ``end``: 1-base system
    """
    cdef int how_many_regions = len(regions_batch_cigar)
    cdef int how_many_pos

    cdef PositionBatchCigarArray position_batch_cigar_array
    cdef BatchInfo batch_info
    cdef bint is_empty = True
    cdef int n = 0, i = 0, j = 0
    for i in range(how_many_regions):

        how_many_pos = len(regions_batch_cigar[i])
        for j in range(how_many_pos):
            position_batch_cigar_array = regions_batch_cigar[i][j]
            batch_info = position_batch_cigar_array.convert_position_batch_cigar_array_to_batchinfo()
            if n % 10000 == 0:
                print ("Have been loading %d lines when hit position %s:%s" %
                       (n if n > 0 else 1, batch_info.chrid, batch_info.position))
            n += 1

            if batch_info.depth == 0:
                continue

            # Not empty
            is_empty = False

            # Calling varaints position one by one and output files.
            _basetypeprocess(batch_info, popgroup, min_af, CVG, VCF)

    return is_empty

cdef void _basetypeprocess(BatchInfo batchinfo, dict pop_group, float min_af, cvg_file_handle, vcf_file_handle):
    """
    
    :param batchinfo: 
    :param popgroup: 
    :param min_af: 
    :param cvg_file_handle: 
    :param vcf_file_handle: 
    :return: 
    """
    _out_cvg_file(batchinfo, pop_group, cvg_file_handle)

    cdef dict pop_group_bt = {}
    cdef bint is_variant = True

    cdef BaseType bt, group_bt
    cdef char ** group_sample_bases
    cdef int *group_sample_base_quals
    cdef int group_sample_size
    cdef int i = 0
    if vcf_file_handle:

        bt = BaseType()
        bt.cinit(batchinfo.ref_base.upper(), batchinfo.sample_bases, batchinfo.sample_base_quals,
                 batchinfo.size, min_af)

        is_variant = bt.lrt(None)  # do not need to set specific_base_combination
        if is_variant:

            pop_group_bt = {}
            for group, index in pop_group.items():

                group_sample_size = len(index)
                group_sample_bases = <char**> (calloc(group_sample_size, sizeof(char*)))
                if group_sample_bases == NULL:
                    fprintf(stderr,"[ERROR]Fail allocate memory for ``group_sample_bases`` in _basetypeprocess.")
                    exit(EXIT_FAILURE)

                group_sample_base_quals = <int*> (calloc(group_sample_size, sizeof(int)))
                if group_sample_base_quals == NULL:
                    fprintf(stderr, "[ERROR]Fail allocate memory for ``group_sample_base_quals`` in _basetypeprocess.")
                    exit(EXIT_FAILURE)

                # for i in index:
                for i in range(group_sample_size):
                    group_sample_bases[i] = batchinfo.sample_bases[index[i]]
                    group_sample_base_quals[i] = batchinfo.sample_base_quals[index[i]]

                group_bt = BaseType()
                group_bt.cinit(batchinfo.ref_base.upper(), group_sample_bases, group_sample_base_quals,
                               group_sample_size, min_af)

                group_bt.lrt([batchinfo.ref_base.upper()] + bt.alt_bases)
                pop_group_bt[group] = group_bt

                free(group_sample_bases)
                free(group_sample_base_quals)

            _out_vcf_line(batchinfo, bt, pop_group_bt, vcf_file_handle)
    return

cdef list _base_depth_and_indel(char ** bases, int size):
    # coverage info for each position
    cdef dict base_depth = {b: 0 for b in BASE}
    cdef dict indel_depth = {}

    cdef int i = 0
    for i in range(size):

        # The size of `base[i]` would be 1 except indel!
        if bases[i][0] == 'N':
            continue

        if bases[i] in base_depth:
            # ignore all bases('*') which not match ``cmm.BASE``
            base_depth[bases[i]] += 1
        else:
            # Indel
            indel_depth[bases[i]] = indel_depth.get(bases[i], 0) + 1

    cdef bytes indels = bytes(','.join(
        [k + '|' + str(v) for k, v in indel_depth.items()]
    ) if indel_depth else ".")

    return [base_depth, indels]

cdef void _out_cvg_file(BatchInfo batchinfo, dict popgroup, out_file_handle):
    """output coverage information into `out_file_handle`"""
    # coverage info for each position
    cdef dict base_depth
    cdef bytes indels
    base_depth, indels = _base_depth_and_indel(batchinfo.sample_bases, batchinfo.size)

    # base depth and indels for each subgroup
    cdef char ** group_sample_bases
    cdef int group_sample_size
    cdef dict group_cvg = {}
    cdef bytes group
    cdef list index

    # Here is one of the two parts which most time consuming!
    cdef int i = 0
    cdef dict sub_bd
    cdef bytes sub_inds
    for group, index in popgroup.items():

        group_sample_size = len(index)
        group_sample_bases = <char**> (calloc(group_sample_size, sizeof(char*)))
        if group_sample_bases == NULL:
            fprintf(stderr, "[ERROR]Fail allocate memory for ``group_sample_bases`` in _out_cvg_file.")
            exit(EXIT_FAILURE)

        for i in range(group_sample_size):
            group_sample_bases[i] = batchinfo.sample_bases[index[i]]

        sub_bd, sub_inds = _base_depth_and_indel(group_sample_bases, group_sample_size)
        group_cvg[group] = [sub_bd, sub_inds]

        free(group_sample_bases)

    cdef double fs, sor
    cdef int ref_fwd, ref_rev, alt_fwd, alt_rev
    fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = 0, -1, 0, 0, 0, 0

    # Here is the other part which most time consuming!
    cdef bytes ref_base = batchinfo.ref_base
    cdef bytes b1, b2
    if sum(base_depth.values()) > 0:
        # could we stop sorting ?
        base_sorted = sorted(base_depth.items(), key=lambda x: x[1], reverse=True)
        b1, b2 = base_sorted[0][0], base_sorted[1][0]

        # This is very efficient.
        fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
            ref_base.upper(),  # reference
            [b1 if b1 != ref_base.upper() else b2],  # alt-allele
            batchinfo.sample_bases,
            batchinfo.strands,
            batchinfo.size
        )

    cdef list group_info
    if sum(base_depth.values()):

        group_info = []
        if group_cvg:
            for k in group_cvg.keys():
                depth, indel = group_cvg[k]
                indel = [indel] if indel != "." else []
                s = ':'.join(map(str, [depth[b] for b in BASE]) + indel)
                group_info.append(s)

        out_file_handle.write(
            '\t'.join(
                [batchinfo.chrid, str(batchinfo.position), ref_base, str(sum(base_depth.values()))] +
                [str(base_depth[b]) for b in BASE] +
                [indels] +
                [str("%.3f" % fs), str("%.3f" % sor), ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev]))] +
                group_info
            ) + '\n'
        )

    return

cdef void _out_vcf_line(BatchInfo batchinfo, BaseType bt, dict pop_group_bt, out_file_handle):
    """output vcf lines into `out_file_handle`"""

    cdef dict alt_gt = {b: './' + str(k + 1) for k, b in enumerate(bt.alt_bases)}
    cdef list samples = []
    cdef int k
    cdef char *b
    # for k, b in enumerate(bases):
    for k in range(batchinfo.size):

        b = batchinfo.sample_bases[k]
        # For sample FORMAT
        if b[0] not in ['N', '-', '+']:
            # For the base which not in bt.alt_bases()
            if b not in alt_gt:
                alt_gt[b] = './.'

            gt = '0/.' if b == batchinfo.ref_base.upper() else alt_gt[b]

            samples.append(gt + ':' + b + ':' + chr(batchinfo.strands[k]) + ':' +
                           str(round(bt.qual_pvalue[k], 6)))
        else:
            samples.append('./.')  # 'N' base or indel

    # Rank Sum Test for mapping qualities of REF versus ALT reads
    mq_rank_sum = ref_vs_alt_ranksumtest(batchinfo.ref_base.upper(), bt.alt_bases, batchinfo.sample_bases,
                                         batchinfo.mapqs, batchinfo.size)

    # Rank Sum Test for variant appear position among read of REF versus ALT
    read_pos_rank_sum = ref_vs_alt_ranksumtest(batchinfo.ref_base.upper(), bt.alt_bases, batchinfo.sample_bases,
                                               batchinfo.read_pos_rank, batchinfo.size)

    # Rank Sum Test for base quality of REF versus ALT
    base_q_rank_sum = ref_vs_alt_ranksumtest(batchinfo.ref_base.upper(), bt.alt_bases, batchinfo.sample_bases,
                                             batchinfo.sample_base_quals, batchinfo.size)

    # Variant call confidence normalized by depth of sample reads
    # supporting a variant.
    ad_sum = sum([bt.depth[b] for b in bt.alt_bases])
    qd = round(float(bt.var_qual / ad_sum), 3)

    # Strand bias by fisher exact test and Strand bias estimated by the
    # Symmetric Odds Ratio test
    fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
        batchinfo.ref_base.upper(),
        bt.alt_bases,
        batchinfo.sample_bases,
        batchinfo.strands,
        batchinfo.size
    )

    # base=>[CAF, allele depth], CAF = Allele frequency by read count
    caf = {b: ['%f' % round(bt.depth[b] / float(bt.total_depth), 6),
               bt.depth[b]] for b in bt.alt_bases}

    info = {'CM_DP': str(int(bt.total_depth)),
            'CM_AC': ','.join(map(str, [caf[b][1] for b in bt.alt_bases])),
            'CM_AF': ','.join(map(str, [bt.af_by_lrt[b] for b in bt.alt_bases])),
            'CM_CAF': ','.join(map(str, [caf[b][0] for b in bt.alt_bases])),
            'MQRankSum': str("%.3f" % mq_rank_sum) if mq_rank_sum != -1 else 'nan',
            'ReadPosRankSum': str("%.3f" % read_pos_rank_sum) if read_pos_rank_sum != -1 else 'nan',
            'BaseQRankSum': str("%.3f" % base_q_rank_sum) if base_q_rank_sum != -1 else 'nan',
            'QD': str(qd),
            'SOR': str("%.3f" % sor),
            'FS': str("%.3f" % fs),
            'SB_REF': str(ref_fwd) + ',' + str(ref_rev),
            'SB_ALT': str(alt_fwd) + ',' + str(alt_rev)}

    cdef BaseType g_bt
    cdef bytes group

    if pop_group_bt:
        for group, g_bt in pop_group_bt.items():
            af = ','.join(map(str, [g_bt.af_by_lrt[bb] if bb in g_bt.af_by_lrt else 0
                                    for bb in bt.alt_bases]))
            info[group] = af

    out_file_handle.write('\t'.join([batchinfo.chrid, str(batchinfo.position), '.', batchinfo.ref_base,
                                     ','.join(bt.alt_bases), str(bt.var_qual),
                                     '.' if bt.var_qual > QUAL_THRESHOLD else 'LowQual',
                                     ';'.join([kk + '=' + vv for kk, vv in sorted(
                                         info.items(), key=lambda x: x[0])]),
                                     'GT:AB:SO:BP'] + samples) + '\n')
    return
