/**
 * @file basetype_utils.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */
#include <sstream>
#include <fstream>
#include <ctime>      // clock
#include <algorithm>  // std::min

#include <htslib/bgzf.h>
#include <htslib/tbx.h>

#include "basetype_caller.h"
#include "external/threadpool.h"

void BaseTypeRunner::set_arguments(int cmd_argc, char *cmd_argv[]) {

    if (cmd_argc < 2) {
        std::cout << usage() << "\n" << std::endl;
        exit(1);
    }

    if (_args) {
        throw std::runtime_error("[basetype.cpp::BaseTypeRunner:args] 'args' must be "
                                 "a NULL pointer before it can be assigned a value.");
    }
    // Inital a new BasTypeARGS and set defaut argument.
    _args = new BaseTypeARGS;
    
    // Parsing the commandline options. 
    char c;
    while((c = getopt_long(cmd_argc, cmd_argv, "I:L:R:m:q:B:t:r:G:h", BASETYPE_CMDLINE_LOPTS, NULL)) >= 0) {
        // 字符流解决命令行参数转浮点等类型的问题
        std::stringstream ss(optarg ? optarg: "");  
        switch (c) {
            case 'I': _args->input_bf.push_back(optarg);         break;  // 恒参 (一直用)
            case 'L': _args->in_bamfilelist = optarg;            break;  /* 临参 */
            case 'R': _args->reference      = optarg;            break;  /* 临参 */

            case 'm': ss >> _args->min_af;                       break;  // 恒参
            case 'q': ss >> _args->mapq;                         break;  // 恒参
            case 'B': ss >> _args->batchcount;                   break;  // 恒参
            case 't': ss >> _args->thread_num;                   break;  // 恒参

            case 'r': _args->regions        = optarg;            break;  /* 临参 */
            case 'G': _args->pop_group_file = optarg;            break;  // 
            case '1': _args->output_vcf     = optarg;            break;  // 恒参
            case '2': _args->output_cvg     = optarg;            break;  // 恒参

            case '3': _args->filename_has_samplename = true;     break;  // 恒参
            case '4': _args->smart_rerun             = true;     break;  // 恒参
            case 'h': std::cout << usage() << std::endl; exit(1);
            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                exit(1);
        }
    }
    /* Make sure we set valid arguments */
    if (_args->input_bf.empty() && _args->in_bamfilelist.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-I/--input' or '-L/--align-file-list'");
    if (_args->reference.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-R/--reference'");

    if (_args->output_vcf.empty())
        throw std::invalid_argument("[ERROR] Missing argument '--output-vcf'");
    if (_args->output_cvg.empty())
        throw std::invalid_argument("[ERROR] Missing argument '--output-cvg'");
    
    if (_args->min_af <= 0)
        throw std::invalid_argument("[ERROR] '-m/--min-af' argument must be > 0");
    if (_args->mapq <= 0)
        throw std::invalid_argument("[ERROR] '-q/--mapq' argument must be > 0");
    if (_args->batchcount <= 0)
        throw std::invalid_argument("[ERROR] '-B/--batch-count' argument must be > 0");
    if (_args->thread_num <= 0)
        throw std::invalid_argument("[ERROR] '-t/--thread' argument must be > 0");

    // recovering the absolute paths of output files
    _args->output_vcf = ngslib::abspath(_args->output_vcf);
    _args->output_cvg = ngslib::abspath(_args->output_cvg);

    // Output the commandline options
    std::cout << 
        "[INFO] BaseVar arguments:\n"
        "basevar basetype -R " + _args->reference        + " \\ \n"  + (_args->input_bf.empty() ? "" : 
        "   -I " + ngslib::join(_args->input_bf, " -I ") + " \\ \n") + (_args->in_bamfilelist.empty() ? "" : 
        "   -L " + _args->in_bamfilelist       + " \\ \n") <<
        "   -q " << _args->mapq               << " \\ \n"
        "   -m " << _args->min_af             << " \\ \n"
        "   -B " << _args->batchcount         << " \\ \n"
        "   -t " << _args->thread_num         << " \\ \n"  << (_args->regions.empty() ? "" : 
        "   -r " + _args->regions              + " \\ \n") << (_args->pop_group_file.empty() ? "" : 
        "   -G " + _args->pop_group_file       + " \\ \n") <<
        "   --output-vcf " + _args->output_vcf + " \\ \n"
        "   --output-vcg " + _args->output_cvg << (_args->filename_has_samplename ? " \\ \n"
        "   --filename-has-samplename" : "")   << (_args->smart_rerun ? " \\ \n"
        "   --smart-rerun": "")                << "\n" << std::endl;
    
    if (_args->smart_rerun) {
        std::cout << "************************************************\n"
                     "******************* WARNING ********************\n"
                     "************************************************\n"
                     ">>>>>>>> You have setted `smart rerun` <<<<<<<<<\n"
                     "Please make sure all the parameters are the same\n"
                     "with your previous commands.\n"
                     "************************************************\n"
                     "************************************************\n\n";
    }

    if (!_args->in_bamfilelist.empty()) {
        std::vector<std::string> filelist = get_firstcolumn_from_file(_args->in_bamfilelist);
        _args->input_bf.insert(_args->input_bf.end(), filelist.begin(), filelist.end());
    }
    std::cout << "[INFO] Finish loading arguments and we have " << _args->input_bf.size()
              << " BAM/CRAM files for variants calling.\n"      << std::endl;

    // Setting the resolution of AF
    _args->min_af = std::min(float(100)/_args->input_bf.size(), _args->min_af);

    // load fasta
    reference = _args->reference;

    _get_calling_interval();
    print_calling_interval();
    _get_sample_id_from_bam();  // keep the order of '_samples_id' as the same as input 'aligne_files'

    // check the input bamfiles have duplicate or not
    std::vector<std::string> duplicate_samples = ngslib::find_duplicates(_samples_id);
    if (!duplicate_samples.empty()) {
        std::cout << "[WARNING] Find " << duplicate_samples.size() << " duplicated samples within " 
                  << "the input bamfiles: " + ngslib::join(duplicate_samples, ",") + "\n" 
                  << std::endl;
    }

    if (!_args->pop_group_file.empty()) 
        _get_popgroup_info();

    return;
}

// Run the processes of calling variant and output files.
void BaseTypeRunner::run() {
    _variant_caller_process();
    return;
}

void BaseTypeRunner::_variant_caller_process() {

    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);

    // Get filepath and stem name first.
    std::string _bname = ngslib::basename(_args->output_vcf);
    size_t si = _bname.find(".vcf");
    std::string stem_bn = (si > 0 && si != std::string::npos) ? _bname.substr(0, si) : _bname;

    std::string outdir = ngslib::dirname(_args->output_vcf);
    std::string cache_outdir = outdir + "/cache_" + stem_bn;

    if (IS_DELETE_CACHE_BATCHFILE && ngslib::path_exists_and_not_empty(cache_outdir)) {
        throw std::runtime_error("[ERROR] [" + cache_outdir + "] must be an empty folder. "
                                 "You can delete it manual before execute BaseVar.");
    }
    ngslib::safe_mkdir(cache_outdir);  // make cache directory for batchfiles

    if (_args->smart_rerun) {
        // Remove and rollback `thread_num` last modification files. 
        // Must do is before calling '_create_batchfiles'
        for (size_t i(0); i < _args->thread_num; ++i)
            ngslib::safe_remove(ngslib::get_last_modification_file(cache_outdir));
    }

    // 以区间为单位进行变异检测, 每个区间里调用多线程
    std::vector<std::string> batchfiles, vcffiles, cvgfiles;
    std::string ref_id; uint32_t reg_start, reg_end;
    for (size_t i(0); i < _calling_intervals.size(); ++i) {

        // Region information
        std::tie(ref_id, reg_start, reg_end) = _calling_intervals[i];
        std::string regstr = ref_id + "_" + ngslib::tostring(reg_start) + "-" + ngslib::tostring(reg_end);

        ///////////////////////////////////////////////////////////
        ///////// Create batchfiles with multiple-thread //////////
        ///////////////////////////////////////////////////////////
        std::string prefix = cache_outdir + "/" + stem_bn + "." + regstr;     
        batchfiles = _create_batchfiles(_calling_intervals[i], prefix);

        // Time information
        time_t now = time(0);
        std::string ct(ctime(&now)); 
        ct.pop_back();
        std::cout << "[INFO] "+ ct +". Done for creating all " << batchfiles.size() << " batchfiles in " 
                  << regstr + " and start to call variants, "  << difftime(now, real_start_time)
                  << " (CPU time: " << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
                  << ") seconds elapsed in total.\n" << std::endl;

        ///////////////////////////////////////////////////////////
        // Calling variants from batchfiles with mulitple thread //
        ///////////////////////////////////////////////////////////
        std::string sub_vcf_fn = prefix + ".vcf.gz";
        std::string sub_cvg_fn = prefix + ".cvg.gz";
        _variants_discovery(batchfiles, _calling_intervals[i], sub_vcf_fn, sub_cvg_fn);

        // Time information
        now = time(0);
        ct  = ctime(&now); 
        ct.pop_back();
        std::cout << "\n" << "[INFO] "+ ct +". Done for calling all variants in " + regstr + ": "
                  << "[" + sub_vcf_fn + "," + sub_cvg_fn + "], " << difftime(now, real_start_time)
                  << " (CPU time: " << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
                  << ") seconds elapsed in total.\n" << std::endl;

        vcffiles.push_back(sub_vcf_fn);
        cvgfiles.push_back(sub_cvg_fn);

        if (IS_DELETE_CACHE_BATCHFILE) {
            for (auto bf: batchfiles) {
                ngslib::safe_remove(bf);
                ngslib::safe_remove(bf+".tbi");
            }
        }
    }

    // Merge all VCF subfiles
    std::vector<std::string> add_group_info, group_name;
    std::map<std::string, std::vector<size_t>>::iterator it = _groups_idx.begin();
    for(; it != _groups_idx.end(); ++it) {
        group_name.push_back(it->first);
        add_group_info.push_back("##INFO=<ID=" + it->first + "_AF,Number=A,Type=Float,Description="
                                 "\"Allele frequency in the " + it->first + " populations calculated "
                                 "base on LRT, in the range (0,1)\">");
    }

    // Merge VCF
    std::string header = vcf_header_define(_args->reference, add_group_info, _samples_id);
    merge_file_by_line(vcffiles, _args->output_vcf, header, true);

    const tbx_conf_t bf_tbx_conf = {1, 1, 2, 0, '#', 0};  // {preset, seq col, beg col, end col, header-char, skip-line}
    if ((ngslib::suffix_name(_args->output_vcf) == ".gz") &&          // create index 
        tbx_index_build(_args->output_vcf.c_str(), 0, &bf_tbx_conf))  // file suffix is ".tbi"
        throw std::runtime_error("tbx_index_build failed: Is the file bgzip-compressed? "
                                 "Check this file: " + _args->output_vcf + "\n");

    // Merge CVG
    header = cvg_header_define(group_name, BASES);
    merge_file_by_line(cvgfiles, _args->output_cvg, header, true);
    if ((ngslib::suffix_name(_args->output_cvg) == ".gz") &&          // create index
        tbx_index_build(_args->output_cvg.c_str(), 0, &bf_tbx_conf))  // file suffix is ".tbi"
        throw std::runtime_error("tbx_index_build failed: Is the file bgzip-compressed? "
                                 "Check this file: " + _args->output_cvg + "\n");

    if (IS_DELETE_CACHE_BATCHFILE) {
        ngslib::safe_remove(cache_outdir);
    }
    return;
}

void BaseTypeRunner::_get_sample_id_from_bam() {
    time_t real_start_time = time(0);

    // Loading sample ID in BAM/CRMA files from RG tag.
    if (_args->filename_has_samplename)
        std::cout << "[INFO] BaseVar'll load samples id from filename directly, becuase you set "
                     "--filename-has-samplename.\n";

    std::string samplename, filename;
    size_t si;
    for (size_t i(0); i < _args->input_bf.size(); ++i) {

        if ((i+1) % 1000 == 0)
            std::cout << "[INFO] loading "   << i+1 << "/" << _args->input_bf.size() 
                      << " alignment files." << std::endl;
        
        if (_args->filename_has_samplename) {
            filename = ngslib::remove_filename_extension(ngslib::basename(_args->input_bf[i]));
            si = filename.find('.');
            samplename = si > 0 && si != std::string::npos ? filename.substr(0, si) : filename;
        } else {
            // Get sampleID from BAM header, a bit time-consuming.
            ngslib::BamHeader bh(_args->input_bf[i]);
            samplename = bh.get_sample_name();
        }

        if (!samplename.empty()) {
            _samples_id.push_back(samplename);
        } else {
            throw std::invalid_argument("[BaseTypeRunner::_load_sample_id_from_bam] " + 
                                        _args->input_bf[i] + " sample ID not found.\n");
        }
    }

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "[INFO] " + ct + ". Done for loading all samples' id from alignment files, " 
              << difftime(real_start_time, now) << " seconds elapsed.\n" 
              << std::endl;

    return;
}

void BaseTypeRunner::_get_calling_interval() {

    // clear
    _calling_intervals.clear();
    if (!_args->regions.empty()) {
        std::vector<std::string> rg_v;
        ngslib::split(_args->regions, rg_v, ",");

        for (size_t i(0); i < rg_v.size(); ++i) {
            if (rg_v[i].length() == 0) continue; // ignore empty string 
            _calling_intervals.push_back(_make_gregion_tuple(rg_v[i]));
        }
    } else {
        // Calling the whole genome
        int n = reference.nseq();
        for (size_t i(0); i < n; ++i) {
            std::string ref_id = reference.iseq_name(i);
            _calling_intervals.push_back(std::make_tuple(ref_id, 1, reference.seq_length(ref_id)));
        } 
    }
    
    return;
}

ngslib::GenomeRegionTuple BaseTypeRunner::_make_gregion_tuple(std::string gregion) {

    // Genome Region, 1-based
    std::string ref_id; 
    uint32_t reg_start, reg_end;

    std::vector<std::string> gr;
    ngslib::split(gregion, gr, ":");     // get reference id
    ref_id = gr[0];

    if (gr.size() == 2) {                // 'start-end' or start
        std::vector<uint32_t> gs;
        ngslib::split(gr[1], gs, "-");   // get position coordinate
        reg_start = gs[0];
        reg_end = (gs.size() == 2) ? gs[1] : reference.seq_length(ref_id);
    } else {
        reg_start = 1;
        reg_end = reference.seq_length(ref_id);  // the whole ``ref_id`` length
    }

    if (reg_start > reg_end) {
        throw std::invalid_argument("[ERROR] start postion is larger than end position in "
                                    "-r/--regions " + gregion);
    }

    return std::make_tuple(ref_id, reg_start, reg_end);  // 1-based
}

void BaseTypeRunner::print_calling_interval() {

    std::string ref_id;
    uint32_t reg_start, reg_end;
    std::cout << "---- Calling Intervals ----\n";
    for (size_t i(0); i < _calling_intervals.size(); ++i) {
        std::tie(ref_id, reg_start, reg_end) = _calling_intervals[i];
        std::cout << i+1 << " - " << ref_id << ":" << reg_start << "-" << reg_end << "\n";
    }
    std::cout << "\n";
    return;
}

void BaseTypeRunner::_get_popgroup_info() {
    // group_id => [index in _samples_id of group_id]
    std::ifstream i_fn(_args->pop_group_file.c_str());
    if (!i_fn) {
        std::cerr << "[ERROR] Cannot open file: " + _args->pop_group_file << std::endl;
        exit(1);
    }

    std::map<std::string, std::string> sample2group;
    std::string skip, sn, gn;
    while (1) {
        // Only two columns: sample_id and group_id
        i_fn >> sn >> gn;
        if (i_fn.eof()) break;
        
        sample2group[sn] = gn;
        std::getline(i_fn, skip, '\n');  // skip the rest information of line.
    }
    i_fn.close();

    _groups_idx.clear();
    std::map<std::string, std::string>::iterator s2g_it;
    
    // follow the order of '_samples_id'
    if (sample2group.size() > 0) {
        for (size_t i(0); i < _samples_id.size(); ++i) {
            s2g_it = sample2group.find(_samples_id[i]);
            
            // ignore all the samples which not found
            if (s2g_it != sample2group.end()) {
                // record sample index of group groups
                // group -> index of _samples_id
                _groups_idx[sample2group[_samples_id[i]]].push_back(i);
            }
        }
    }

    return;
}

std::vector<std::string> BaseTypeRunner::_create_batchfiles(const ngslib::GenomeRegionTuple &genome_region, 
                                                            const std::string bf_prefix) 
{
    std::string ref_id; uint32_t reg_beg, reg_end;
    std::tie(ref_id, reg_beg, reg_end) = genome_region;
    std::string fa_seq = reference[ref_id];  // use the whole sequence of ``ref_id`` for simply

    int bn = _args->input_bf.size() / _args->batchcount;  // number of batchfiles
    if (_args->input_bf.size() % _args->batchcount > 0)
        bn++;

    ThreadPool thread_pool(_args->thread_num);  // set multiple-thread
    std::vector<std::future<bool>> create_batchfile_processes;

    std::vector<std::string> batchfiles;
    for (size_t i(0), j(1); i < _args->input_bf.size(); i+=_args->batchcount, ++j) {
        // set name of batchfile and must be compressed by BGZF.
        std::string batchfile = bf_prefix + "." + ngslib::tostring(j) + "_" + ngslib::tostring(bn) + ".bf.gz";
        batchfiles.push_back(batchfile);   // Store the name of batchfile into a vector

        if (_args->smart_rerun && ngslib::is_readable(batchfile)) {
            // do not have to create the exists batchfiles again if set `--smart-rerun`
            std::cout << "[INFO] " + batchfile + " already exists, we don't have to "
                         "create it again, when we set `--smart-rerun`.\n";
            continue;
        }

        // slicing bamfiles for a batchfile.
        size_t x(i), y(i + _args->batchcount);
        std::vector<std::string> batch_align_files = ngslib::vector_slicing(_args->input_bf, x, y);
        std::vector<std::string> batch_sample_ids  = ngslib::vector_slicing(_samples_id, x, y);

        // make Thread Pool
        create_batchfile_processes.emplace_back(
            thread_pool.enqueue(__create_a_batchfile, 
                                batch_align_files,  // 局部变量，会变，必拷贝，不可传引用，否则线程执行时将丢失该值
                                batch_sample_ids,   // 局部变量，会变，必拷贝，不可传引用，否则线程执行时将丢失该值
                                std::cref(fa_seq),  // 外部变量，不变，传引用，省内存
                                std::cref(genome_region),
                                _args->mapq,
                                batchfile));
    }
    
    for (auto && p: create_batchfile_processes) {
        // Run and make sure all processes are finished
        // 一般来说，只有当 valid() 返回 true 的时候才调用 get() 去获取结果，这也是 C++ 文档推荐的操作。
        if (p.valid()) {
            // get() 调用会改变其共享状态，不再可用，也就是说 get() 只能被调用一次，多次调用会触发异常。
            // 如果想要在多个线程中多次获取产出值需要使用 shared_future。
            bool x = p.get(); // retrieve the return value of `__create_a_batchfile`
        }
    }
    create_batchfile_processes.clear();  // release the thread

    return batchfiles;  // done for batchfiles created and return
}

void BaseTypeRunner::_variants_discovery(const std::vector<std::string> &batchfiles, 
                                         const ngslib::GenomeRegionTuple &genome_region,
                                         const std::string out_vcf_fn,
                                         const std::string out_cvg_fn) 
{
    static const uint32_t STEP_REGION_LEN = 100000; // 10 for test, should set to be large than 100000

    // get region information
    std::string ref_id; uint32_t reg_beg, reg_end, sub_reg_start, sub_reg_end;
    std::tie(ref_id, reg_beg, reg_end) = genome_region;

    int bn = (reg_end - reg_beg + 1) / STEP_REGION_LEN;
    if ((reg_end - reg_beg + 1) % STEP_REGION_LEN > 0)
        bn++;

    // prepare multiple-thread
    ThreadPool thread_pool(_args->thread_num);  
    std::vector<std::future<bool>> call_variants_processes;

    std::vector<std::string> subvcfs, subcvgs;
    for (uint32_t i(reg_beg), j(1); i < reg_end + 1; i += STEP_REGION_LEN, ++j) {

        std::string tmp_vcf_fn = out_vcf_fn + "." + ngslib::tostring(j) + "_" + ngslib::tostring(bn);
        std::string tmp_cvg_fn = out_cvg_fn + "." + ngslib::tostring(j) + "_" + ngslib::tostring(bn);
        subvcfs.push_back(tmp_vcf_fn);
        subcvgs.push_back(tmp_cvg_fn);

        sub_reg_start = i;
        sub_reg_end = (sub_reg_start+STEP_REGION_LEN-1) > reg_end ? reg_end : sub_reg_start+STEP_REGION_LEN-1;
        std::string regstr = ref_id + ":" + ngslib::tostring(sub_reg_start) + "-" + ngslib::tostring(sub_reg_end);

        // Performance multi-thread here.
        call_variants_processes.emplace_back(
            thread_pool.enqueue(_variant_calling_unit, 
                                std::cref(batchfiles), 
                                std::cref(_samples_id),
                                std::cref(_groups_idx),
                                _args->min_af,
                                regstr,        // 局部变量必须拷贝，会变 
                                tmp_vcf_fn,    // 局部变量必须拷贝，会变 
                                tmp_cvg_fn));  // 局部变量必须拷贝，会变 
    }

    // Run and make sure all processes could be finished.
    for (auto & p: call_variants_processes) {
        if (p.valid()) {
            bool x = p.get();
        }
    }
    call_variants_processes.clear();  // release the thread

    std::string header = "## No need header here";
    merge_file_by_line(subvcfs, out_vcf_fn, header, true);
    merge_file_by_line(subcvgs, out_cvg_fn, header, true);

    return;
}

/// Functions for calling variants outside of 'BaseTypeRunner' class 
// A unit for calling variants and let it run in a thread.
bool _variant_calling_unit(const std::vector<std::string> &batchfiles, 
                           const std::vector<std::string> &sample_ids,  // total samples
                           const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                           const double min_af,
                           const std::string region,    // genome region format like samtools
                           const std::string tmp_vcf_fn,
                           const std::string tmp_cvg_fn) 
{
    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);

    /*************** Preparing for reading data **************/
    std::vector<std::string> bf_smp_ids = _get_sampleid_from_batchfiles(batchfiles);
    // Check smaples id
    if (ngslib::join(bf_smp_ids, ",") != ngslib::join(sample_ids, ","))
        throw std::runtime_error("[BUG] The order of sample ids in batchfiles must be the same as "
                                 "input bamfiles.\n" 
                                 "Sample ids in batchfiles: " + ngslib::join(bf_smp_ids, ",") + "\n" 
                                 "Sample ids in bamfiles  : " + ngslib::join(sample_ids, ",") + "\n");

    std::vector<BGZF*> batch_file_hds;
    std::vector<tbx_t*> batch_file_tbx;
    std::vector<hts_itr_t*> batch_file_itr;
    for (size_t i(0); i < batchfiles.size(); ++i) {

        // Get and check the sample id from batchfiles header.
        BGZF *f = bgzf_open(batchfiles[i].c_str(), "r");
        f = bgzf_open(batchfiles[i].c_str(), "r"); // open again
        batch_file_hds.push_back(f);               // record the file handle

        tbx_t *tbx = tbx_index_load(batchfiles[i].c_str());
        if (!tbx) {
            throw std::runtime_error("[ERROR] " + batchfiles[i] + " tabix index load failed.");
        }
        batch_file_tbx.push_back(tbx);

        hts_itr_t *itr =NULL;
        if ((itr = tbx_itr_querys(tbx, region.c_str())) == 0) {
            throw std::runtime_error("[ERROR] load data in " + region + " failed. "
                                     "Check input file: " + batchfiles[i]);
        }
        batch_file_itr.push_back(itr);
    }
    /************* Done for reading data prepare *************/

    // Start calling variants
    BGZF *VCF = bgzf_open(tmp_vcf_fn.c_str(), "w"); // output vcf file
    if (!VCF) throw std::runtime_error("[ERROR] " + tmp_vcf_fn + " open failure.");
    
    BGZF *CVG = bgzf_open(tmp_cvg_fn.c_str(), "w"); // output coverage file
    if (!CVG) throw std::runtime_error("[ERROR] " + tmp_cvg_fn + " open failure.");

    std::vector<std::string> smp_bf_line_vector;
    smp_bf_line_vector.reserve(batchfiles.size());  // size number is as the same as batchfiles.

    bool is_eof(false), has_data(false);
    uint32_t n = 0;
    while (!is_eof) {
        // clear the vector before next loop.
        smp_bf_line_vector.clear();
        for (size_t i(0); i < batchfiles.size(); ++i) {
            
            // Fetch one line from each batchfile per loop and the row number of batchfiles must be the same.
            kstring_t s; s.s = NULL; s.l = s.m = 0; // must be refreshed in loop
            int eof = tbx_bgzf_itr_next(batch_file_hds[i], batch_file_tbx[i], batch_file_itr[i], &s);
            if (eof < 0) {
                is_eof = true;
                free(s.s);
                break;
            }
            smp_bf_line_vector.push_back(s.s);
            free(s.s);
        }
        
        ++n;  // line number
        if (n % 10000 == 0) {
            std::cout << "[INFO] Have been loaded " << n << " lines.\n";
        }

        // Samples' data for each poistion is ready and calling variants now
        bool cc = _basevar_caller(smp_bf_line_vector, group_smp_idx, min_af, sample_ids.size(), VCF, CVG);
        if (cc && !has_data) has_data = cc;
    }

    // close all batchfiles
    for (size_t i(0); i < batchfiles.size(); ++i) {
        tbx_destroy(batch_file_tbx[i]);
        tbx_itr_destroy(batch_file_itr[i]);
        if ((bgzf_close(batch_file_hds[i]) < 0))
            throw std::runtime_error("[ERROR] " + batchfiles[i] + " fail close.");
    }

    // close VCF and CVG file
    if (bgzf_close(VCF) < 0) throw std::runtime_error("[ERROR] " + tmp_vcf_fn + " fail close.");
    if (bgzf_close(CVG) < 0) throw std::runtime_error("[ERROR] " + tmp_cvg_fn + " fail close.");

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();
    std::cout << "[INFO] " + ct + ". Done for creating [" + tmp_vcf_fn + ", " + tmp_cvg_fn + "], "
              << difftime(now, real_start_time) 
              << " (CPU time: " << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
              << ") seconds elapsed." << std::endl;

    return has_data;
}

std::vector<std::string> _get_sampleid_from_batchfiles(const std::vector<std::string> &batchfiles) {
    // Get sample id from batchfiles header.

    std::vector<std::string> bf_smp_ids;
    for (size_t i(0); i < batchfiles.size(); ++i) {

        BGZF *f = bgzf_open(batchfiles[i].c_str(), "r");
        kstring_t s; s.s = NULL; s.l = s.m = 0;

        // get all the samples' id from the header of batchfiles.
        while (bgzf_getline(f, '\n', &s) >= 0) {
            if (s.s[0] != '#') { // head char is '#' in batchfile.
                break;
            } else if (strncmp(s.s, "##SampleIDs=", 12) == 0) {
                // Header looks like: ##SampleIDs=smp1,smp2,smp3,...
                std::vector<std::string> h; ngslib::split(s.s, h, "=");

                // `h[1]` is a string of sample ids, like: 'smp1,smp2,smp3,smp4,...'
                // set 'is_append' to be 'true', keep pushing data back in 'bf_smp_ids' vector.
                ngslib::split(h[1], bf_smp_ids, ",", true); 
                break;  // complete fetching the sample ids, end the loop
            }
        }
        free(s.s);
        if ((bgzf_close(f)) < 0) throw std::runtime_error("[ERROR] " + batchfiles[i] + " fail close.");
    }

    return bf_smp_ids;
}

bool _basevar_caller(const std::vector<std::string> &smp_bf_line_vector, 
                     const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                     double min_af, 
                     size_t n_sample,
                     BGZF *vcf_hd, 
                     BGZF *cvg_hd) 
{
    // Initial
    BatchInfo samples_bi;

    // Prepare the capacity for saving time
    samples_bi.align_bases.reserve(n_sample);
    samples_bi.align_base_quals.reserve(n_sample);
    samples_bi.mapqs.reserve(n_sample);
    samples_bi.map_strands.reserve(n_sample);
    samples_bi.base_pos_ranks.reserve(n_sample);
    samples_bi.n     = n_sample;
    samples_bi.depth = 0;

    std::vector<std::string> col_info;
    bool keep_push_back = true;
    for (size_t i(0); i < smp_bf_line_vector.size(); ++i) {
        // smp_bf_line_vector[i] format is: 
        // [CHROM, POS, REF, Depth, MappingQuality, Readbases, ReadbasesQuality, ReadPositionRank, Strand]
        // Looks like: chr11   5246595   C   10 37,0,0,0,0    ...
        ngslib::split(smp_bf_line_vector[i], col_info, "\t");
        if (col_info.size() != 9) {  // 9 is the column number of batchfile 
            throw std::runtime_error("[ERROR] batchfile has invalid data:\n" + smp_bf_line_vector[i]);
        }

        if (i == 0) {
            samples_bi.ref_id   = col_info[0];             // chromosome id
            samples_bi.ref_pos  = std::stoi(col_info[1]);  // string to int type
            samples_bi.ref_base = col_info[2];
            
        } else if ((samples_bi.ref_id   != col_info[0])            || 
                   (samples_bi.ref_pos  != std::stoi(col_info[1])) || 
                   (samples_bi.ref_base != col_info[2])) 
        {
            throw std::runtime_error("[ERROR] Batchfiles must have the same genome coordinate in each line.");
        }
         
        samples_bi.depth += std::stoi(col_info[3]);
        ngslib::split(col_info[4], samples_bi.mapqs,            " ", keep_push_back);  // keep adding data to vector
        ngslib::split(col_info[5], samples_bi.align_bases,      " ", keep_push_back);
        ngslib::split(col_info[6], samples_bi.align_base_quals, " ", keep_push_back);
        ngslib::split(col_info[7], samples_bi.base_pos_ranks,   " ", keep_push_back);
        ngslib::split(col_info[8], samples_bi.map_strands,      " ", keep_push_back);
    }

    // Coverage is 0 of all samples on this position. skip
    if (samples_bi.depth == 0) return false;
    
    // check data
    if ((samples_bi.mapqs.size()            != n_sample) || 
        (samples_bi.align_bases.size()      != n_sample) || 
        (samples_bi.align_base_quals.size() != n_sample) || 
        (samples_bi.map_strands.size()      != n_sample) || 
        (samples_bi.base_pos_ranks.size()   != n_sample))
    {
        std::cerr << "Total samples size is: " << n_sample            << "\n" + samples_bi.ref_id    << " " 
                  << samples_bi.ref_pos << " " << samples_bi.ref_base << " " << samples_bi.depth     << "\n" 
                  << "bt.samples_bi.mapqs.size():            " << samples_bi.mapqs.size()            << "\t" << ngslib::join(samples_bi.mapqs, " ")            << "\n"
                  << "bt.samples_bi.align_bases.size():      " << samples_bi.align_bases.size()      << "\t" << ngslib::join(samples_bi.align_bases, " ")      << "\n"
                  << "bt.samples_bi.align_base_quals.size(): " << samples_bi.align_base_quals.size() << "\t" << ngslib::join(samples_bi.align_base_quals, " ") << "\n"
                  << "bt.samples_bi.map_strands.size():      " << samples_bi.map_strands.size()      << "\t" << ngslib::join(samples_bi.map_strands, " ")      << "\n"
                  << "bt.samples_bi.base_pos_ranks.size():   " << samples_bi.base_pos_ranks.size()   << "\t" << ngslib::join(samples_bi.base_pos_ranks, " ")   << "\n"
                  << std::endl;
        throw std::runtime_error("[ERROR] Something is wrong in batchfiles.");
    }
    
    // output coverage first
    _out_cvg_line(&samples_bi, group_smp_idx, cvg_hd);

    // Detect variant and output
    BaseType bt(&samples_bi, min_af);
    bt.lrt();

    if (bt.get_alt_bases().size() && bt.is_only_snp()) { // only output SNPs in VCF file
        
        std::map<std::string, BaseType> popgroup_bt;
        if (!group_smp_idx.empty()) { // group is not empty

            std::vector<char> basecombination;
            basecombination.push_back(toupper(bt.get_ref_base()[0]));  // push upper reference base. 
            // do not sort, keep the original order as the same as 'alt_bases'
            basecombination.insert(basecombination.end(), bt.get_alt_bases().begin(), bt.get_alt_bases().end());
            
            // Call BaseType for each group
            std::map<std::string, std::vector<size_t>>::const_iterator it = group_smp_idx.begin();
            for (; it != group_smp_idx.end(); ++it) {
                popgroup_bt[it->first] = __gb(&samples_bi, it->second, basecombination, min_af);
            }
        }
        _out_vcf_line(bt, popgroup_bt, &samples_bi, vcf_hd);
    }

    return !bt.get_alt_bases().empty() && bt.is_only_snp();  // has SNP
}

const BaseType __gb(const BatchInfo *smp_bi, 
                    const std::vector<size_t> group_idx, 
                    const std::vector<char> &basecombination, 
                    double min_af) 
{
    BatchInfo g_smp_bi = __get_group_batchinfo(smp_bi, group_idx);
    BaseType bt(&g_smp_bi, min_af);
    bt.lrt(basecombination);

    return bt;
}

const BatchInfo __get_group_batchinfo(const BatchInfo *smp_bi, const std::vector<size_t> group_idx) {

    BatchInfo g_smp_bi;
    g_smp_bi.ref_id   = smp_bi->ref_id;
    g_smp_bi.ref_pos  = smp_bi->ref_pos;
    g_smp_bi.ref_base = smp_bi->ref_base;
    g_smp_bi.depth    = smp_bi->depth;    // 这个深度肯定不是 group samples 的深度，但没关系，不会用到
    g_smp_bi.n        = group_idx.size(); // sample size for speific group

    for (auto i : group_idx) {
        g_smp_bi.align_bases.push_back(smp_bi->align_bases[i]);
        g_smp_bi.align_base_quals.push_back(smp_bi->align_base_quals[i]);
        g_smp_bi.mapqs.push_back(smp_bi->mapqs[i]);
        g_smp_bi.map_strands.push_back(smp_bi->map_strands[i]);
        g_smp_bi.base_pos_ranks.push_back(smp_bi->base_pos_ranks[i]);
    }

    return g_smp_bi;
}


bool __create_a_batchfile(const std::vector<std::string> batch_align_files,  
                          const std::vector<std::string> batch_sample_ids,   
                          const std::string &fa_seq,                         
                          const ngslib::GenomeRegionTuple &genome_region,  // [chr, start, end]
                          const int mapq_thd,                              // mapping quality threshold
                          const std::string output_batch_file)             // output batchfile name
{   // 原为 BaseTypeRunner 的成员函数，未掌握如何将该函数指针传入 ThreadPool，遂作罢，后再改。
    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);

    // This value affected the computing memory, could be set larger than 500000, 20 just for test
    // 这个参数是为了限制存入 `batchsamples_posinfomap_vector` 的最大读取区间，从而控制内存消耗不要太大
    static const uint32_t STEP_REGION_LEN = 500000;  

    std::string ref_id; uint32_t reg_beg, reg_end;
    std::tie(ref_id, reg_beg, reg_end) = genome_region;  // 1-based

    BGZF *obf = bgzf_open(output_batch_file.c_str(), "w"); // output file handle of output_batch_file
    if (!obf) throw std::runtime_error("[ERROR] " + output_batch_file + " open failure.");

    // Header of batchfile
    std::string bf_header = "##fileformat=BaseVarBatchFile_v1.0\n" 
                            "##SampleIDs=" + ngslib::join(batch_sample_ids, ",") + "\n" + 
                            "#CHROM\tPOS\tREF\tDepth(CoveredSample)\tMappingQuality\t"
                            "Readbases\tReadbasesQuality\tReadPositionRank\tStrand\n";
    if (bgzf_write(obf, bf_header.c_str(), bf_header.length()) != bf_header.length())
        throw std::runtime_error("[ERROR] fail to write data");

    PosMapVector batchsamples_posinfomap_vector;
    batchsamples_posinfomap_vector.reserve(batch_align_files.size());  //  pre-set the capacity

    bool is_empty = true, has_data = false;
    uint32_t sub_reg_beg, sub_reg_end;
    for (uint32_t i(reg_beg), j(0); i < reg_end + 1; i += STEP_REGION_LEN, ++j) {
        // Cut smaller regions to save computing memory.
        sub_reg_beg = i;
        sub_reg_end = sub_reg_beg + STEP_REGION_LEN - 1 > reg_end ? reg_end : sub_reg_beg + STEP_REGION_LEN - 1;
        is_empty = __fetch_base_in_region(batch_align_files, fa_seq, mapq_thd, 
                                          std::make_tuple(ref_id, sub_reg_beg, sub_reg_end),
                                          batchsamples_posinfomap_vector);  // 传引用，省内存，得数据

        if (!has_data && !is_empty) {
            has_data = true;
        }

        /* Output batchfile, no matter 'batchsamples_posinfomap_vector' is empty or not. */
        __write_record_to_batchfile(batchsamples_posinfomap_vector, fa_seq,
                                    std::make_tuple(ref_id, sub_reg_beg, sub_reg_end), 
                                    obf);
        batchsamples_posinfomap_vector.clear();  // 必须清空，为下个循环做准备
    }

    int is_cl = bgzf_close(obf);  // close file
    if (is_cl < 0) {
        throw std::runtime_error("[ERROR] " + output_batch_file + " fail close.");
    }

    // Create a Tabix index for 'output_batch_file'
    // conf: {preset, seq col, beg col, end col, header-char, skip-line}
    const tbx_conf_t bf_tbx_conf = {1, 1, 2, 0, '#', 0};
    if (tbx_index_build(output_batch_file.c_str(), 0, &bf_tbx_conf))  // file suffix is ".tbi"
        throw std::runtime_error("tbx_index_build failed: Is the file bgzip-compressed? "
                                 "Check this file: " + output_batch_file + "\n");

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now)); 
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "[INFO] " + ct + ". Done for creating batchfile " 
              << output_batch_file << ", " << difftime(now, real_start_time) 
              << " (CPU time: " << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
              << ") seconds elapsed." << std::endl;

    return has_data;
}

bool __fetch_base_in_region(const std::vector<std::string> &batch_align_files,  
                            const std::string &fa_seq,  // must be the whole chromosome sequence  
                            const int mapq_thd,         // mapping quality threshold
                            const ngslib::GenomeRegionTuple genome_region,
                            PosMapVector &batchsamples_posinfomap_vector)  
{
    // only using here, In case of missing the overlap reads on side position, 200bp would be enough
    static const uint32_t REG_EXPEND_SIZE = 200;

    std::string ref_id; uint32_t reg_start, reg_end;
    std::tie(ref_id, reg_start, reg_end) = genome_region;  // 1-based

    uint32_t exp_reg_start = reg_start > REG_EXPEND_SIZE ? reg_start - REG_EXPEND_SIZE : 1;
    uint32_t exp_reg_end   = reg_end + REG_EXPEND_SIZE;
    std::string exp_regstr = ref_id + ":" + ngslib::tostring(exp_reg_start) + "-" + ngslib::tostring(exp_reg_end);

    // Loop all alignment files
    bool is_empty = true;
    for(size_t i(0); i < batch_align_files.size(); ++i) {
        ngslib::Bam bf(batch_align_files[i], "r");  // open bamfile in reading mode (one sample, one bamfile)

        // 位点信息存入该变量, 且由于是按区间读取比对数据，key 值无需再包含 ref_id，因为已经不言自明。
        PosMap sample_posinfo_map;

        if (bf.fetch(exp_regstr)) { // Set 'bf' only fetch alignment reads in 'exp_regstr'.
            hts_pos_t map_ref_start, map_ref_end;  // hts_pos_t is uint64_t
            std::vector<ngslib::BamRecord> sample_target_reads; 
            ngslib::BamRecord al;       // alignment read

            while (bf.next(al) >= 0) {  // -1 => hit the end of alignement file.
                if (al.mapq() < mapq_thd || al.is_duplicate() || al.is_qc_fail()) continue;
                map_ref_start = al.map_ref_start_pos() + 1;  // al.map_ref_start_pos() is 0-based, convert to 1-based
                map_ref_end   = al.map_ref_end_pos();        // al.map_ref_end_pos() is 1-based

                // Only fetch reads which in [reg_start, reg_end]
                if (reg_start > map_ref_end) continue;
                if (reg_end < map_ref_start) break;

                sample_target_reads.push_back(al);  // record the proper reads of sample
            }

            sample_posinfo_map.clear();  // make sure it's empty
            if (sample_target_reads.size() > 0) {
                // get alignment information of [i] sample.
                __seek_position(sample_target_reads, fa_seq, genome_region, sample_posinfo_map);
            }
        }

        if (is_empty && !sample_posinfo_map.empty()) { 
            // at least one sample has data in this region
            is_empty = false; 
        }

        // Push it into 'batchsamples_posinfomap_vector' even if 'sample_posinfo_map' is empty, 
        // make sure 'batchsamples_posinfomap_vector' has the same size as `batch_align_files`
        batchsamples_posinfomap_vector.push_back(sample_posinfo_map);
    }

    if (batchsamples_posinfomap_vector.size() != batch_align_files.size())
        throw std::runtime_error("[basetype.cpp::__fetch_base_in_region] 'pos_batchinfo_vector.size()' "
                                 "should be the same as 'batch_align_files.size()'");

    return is_empty;  // no cover reads in 'genome_region' if empty.
}

void __seek_position(const std::vector<ngslib::BamRecord> &sample_map_reads,
                     const std::string &fa_seq,   // must be the whole chromosome sequence
                     const ngslib::GenomeRegionTuple genome_region,
                     PosMap &sample_posinfo_map)
{
    if (!sample_posinfo_map.empty())
        throw std::runtime_error("[basetype.cpp::__seek_position] 'sample_posinfo_map' must be empty.");

    std::string ref_id; uint32_t reg_start, reg_end;
    std::tie(ref_id, reg_start, reg_end) = genome_region;  // 1-based

    AlignBaseInfo align_base_info;
    align_base_info.ref_id = ref_id;

    // A vector of: (cigar_op, read position, reference position, read base, read_qual, reference base)
    std::vector<ngslib::ReadAlignedPair> aligned_pairs;
    for(size_t i(0); i < sample_map_reads.size(); ++i) {

        align_base_info.map_strand = sample_map_reads[i].map_strand();  // '*', '-' or '+'
        align_base_info.mapq = sample_map_reads[i].mapq();

        aligned_pairs = sample_map_reads[i].get_aligned_pairs(fa_seq);
        char mean_qqual_char = int(sample_map_reads[i].mean_qqual()) + 33; // 33 is the offset of base QUAL
        uint32_t map_ref_pos;
        for (size_t i(0); i < aligned_pairs.size(); ++i) {
            // Todo: data of 'align_base_info' and 'aligned_pairs[i]' is similar, 
            // just use 'aligned_pairs[i]' to replace 'align_base_info'?
            map_ref_pos = aligned_pairs[i].ref_pos + 1;  // ref_pos is 0-based, convert to 1-based;

            if (reg_end < map_ref_pos) break;
            if (reg_start > map_ref_pos) continue;

            // 'BAM_XXX' are macros for CIGAR, which defined in 'htslib/sam.h'
            if (aligned_pairs[i].op == BAM_CMATCH ||  /* CIGAR: M */ 
                aligned_pairs[i].op == BAM_CEQUAL ||  /* CIGAR: = */
                aligned_pairs[i].op == BAM_CDIFF)     /* CIGAR: X */
            {
                // One character
                align_base_info.ref_base       = aligned_pairs[i].ref_base[0];
                align_base_info.read_base      = aligned_pairs[i].read_base[0];
                align_base_info.read_base_qual = aligned_pairs[i].read_qual[0];
            } else if (aligned_pairs[i].op == BAM_CINS) {  /* CIGAR: I */
                if (!aligned_pairs[i].ref_base.empty()) {
                    std::cerr << sample_map_reads[i] << "\n";
                    throw std::runtime_error("[ERROR] We got reference base in insertion region.");
                }

                // roll back one position to the left side of insertion break point.
                --map_ref_pos;
                align_base_info.ref_base       = fa_seq[aligned_pairs[i].ref_pos-1]; // break point's ref base
                align_base_info.read_base      = fa_seq[aligned_pairs[i].ref_pos-1] + aligned_pairs[i].read_base;
                align_base_info.read_base_qual = mean_qqual_char;  // set to be mean quality of the whole read
            } else if (aligned_pairs[i].op == BAM_CDEL) {  /* CIGAR: D */
                if (!aligned_pairs[i].read_base.empty()) {
                    std::cerr << sample_map_reads[i] << "\n";
                    throw std::runtime_error("[ERROR] We got read bases in deletion region.");
                }

                // roll back one position to the left side of deletion break point.
                --map_ref_pos;
                align_base_info.ref_base       = fa_seq[aligned_pairs[i].ref_pos-1] + aligned_pairs[i].ref_base;
                align_base_info.read_base      = fa_seq[aligned_pairs[i].ref_pos-1]; // break point's ref base
                align_base_info.read_base_qual = mean_qqual_char;  // set to be mean quality of the whole read
            } else { 
                // Skip in basevar for other kind of CIGAR symbals.
                continue;
            }
            align_base_info.ref_pos = map_ref_pos;

            // qpos is 0-based, conver to 1-based to set the rank of base on read.
            align_base_info.rpr = aligned_pairs[i].qpos + 1;

            if (sample_posinfo_map.find(map_ref_pos) == sample_posinfo_map.end()) {
                // Just get the base from first read which aligned on this ref_pos,
                // no matter the first one it's indel or not.

                // {ref_pos (no need to add ref_id in the key) => map info}
                sample_posinfo_map.insert({map_ref_pos, align_base_info});
            }
        }
    }

    return;
}

// Create batch file for variant discovery
void __write_record_to_batchfile(const PosMapVector &batchsamples_posinfomap_vector,
                                 const std::string &fa_seq,
                                 const ngslib::GenomeRegionTuple genome_region, 
                                 BGZF *obf) 
{
    const static char BASE_Q0_ASCII = '!';  // The ascii code of '!' character is 33

    std::string ref_id; uint32_t reg_start, reg_end;
    std::tie(ref_id, reg_start, reg_end) = genome_region;  // 1-based

    // Output columns and set zhe capacity for each vector: 
    // [CHROM, POS, REF, Depth(CoveredSample), MappingQuality, 
    //  Readbases, ReadbasesQuality, ReadPositionRank, Strand]
    size_t sn = batchsamples_posinfomap_vector.size();
    std::vector<int> mapq;                     mapq.reserve(sn);
    std::vector<std::string> map_read_bases;   map_read_bases.reserve(sn);
    std::vector<char> map_read_base_qualities; map_read_base_qualities.reserve(sn);
    std::vector<int> read_pos_ranks;           read_pos_ranks.reserve(sn);
    std::vector<char> map_strands;             map_strands.reserve(sn);

    for (uint32_t pos(reg_start); pos < reg_end+1; ++pos) {

        PosMap::const_iterator smp_pos_it;  // specifi sample position map
        uint32_t depth = 0;
        for (size_t i = 0; i < sn; i++) {
            smp_pos_it = batchsamples_posinfomap_vector[i].find(pos);
            if (smp_pos_it != batchsamples_posinfomap_vector[i].end()) {

                ++depth;
                if (smp_pos_it->second.ref_id != ref_id || smp_pos_it->second.ref_pos != pos)
                    throw std::runtime_error("[ERROR] reference id or position not match.");
                
                mapq.push_back(smp_pos_it->second.mapq);
                if (smp_pos_it->second.ref_base.size() == smp_pos_it->second.read_base.size()) { // SNV
                    map_read_bases.push_back(smp_pos_it->second.read_base);
                } else if (smp_pos_it->second.ref_base.size() < smp_pos_it->second.read_base.size()) {  // INS
                    map_read_bases.push_back("+" + smp_pos_it->second.read_base);
                } else { // DEL
                    map_read_bases.push_back("-" + smp_pos_it->second.ref_base);
                }
                map_read_base_qualities.push_back(smp_pos_it->second.read_base_qual);
                read_pos_ranks.push_back(smp_pos_it->second.rpr);
                map_strands.push_back(smp_pos_it->second.map_strand);

            } else {
                mapq.push_back(0);
                map_read_bases.push_back("N");
                map_read_base_qualities.push_back(BASE_Q0_ASCII);
                read_pos_ranks.push_back(0);
                map_strands.push_back('.');
            }
        }

        std::string out = ref_id + "\t" + std::to_string(pos) + "\t" + fa_seq[pos-1] + "\t" + 
                          std::to_string(depth)                      + "\t" + 
                          ngslib::join(mapq, " ")                    + "\t" + 
                          ngslib::join(map_read_bases, " ")          + "\t" + 
                          ngslib::join(map_read_base_qualities, " ") + "\t" +
                          ngslib::join(read_pos_ranks, " ")          + "\t" + 
                          ngslib::join(map_strands, " ")             + "\n";
        
        // write to file and check is successful or not.
        if (bgzf_write(obf, out.c_str(), out.length()) != out.length())
            throw std::runtime_error("[ERROR] fail to write data");

        // clear up 
        mapq.clear();
        map_read_bases.clear();
        map_read_base_qualities.clear();
        read_pos_ranks.clear();
        map_strands.clear();
    }

    return;
}

void _out_vcf_line(const BaseType &bt, 
                   const std::map<std::string, BaseType> &group_bt, 
                   const BatchInfo *smp_bi, 
                   BGZF *vcf_hd) 
{

    std::map<char, std::string> alt_gt;
    std::vector<int> cm_ac; 
    std::vector<double> cm_af, cm_caf;
    double ad_sum = 0;
    for (size_t i(0); i < bt.get_alt_bases().size(); ++i) {

        char b = bt.get_alt_bases()[i];
        alt_gt[b] = "./" + std::to_string(i+1);

        ad_sum = ad_sum + bt.get_base_depth(b);
        cm_ac.push_back(bt.get_base_depth(b));
        cm_af.push_back(bt.get_lrt_af(b));
        cm_caf.push_back(bt.get_base_depth(b) / bt.get_total_depth());  // Allele frequency by read count
    }

    // Sample genotype information
    std::vector<std::string> samples;
    std::vector<char> align_bases;  // used in ranksumtest
    std::string gt;
    char upper_ref_base = toupper(bt.get_ref_base()[0]);  // only for SNPs
    for (size_t i(0); i < smp_bi->n; ++i) {

        char fb = smp_bi->align_bases[i][0];  // a string, get first base
        align_bases.push_back(fb);            // only get the first base, suit for SNPs

        if (fb != 'N' && fb != '+' && fb != '-') {
            // For the base which not in bt.get_alt_bases()
            if (alt_gt.find(fb) == alt_gt.end()) alt_gt[fb] = "./.";
            gt = (fb == upper_ref_base) ? "0/." : alt_gt[fb];

            double qual = bt.get_qual_pvalue()[i];
            // GT:AB:SO:BP
            samples.push_back(gt + ":" + fb + ":" + smp_bi->map_strands[i] + ":" + std::to_string(qual));
        } else {
            samples.push_back("./.");  // 'N' or Indel.
        }
    }

    // For INFO
    std::string alt_bases_string = ngslib::join(bt.get_alt_bases(), "");

    // Rank Sum Test for mapping qualities of REF versus ALT reads
    int mq_rank_sum = ref_vs_alt_ranksumtest(upper_ref_base, alt_bases_string, align_bases, smp_bi->mapqs);

    // Rank Sum Test for variant appear position among read of REF versus ALT
    int read_pos_rank_sum = ref_vs_alt_ranksumtest(upper_ref_base, alt_bases_string, align_bases, smp_bi->base_pos_ranks);
    
    // Rank Sum Test for base quality of REF versus ALT
    int base_q_rank_sum = ref_vs_alt_ranksumtest(upper_ref_base, alt_bases_string, align_bases, smp_bi->align_base_quals);

    // Variant call confidence normalized by depth of sample reads supporting a variant.
    double qd = bt.get_var_qual() / ad_sum; 
    if (qd == 0) qd = 0.0; // make sure -0.0 => 0.0

    // Strand bias by fisher exact test and Strand bias estimated by the Symmetric Odds Ratio test.
    StrandBiasInfo sbi = strand_bias(upper_ref_base, alt_bases_string, align_bases, smp_bi->map_strands);

    // construct the INFO field.
    std::vector<std::string> info = {
        "CM_DP="  + std::to_string(bt.get_total_depth()),
        "CM_AC="  + ngslib::join(cm_ac, ","),
        "CM_AF="  + ngslib::join(cm_af, ","), 
        "CM_CAF=" + ngslib::join(cm_caf, ","),
        "MQRankSum="      + std::to_string(mq_rank_sum),
        "ReadPosRankSum=" + std::to_string(read_pos_rank_sum),
        "BaseQRankSum="   + std::to_string(base_q_rank_sum),
        "QD="     + std::to_string(qd),  // use 'roundf(qd * 1000) / 1000' to round the value to 3 decimal places
        "SOR="    + std::to_string(sbi.sor),
        "FS="     + std::to_string(sbi.fs), 
        "SB_REF=" + std::to_string(sbi.ref_fwd) + "," + std::to_string(sbi.ref_rev),
        "SB_ALT=" + std::to_string(sbi.alt_fwd) + "," + std::to_string(sbi.alt_rev)
    };

    if (!group_bt.empty()) { // not empty
        std::vector<std::string> group_af_info;
        for (std::map<std::string, BaseType>::const_iterator it(group_bt.begin()); it != group_bt.end(); ++it) {

            std::vector<double> af;
            char b;
            for (auto b : it->second.get_alt_bases()) {
                af.push_back(it->second.get_lrt_af(b));
            }

            if (!af.empty())  // has group AF
                group_af_info.push_back(it->first + "_AF=" + ngslib::join(af, ",")); // groupID_AF=xxx,xxx
        }
        info.insert(info.end(), group_af_info.begin(), group_af_info.end());
    }

    std::string sample_format = "GT:AB:SO:BP";
    std::string qs  = (bt.get_var_qual() > QUAL_THRESHOLD) ? "." : "LowQual";
    std::string out = bt.get_ref_id() + "\t" + std::to_string(bt.get_ref_pos()) + "\t.\t" + bt.get_ref_base() + "\t" + 
                      ngslib::join(bt.get_alt_bases(), ",") + "\t" + std::to_string(bt.get_var_qual()) + "\t" + 
                      qs + "\t" + ngslib::join(info, ";")   + "\t" + sample_format + "\t" + 
                      ngslib::join(samples, "\t") + "\n";
    // write to file and check is successful or not.
    if (bgzf_write(vcf_hd, out.c_str(), out.length()) != out.length())
        throw std::runtime_error("[ERROR] fail to write data");

    return;
}

void _out_cvg_line(const BatchInfo *smp_bi, 
                   const std::map<std::string, std::vector<size_t>> & group_smp_idx, 
                   BGZF *cvg_hd) 
{
    // base depth and indels for each subgroup
    std::map<std::string, IndelTuple> group_cvg;
    // Call BaseType for each group
    std::map<std::string, std::vector<size_t>>::const_iterator it = group_smp_idx.begin();
    for (; it != group_smp_idx.end(); ++it) {
        const BatchInfo g_smp_bi = __get_group_batchinfo(smp_bi, it->second);
        group_cvg[it->first] = __base_depth_and_indel(g_smp_bi.align_bases);
    } // group_cvg 信息计算了，因为有 bug 未用上（2024-09-05）

    std::map<char, int> base_depth;
    std::string indel_string;
    // coverage info for each position
    std::tie(indel_string, base_depth) = __base_depth_and_indel(smp_bi->align_bases);

    std::vector<char> align_bases;  // used in ranksumtest
    for (size_t i(0); i < smp_bi->n; ++i) {
        char fb = smp_bi->align_bases[i][0];  // a string, get first base
        align_bases.push_back(fb); 
    }

    char upper_ref_base = toupper(smp_bi->ref_base[0]);  // only for SNPs
    std::vector<char> alt_bases;  // used in ranksumtest
    int total_depth = 0;
    for (std::map<char, int>::iterator it(base_depth.begin()); it != base_depth.end(); ++it) {
        total_depth += it->second;
        if (it->first != upper_ref_base) {
            alt_bases.push_back(it->first);
        }
    }

    StrandBiasInfo sbi = strand_bias(upper_ref_base, ngslib::join(alt_bases, ""), align_bases, smp_bi->map_strands);
    if (total_depth > 0) {
        std::vector<int> dd;
        for (auto b : BASES) dd.push_back(base_depth[b]);
        std::string out = smp_bi->ref_id              + "\t" + std::to_string(smp_bi->ref_pos) + "\t" + 
                          smp_bi->ref_base            + "\t" + std::to_string(total_depth)     + "\t" + 
                          ngslib::join(dd, "\t")      + "\t" + indel_string                    + "\t" +
                          std::to_string(sbi.fs)      + "\t" + std::to_string(sbi.sor)         + "\t" +
                          std::to_string(sbi.ref_fwd) + ","  + std::to_string(sbi.ref_rev)     + ","  +
                          std::to_string(sbi.alt_fwd) + ","  + std::to_string(sbi.alt_rev)     + "\n";
        if (bgzf_write(cvg_hd, out.c_str(), out.length()) != out.length())
            throw std::runtime_error("[ERROR] fail to write data");
    }

    return;
}


IndelTuple __base_depth_and_indel(const std::vector<std::string> &align_bases) {
// bug: indel 没被正常输出，而是被跳过了 (原因是在 batchfile 时，优先输出了 snp，而不是 indels 怎被设为 0 覆盖)
    std::map<char, int> base_depth;
    std::string indel_string;

    for (auto b : BASES) base_depth[b] = 0;

    std::map<std::string, int> indel_depth;
    for (auto bs: align_bases) {
        if (bs[0] == 'N') continue;

        if (base_depth.find(bs[0]) != base_depth.end()) {
            // ignore all bases('*') which not match ``BASE``
            base_depth[bs[0]]++;
        } else {
            // indel
            indel_depth[bs]++;
        }
    }

    std::vector<std::string> indels;
    for (std::map<std::string, int>::iterator it(indel_depth.begin()); it != indel_depth.end(); ++it) {
        indels.push_back(it->first + "|" + std::to_string(it->second));
    }
    indel_string = (!indels.empty()) ? ngslib::join(indels, ",") : ".";
    return std::make_tuple(indel_string, base_depth);
}
