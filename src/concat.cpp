#include "concat.h"

int _concat_basevar_outfile(const std::vector<std::string> &infiles, const std::string outfile) {

    if (infiles.empty()) return 1;

    BGZF *f = bgzf_open(infiles[0].c_str(), "r");
    if (!f) throw std::runtime_error("[ERROR] " + infiles[0] + " open failure.");

    // Get header information from the first file would be enough.
    std::vector<std::string> h;
    kstring_t s; s.s = NULL; s.l = s.m = 0;
    while (bgzf_getline(f, '\n', &s) >= 0) {
        if (s.s[0] != '#') { // head char is '#'
            break;
        }
        h.push_back(s.s);
    }
    free(s.s);
    if ((bgzf_close(f)) < 0) throw std::runtime_error("[ERROR] " + infiles[0] + " fail close.");

    merge_file_by_line(infiles, outfile, ngslib::join(h, "\n"));

    return 0;
}

int concat_runner(int argc, char *argv[]) {

    if (argc < 2) {
        std::cout << __CONCAT_USAGE << "\n" << std::endl;
        exit(1);
    }

    std::vector<std::string> input_files;
    std::string in_filelist, output_file;
    // Parsing the commandline options. 
    char c;
    while((c = getopt_long(argc, argv, "I:L:O:h", BASETYPE_CMDLINE_LOPTS, NULL)) >= 0) {
        // 字符流解决命令行参数转浮点等类型的问题
        std::stringstream ss(optarg ? optarg: "");  
        switch (c) {
            case 'I': input_files.push_back(optarg);   break;
            case 'L': in_filelist = optarg;            break;
            case 'O': output_file = optarg;            break;
            case 'h': std::cout << __CONCAT_USAGE << std::endl; exit(1);
            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                exit(1);
        }
    }

    /* Make sure we set valid arguments */
    if (input_files.empty() && in_filelist.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-I/--input' or '-L/--file-list'");
    if (output_file.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-O/--output'");
    
    if (!in_filelist.empty()) {
        std::vector<std::string> filelist = get_firstcolumn_from_file(in_filelist);
        input_files.insert(input_files.end(), filelist.begin(), filelist.end());
    }
    std::cout << "[INFO] Finish loading arguments and we have " << input_files.size()
              << " files to concat.\n" << std::endl;
    
    return _concat_basevar_outfile(input_files, output_file);
}

