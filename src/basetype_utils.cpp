#include <iostream>
#include <fstream>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "basetype_utils.h"
#include "fasta.h"
#include "utils.h"

std::vector<std::string> get_firstcolumn_from_file(const std::string fn) {

    std::ifstream i_fn(fn.c_str());
    if (!i_fn) {
        std::cerr << "[ERROR] Cannot open file: " + fn << std::endl;
        exit(1);
    }

    std::vector<std::string> first_column;
    std::string skip, first_col_str;
    while (1) {
        i_fn >> first_col_str;
        if (i_fn.eof()) break;
        std::getline(i_fn, skip, '\n');  // skip the rest information of line.

        first_column.push_back(first_col_str);
    }
    i_fn.close();

    return first_column;
}

std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &addition_info, 
                              const std::vector<std::string> &samples)
{
    std::vector<std::string> header = {
        "##fileformat=VCFv4.2",
        "##FILTER=<ID=LowQual,Description=\"Low quality (QUAL < 60)\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=AB,Number=1,Type=String,Description=\"Allele Base\">",
        "##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Strand orientation of the mapping base. Marked as + or -\">",
        "##FORMAT=<ID=BP,Number=1,Type=String,Description=\"Base Probability which calculate by base quality\">",
        "##INFO=<ID=CM_AF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies base on LRT algorithm\">",
        "##INFO=<ID=CM_CAF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies just base on read count\">",
        "##INFO=<ID=CM_AC,Number=A,Type=Integer,Description=\"An ordered, comma delimited allele depth in CMDB\">",
        "##INFO=<ID=CM_DP,Number=A,Type=Integer,Description=\"Total Depth in CMDB\">",
        "##INFO=<ID=SB_REF,Number=A,Type=Integer,Description=\"Read number support REF: Forward,Reverse\">",
        "##INFO=<ID=SB_ALT,Number=A,Type=Integer,Description=\"Read number support ALT: Forward,Reverse\">",
        "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">",
        "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Phred-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">",
        "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">",
        "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Phred-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">",
        "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Phred-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">",
        "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence Quality by Depth\">"
    };  // initial by common information of header
    if (!addition_info.empty()) header.insert(header.end(), addition_info.begin(), addition_info.end());

    ngslib::Fasta fa = ref_file_path;
    std::vector<std::string> contigs;
    for (size_t i(0); i < fa.nseq(); ++i) {
        std::string seqname = fa.iseq_name(i);
        uint32_t seqlen = fa.seq_length(seqname);
        contigs.push_back("##contig=<ID=" + seqname + ",length=" + std::to_string(seqlen) + 
                          ",assembly=" + ref_file_path + ">");
    }
    header.insert(header.end(), contigs.begin(), contigs.end());

    header.push_back("##reference=file://" + ngslib::abspath(ref_file_path));
    header.push_back("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + ngslib::join(samples, "\t"));

    return ngslib::join(header, "\n");
}

std::string cvg_header_define(const std::vector<std::string> &group_info, const std::vector<char> &BASES) {

    std::string h = "#CHROM\tPOS\tREF\tDepth\t" + ngslib::join(BASES, "\t") + "\t" +
                    "Indels\tFS\tSOR\tStrand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)";
    // group cvg 的计算有 bug， 暂时不输出
    // if (!group_info.empty())
    //     h += "\t" + ngslib::join(group_info, "\t");

    std::vector<std::string> header = {
        "##fileformat=CVGv1.0",
        "##Group information is the depth of A:C:G:T:Indel", 
        h
    };

    return ngslib::join(header, "\n");
}

void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header, bool is_remove_tempfile) 
{
    if (infiles.empty()) return;

    bool is_compress_out = (ngslib::suffix_name(outfile) == ".gz") ? true : false;
    BGZF *OUT = bgzf_open(outfile.c_str(), is_compress_out ? "w" : "uw");  // output file
    if (!OUT) throw std::runtime_error("[ERROR] " + outfile + " open failure.");

    header += "\n";
    if (bgzf_write(OUT, header.c_str(), header.length()) != header.length())
        throw std::runtime_error("[ERROR] fail to write data");

    /* Merge all files here */
    for (auto fn: infiles) {
        BGZF *f = bgzf_open(fn.c_str(), "r");
        kstring_t s; s.s = NULL; s.l = s.m = 0;
        while (bgzf_getline(f, '\n', &s) >= 0) {
            if (s.s[0] == '#') continue;  // ignore the header of subfiles.
            std::string out(s.s); out += "\n";

            if (bgzf_write(OUT, out.c_str(), out.length()) != out.length())
                throw std::runtime_error("[ERROR] fail to write data");
        }
        bgzf_close(f);

        if (is_remove_tempfile) ngslib::safe_remove(fn);
    }
    
    int is_cl = bgzf_close(OUT);
    if (is_cl < 0) throw std::runtime_error("[ERROR] " + outfile + " fail close.");
    
    return;
}

