/**
 * @file utils.cpp
 * @author Shujia Huang
 * @brief 
 * @date 2021-08-18
 * 
 */
#include <unistd.h>
#include <filesystem>  // Available since C++17

#include "utils.h"


namespace ngslib {

    // http://c.biancheng.net/cpp/html/303.html
    bool is_readable(const char *name) {
        return (access(name, R_OK) == 0);
    }

    std::string abspath(const std::string file_path) {
        std::filesystem::path fp(file_path);
        return std::filesystem::absolute(fp).string();
    }

    std::string dirname(const std::string file_path) {
        size_t p(file_path.find_last_of("/\\"));
        return p > 0 && p != std::string::npos ? file_path.substr(0, p) : file_path;
    }

    std::string basename(const std::string file_path) {
        return file_path.substr(file_path.find_last_of("/\\") + 1);
    }

    std::string suffix_name(const std::string file_path) {
        size_t si(file_path.find_last_of('.'));
        std::string suffix = (si != std::string::npos) ? file_path.substr(si) : file_path;
        return suffix;
    }
    
    std::string remove_filename_extension(const std::string filename) {
        size_t p(filename.find_last_of('.'));
        return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
    }

    bool path_exists_and_not_empty(std::string folder_path) {
        return std::filesystem::exists(folder_path) && !std::filesystem::is_empty(folder_path);
    }

    bool safe_mkdir(std::string folder_path) {
        // set the folder_path tobe 'filesystem::path'
        return std::filesystem::create_directories(folder_path);
    }

    bool safe_remove(std::string filepath) {
        return std::filesystem::remove(filepath);
    }

    std::string get_last_modification_file(std::string directory_path) {
        // Find the last modification file in a directory and return it.
        std::filesystem::path lmf;        
        std::filesystem::file_time_type lmf_time, cf_time;  

        bool flag(true);
        for (const auto &fn : std::filesystem::directory_iterator(directory_path)) {
            cf_time = std::filesystem::last_write_time(fn);
            if (flag) {
                flag = false;
                lmf = fn;
                lmf_time = cf_time;
            } else if (cf_time > lmf_time) {
                // find last modification file.
                lmf = fn;
                lmf_time = cf_time;
            }
        }

        return lmf.string();
    }

    void split(const std::string &in_str, std::vector<std::string> &out, const char *delim, bool is_append) {

        if (!is_append) { out.clear(); }

        size_t i(0), find_start_idx(0), delim_len(std::strlen(delim)), len;
        std::string item;
        while(i != std::string::npos) {

            i = in_str.find(delim, find_start_idx);
            len = (i==std::string::npos) ? in_str.length() - find_start_idx : i - find_start_idx;
            
            item = in_str.substr(find_start_idx, len);
            out.push_back(item);

            find_start_idx = i + delim_len;  // skip 'delim' character and set to find next 'delim' string
        }

        return;
    }

    std::vector<GenomeRegionTuple> region_slice(const GenomeRegionTuple &genome_region, int num) {

        if (num <= 0) {
            throw std::runtime_error("[ERROR] slice number must be bigger than 0.");
        }

        std::string ref_id; uint32_t reg_start, reg_end;
        std::tie(ref_id, reg_start, reg_end) = genome_region;  // 1-based

        size_t reg_len = reg_end - reg_start + 1;
        size_t delta = (reg_len % num) ? int(reg_len / num) + 1 : int(reg_len / num);

        std::vector<GenomeRegionTuple> regions;
        uint32_t sub_reg_beg, sub_reg_end;
        for (size_t i(reg_start); i < reg_end + 1; i += delta) {
            sub_reg_beg = i;
            sub_reg_end   = sub_reg_beg + delta - 1 > reg_end ? reg_end : sub_reg_beg + delta - 1;
            regions.push_back(std::make_tuple(ref_id, sub_reg_beg, sub_reg_end));
        }

        return regions;
    }
    
}  // namespae ngslib
