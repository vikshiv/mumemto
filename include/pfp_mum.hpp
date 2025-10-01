/*
 * File: pfp_mum.hpp
 * Description: Header file for pfp_mum, contains the
 *              function headers and other needed structs 
 *              for pfp_mum.cpp
 *              Heavily adapted from docprofiles by Omar Ahmed: https://github.com/oma219/docprofiles/tree/main
 * Date: December 20th, 2023
 */

#ifndef PFP_MUM_H
#define PFP_MUM_H

#include <string>
#include <iostream>
#include <filesystem>
#include <vector>
#include <ref_builder.hpp>

/* Useful MACROs */
#define FATAL_ERROR(...) do {std::fprintf(stderr, "\n\033[31mError: \033[m"); std::fprintf(stderr, __VA_ARGS__);\
                              std::fprintf(stderr, "\n\n"); std::exit(1);} while(0)
#define ASSERT(condition, msg) do {if (!condition){std::fprintf(stderr, "\n\n\033[31mAssertion Failed:\033[m %s\n\n", msg); \
                                                   std::exit(1);}} while(0)
#define STATUS_LOG(x, ...) do {std::fprintf(stderr, "\033[32m[%s] \033[0m", x); std::fprintf(stderr, __VA_ARGS__ ); \
                               std::fprintf(stderr, " ... ");} while(0)
#define DONE_LOG(x) do {auto sec = std::chrono::duration<double>(x); \
                        std::fprintf(stderr, "done.  (%.3f sec)\n", sec.count());} while(0)
#define FORCE_LOG(func, ...)  do {std::fprintf(stderr, "\033[32m[%s] \033[m", func); \
                                  std::fprintf(stderr, __VA_ARGS__); \
                                  std::fprintf(stderr, "\n");} while (0)

// Defintions
#define PFPMUM_VERSION "1.3.3"

#define MAXLCPVALUE 65535 // 2^16 - 1
#define MAXDOCS 65535

#define AVX2_PRESENT __AVX2__ 
#define AVX512BW_PRESENT __AVX512BW__ 


/* Function declations */
int mumemto_usage();
int build_main(int argc, char** argv);
int mumemto_short_usage();
int is_file(std::string path);
int is_dir(std::string path);
std::string  make_filelist(std::vector<std::string> files, std::string output_prefix);
void remove_temp_files(std::string filename);
std::vector<std::string> split(std::string input, char delim);
bool is_integer(const std::string& str);
bool endsWith(const std::string& str, const std::string& suffix);
std::string execute_cmd(const char* cmd);

struct BuildOptions {
    public:
        std::string input_list = "";
        std::string output_prefix = "output";
        std::string output_ref = "";
        std::vector<std::string> files;
        bool use_rcomp = true;
        size_t pfp_w = 10;
        size_t hash_mod = 100;
        // size_t threads = 0;
        bool arrays_out = false;
        std::string  arrays_in = "";
        bool arrays_in_flag = false;
        bool keep_temp = false;
        int num_distinct_docs = 0;
        bool overlap = true;
        std::string parse_prefix = "";
        bool from_parse_flag = false;
        size_t min_match_len = 20;
        int max_mem_freq = 0;
        int rare_freq = 1;
        bool binary = false;
        bool merge = false;
        bool anchor_merge = false;
        bool use_gsacak = false;
        bool only_parse = false;

        bool validate() {
            /* checks the arguments and make sure they are valid 
               returns MUM vs MEM designation based on input arguments*/
            if (input_list.length() && !is_file(input_list)) // provided a file-list
                FATAL_ERROR("The provided file-list is not valid.");
            else if (input_list.length() && is_file(input_list) && (files.size() > 0)) {
                FORCE_LOG("build_main", "Using filelist, ignoring positional args");
                files.clear();
            }
            else if (input_list.length() == 0 && (files.size() == 0) && (!from_parse_flag && !arrays_in_flag))
                FATAL_ERROR("Need to provide a file-list or files as positional args for processing.");
            
            for (auto f : files) {
                if (!is_file(f)) {
                    FATAL_ERROR(("The following file path is not valid: " + f).c_str());
                }
            }
            std::filesystem::path p (output_prefix);
            if (p.parent_path().string().empty() && !p.string().empty())
                output_prefix = "./" + output_prefix;
            else if (!is_dir(p.parent_path().string()))
                std::filesystem::create_directories(p.parent_path());

            if (from_parse_flag)
                parse_prefix = parse_prefix + ".fna";
            // if (max_mem_freq < -1)
            //     FATAL_ERROR("Maximum MEM frequency cannot be negative (-1 indicates no limit on MEM frequency)"); 

            if (rare_freq < 0)
                FATAL_ERROR("Per-sequence MEM frequency must be > 0 (or 0 for no limit)."); 

            if (binary && rare_freq != 1) {
                FORCE_LOG("build_main", "binary output is not supported for multi-MEMs, ignoring flag");
                binary = false;
            }

            if ((only_parse && use_gsacak) || (only_parse && arrays_in_flag) || (only_parse && from_parse_flag)) {
                only_parse = false;
                FORCE_LOG("build_main", "only-parse flag is not supported with use-gsacak, arrays-in, or from-parse, ignoring flag");
            }

            return (rare_freq == 1);
        }

        void set_parameters(size_t num_docs, bool mum_mode) {
            /* Set the main parameters:
            1) number of unique documents
            2) max frequence per doc
            3) and max total frequency */
            std::string match_type = mum_mode ? "MUMs" : "MEMs";
            // Set number of unique documents, based on valid ranges
            if (num_distinct_docs < -static_cast<int>(num_docs))
            {
                std::string message = "Too few number of sequences, defaulting to multi-" + match_type + " in 2 or more sequences";
                FORCE_LOG("build_main", message.c_str());
                num_distinct_docs = 2;
            }
            else if (num_distinct_docs <= 0) {num_distinct_docs = num_docs + num_distinct_docs;}
            else if (num_distinct_docs == 1)
            {
                std::string match_type = mum_mode ? "MUMs" : "MEMs";
                std::string message = "Too few number of sequences, defaulting to multi-" + match_type + " in 2 or more sequences";
                FORCE_LOG("build_main", message.c_str());
                num_distinct_docs = 2;
            }
            else if (num_distinct_docs >= num_docs)
            {
                std::string match_type = mum_mode ? "MUMs" : "MEMs";
                std::string message = "Too large number of sequences, defaulting to multi-" + match_type + " in all sequences";
                FORCE_LOG("build_main", message.c_str());
                num_distinct_docs = num_docs;
            }

            // Set max total frequency, based on valid ranges
            if (max_mem_freq < -static_cast<int>(num_docs) || max_mem_freq == 1)
            {
                FORCE_LOG("build_main", "Invalid maximum total MEM frequency, defaulting to no upper limit");
                max_mem_freq = 0;
            }
            else if (max_mem_freq < 0) {max_mem_freq = num_docs + max_mem_freq;}

            // max per doc frequency overrides total frequency
            if (rare_freq > 0 && (max_mem_freq == 0 || max_mem_freq > rare_freq * num_docs)) {
                max_mem_freq = rare_freq * num_docs;
            }

        }
};

struct HelperPrograms {
  /* Contains paths to run helper programs */
  std::string base_path = "";
  std::string parseNT_bin = "newscanNT.x";
//   std::string parse_fasta_bin = "newscan.x";
//   std::string parse_bin = "pscan.x";
  
public:
  void build_paths(std::string base) {
      /* Takes the base path, and combines it with names to generate executable paths */
      base_path.assign(base);
      parseNT_bin.assign(base_path + parseNT_bin);
  }

  void validate() const {
      /* Makes sure that each path for an executable is valid */
      bool invalid_path = !is_file(parseNT_bin);
      if (invalid_path) {FATAL_ERROR("One or more of helper program paths are invalid.");}
  }
};

/* Function Declartions involving structs */
void parse_build_options(int argc, char** argv, BuildOptions* opts);
void print_build_status_info(BuildOptions& opts, RefBuilder&ref_build, bool mum_mode);
void run_build_parse_cmd(BuildOptions* build_opts, HelperPrograms* helper_bins);

const std::string SKULL =
"                            ,--.   \n"
"                           {    }  \n"
"                           K,   }  \n"
"                          /  ~Y`   \n"
"                     ,   /   /     \n"
"                    {_'-K.__/      \n"
"                      `/-.__L._    \n"
"                      /  ' /`\\_}   \n"
"                     /  ' /        \n"
"             ____   /  ' /         \n"
"      ,-'~~~~    ~~/  ' /_         \n"
"    ,'             ``~~~  ',       \n"
"   (                        Y      \n"
"  {                         I      \n"
" {      -                    `,    \n"
" |       ',                   )    \n"
" |        |   ,..__      __. Y     \n"
" |    .,_./  Y ' / ^Y   J   )|     \n"
" \\           |' /   |   |   ||    \n"
"  \\          L_/    . _ (_,.'(     \n"
"   \\,   ,      ^^""' / |      )    \n"
"     \\_  \\          /,L]     /     \n"
"       '-_~-,       ` `   ./`      \n"
"          `'{_            )       \n"
"              ^^\\..___,.--`      \n";

#endif /* End of PFP_MUM_H */