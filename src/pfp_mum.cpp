/*
 * File: pfp_mum.cpp
 * Description: Structure adapted from docprofiles: https://github.com/oma219/docprofiles/tree/main
 *              Main file for building the prefix-free parse, and outputing MUMs along with r-index data structures (BWT/SA/LCP). 
 *              This workflow is based on doc profiles by Omar Ahmed, 
 *              which is in turn based on the pfp_lcp.hpp developed by Massimiliano Rossi.
 *
 * Date: December 20th, 2023
 */

#include <pfp_mum.hpp> 
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <filesystem>
#include <ref_builder.hpp>
#include <vector>
#include <pfp.hpp>
#include <pfp_lcp_mum.hpp>
#include <mem_finder.hpp>
// #include <mum_finder.hpp>
#include <getopt.h>
#include <queue>
#include <read_arrays.hpp>

int build_main(int argc, char** argv) {
    /* main method for finding matches */
    if (argc == 1) return mumemto_short_usage();

    // grab the command-line options, and validate them
    BuildOptions build_opts;
    parse_build_options(argc, argv, &build_opts);
    // print_build_status_info(&build_opts);
    bool mum_mode = build_opts.validate();

    // determine output path for reference, generate and store filelist
    build_opts.output_ref.assign(build_opts.output_prefix + ".fna");

    if ((build_opts.input_list.length() == 0) && !(build_opts.from_parse_flag || build_opts.arrays_in_flag)) {build_opts.input_list = make_filelist(build_opts.files, build_opts.output_prefix);}

    // Declare ref_build first
    RefBuilder ref_build = (build_opts.from_parse_flag || build_opts.arrays_in_flag)
        ? RefBuilder(build_opts.from_parse_flag ? build_opts.parse_prefix : build_opts.arrays_in, build_opts.use_rcomp)
        : RefBuilder(build_opts.input_list, build_opts.output_prefix, build_opts.use_rcomp);

    // normalize and reconcile the input parameters
    build_opts.set_parameters(ref_build.num_docs, mum_mode);

    // print parameters and begin pipeline
    print_build_status_info(build_opts, ref_build, mum_mode);

    // Build the input reference file, and bitvector labeling the end for each doc
    STATUS_LOG("build_main", "building the reference file based on file-list");
    auto start = std::chrono::system_clock::now();
    int input_file_status = ref_build.build_input_file();
    if (input_file_status == 1) {
        remove_temp_files(build_opts.output_prefix);
        FATAL_ERROR("Please check the input files and ensure that it contains valid FASTA files. Cleaning up...");
    }
    DONE_LOG((std::chrono::system_clock::now() - start));

    // Determine the paths to the BigBWT executables
    HelperPrograms helper_bins;
    // if (!std::getenv("PFPMUM_BUILD_DIR")) {FATAL_ERROR("Need to set PFPMUM_BUILD_DIR environment variable.");}
    std::filesystem::path path;
    if (!std::getenv("PFPMUM_BUILD_DIR")) {
        path = std::filesystem::canonical("/proc/self/exe").parent_path();
    }
    else {
        path = std::filesystem::path(std::string(std::getenv("PFPMUM_BUILD_DIR")));
    }
    helper_bins.build_paths((path / "").string());
    helper_bins.validate();

    // just read from files if provided
    if (build_opts.arrays_in_flag) {
        file_lcp input_lcp(build_opts.arrays_in, &ref_build);
        start = std::chrono::system_clock::now();
        STATUS_LOG("build_main", "finding multi-%ss from pfp", mum_mode ? "MUM" : "MEM");
        mem_finder match_finder(build_opts.output_prefix, ref_build, build_opts.min_match_len, build_opts.num_distinct_docs, build_opts.rare_freq, build_opts.max_mem_freq, build_opts.binary, build_opts.merge);
        input_lcp.process(match_finder);
        match_finder.close();
        input_lcp.close();

        auto sec = std::chrono::duration<double>((std::chrono::system_clock::now() - start)); std::fprintf(stderr, " done.  (%.3f sec)\n", sec.count());

        FORCE_LOG("build_main", "finished computing matches");
        return 0;
    }

    if (!build_opts.from_parse_flag){
        // Parse the input text with BigBWT, and load it into pf object
        STATUS_LOG("build_main", "generating the prefix-free parse for given reference");
        start = std::chrono::system_clock::now();

        run_build_parse_cmd(&build_opts, &helper_bins);
        DONE_LOG((std::chrono::system_clock::now() - start));
    }
    STATUS_LOG("build_main", "building the parse and dictionary objects");
    start = std::chrono::system_clock::now();
    
    pf_parsing pf = build_opts.from_parse_flag
        ? pf_parsing(build_opts.parse_prefix, build_opts.pfp_w)
        : pf_parsing(build_opts.output_ref, build_opts.pfp_w);

    DONE_LOG((std::chrono::system_clock::now() - start));

    std::cerr << "\n";

    // Builds the BWT, SA, LCP, and MUMs and writes to a file
    
    start = std::chrono::system_clock::now();
    
    pfp_lcp lcp(pf, build_opts.output_prefix, &ref_build, build_opts.arrays_out);
    size_t count = 0;

    STATUS_LOG("build_main", "finding multi-%ss from pfp", mum_mode ? "MUM" : "MEM");
    
    mem_finder match_finder(build_opts.output_prefix, ref_build, build_opts.min_match_len, build_opts.num_distinct_docs, build_opts.rare_freq, build_opts.max_mem_freq, build_opts.binary, build_opts.merge);
    count = lcp.process(match_finder);
    match_finder.close();
    lcp.close();
    
    auto sec = std::chrono::duration<double>((std::chrono::system_clock::now() - start)); std::fprintf(stderr, " done.  (%.3f sec)\n", sec.count());

    FORCE_LOG("build_main", "Found %d matches!", count);
    FORCE_LOG("build_main", "Run stats:");
    double n_r = static_cast<double>(pf.n) / lcp.num_run;
    std::cerr << "n/r = " << pf.n << " / " << lcp.num_run << " = " << std::fixed << std::setprecision(3) << std::round(n_r * 1000) / 1000 << std::endl;

    

    if (!build_opts.keep_temp && !build_opts.from_parse_flag)
        remove_temp_files(build_opts.output_prefix);
    std::cerr << "\n";
    
    return 0;
}

void run_build_parse_cmd(BuildOptions* build_opts, HelperPrograms* helper_bins) {
    // Generates and runs the command-line for executing the PFP of the reference 
    std::ostringstream command_stream;
    // if (build_opts->threads > 0) {
    //     std::string curr_exe = "";
    //     if (build_opts->is_fasta) {curr_exe.assign(helper_bins->parse_fasta_bin);}
    //     else {curr_exe.assign(helper_bins->parse_bin);}

    //     command_stream << curr_exe << " "; //<< " -i ";
    //     command_stream << build_opts->output_ref << " ";
    //     command_stream << "-w " << build_opts->pfp_w;
    //     command_stream << " -p " << build_opts->hash_mod;
    //     // command_stream << " -t " << build_opts->threads;
    // }
    // else {
        std::string curr_exe = "";
        command_stream << helper_bins->parseNT_bin << " "; // << " -i ";
        command_stream << build_opts->output_ref << " ";
        command_stream << "-w " << build_opts->pfp_w;
        command_stream << " -p " << build_opts->hash_mod;
    // }
    if (build_opts->is_fasta) {command_stream << " -f";}

    // std::cout << command_stream.str() << std::endl;
    // std::cout << "Executing this command: " << command_stream.str().c_str() << std::endl;
    auto parse_log = execute_cmd(command_stream.str().c_str());
    //OTHER_LOG(parse_log.data());
}

std::string execute_cmd(const char* cmd) {
    std::array<char, 256> buffer{};
    std::string output = "";

    std::string cmd_plus_stderr = std::string(cmd) + " 2>&1";
    FILE* pipe = popen(cmd_plus_stderr.data(), "r"); // Extract stderr as well
    if (!pipe) {FATAL_ERROR("popen() failed!");}

    try {
        std::size_t bytes;
        while ((bytes = fread(buffer.data(), sizeof(char), sizeof(buffer), pipe))) {
            output += std::string(buffer.data(), bytes);
        }
    } catch (...) {
        pclose(pipe);
        FATAL_ERROR("Error occurred while reading popen() stream.");
    }

    // Correct the usage of pclose and WEXITSTATUS
    int status = pclose(pipe);
    size_t exit_code = WEXITSTATUS(status);

    // Check if the command failed
    if (exit_code != 0) {
        std::cout << "\n";
        FATAL_ERROR("external command failed ... here is the error message:\n%s", output.data());
    }

    return output;
}

bool endsWith(const std::string& str, const std::string& suffix) {
    // Checks if the string ends the suffix
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

bool is_integer(const std::string& str) {
    /* Checks if string passed is an integer */
    std::string::const_iterator iter = str.begin();
    while (iter != str.end() && std::isdigit(*iter)) {++iter;}
    return !str.empty() && iter == str.end();
}

std::vector<std::string> split(std::string input, char delim) {
    /* Takes in a string, and splits it based on delimiters */
    std::vector<std::string> word_list;
    std::string curr_word = "";

    for (char ch: input) {
        if (ch == delim && curr_word.length()) {word_list.push_back(curr_word); curr_word = "";}
        else if (ch != delim) {curr_word += ch;}
    }
    if (curr_word.length()) {word_list.push_back(curr_word);}
    return word_list;
}

int is_file(std::string path) {
    /* Checks if the path is a valid file-path */
    std::ifstream test_file(path.data());
    if (test_file.fail()) {return 0;}
    
    test_file.close();
    return 1;
}

int is_dir(std::string path) {
    /* Checks if the directory is a valid path */
    return std::filesystem::exists(path);
}

void print_build_status_info(BuildOptions& opts, RefBuilder& ref_build, bool mum_mode) {
    /* prints out the information being used in the current run */
    std::fprintf(stderr, "\nOverview of Parameters:\n");
    if (opts.from_parse_flag) {
        std::fprintf(stderr, "\tUsing pre-computed PFP files with prefix: %s\n", opts.parse_prefix.data());
    }
    else if (opts.arrays_in.length() > 0)
        std::fprintf(stderr, "\tUsing pre-computed LCP/BWT/SA arrays from files with prefix: %s\n", opts.arrays_in.data());
    else if (opts.files.size() > 5) {
        std::fprintf(stderr, "\tInput files (N = %d): ", ref_build.num_docs);
        std::fprintf(stderr, "%s,", opts.files.at(0).data());
        std::fprintf(stderr, "%s,", opts.files.at(1).data());
        std::fprintf(stderr, " ... ,", opts.files.at(1).data());
        std::fprintf(stderr, "%s,", opts.files.at(opts.files.size() - 2).data());
        std::fprintf(stderr, "%s\n", opts.files.at(opts.files.size() - 1).data());
    }
    else if (opts.files.size() > 0) {
        std::fprintf(stderr, "\tInput files (N = %d): ", ref_build.num_docs);
        for (int i = 0; i < opts.files.size() - 1; i++)
            {
                std::fprintf(stderr, "%s,", opts.files.at(i).data());
            }
            std::fprintf(stderr, "%s\n", opts.files.at(opts.files.size() - 1).data());
    }
    else
        std::fprintf(stderr, "\tInput file-list (N = %d): %s\n", ref_build.num_docs, opts.input_list.data());

    std::string match_type = mum_mode ? "MUM" : "MEM";
    std::fprintf(stderr, "\tOutput path: %s\n", opts.output_prefix.data());
    if (!opts.from_parse_flag && !opts.arrays_in_flag)
        std::fprintf(stderr, "\tPFP window size: %d\n", opts.pfp_w);
    if (opts.arrays_out) {std::fprintf(stderr, "\tWriting LCP, BWT and suffix arrays\n");}
    std::fprintf(stderr, "\tMinimum %s length: %d\n", match_type.data(), opts.min_match_len);
    std::fprintf(stderr, "\tInclude reverse complement?: %s\n", opts.use_rcomp ? "True" : "False");
    std::fprintf(stderr, "\tMax occurences per sequence: %d\n", opts.rare_freq);
    if (opts.max_mem_freq < opts.rare_freq * ref_build.num_docs)
        std::fprintf(stderr, "\t\t- excluding multi-MEMs that occur more than %d times in total\n", opts.max_mem_freq);
    if (opts.num_distinct_docs == ref_build.num_docs)
        std::fprintf(stderr, "\tfinding multi-%ss present in all genomes\n", match_type.data());
    else
        std::fprintf(stderr, "\tfinding multi-%ss present in %d genomes\n", match_type.data(), opts.num_distinct_docs);
    if (opts.merge) {
        std::fprintf(stderr, "\twriting mergeable output\n");
    }
    std::fprintf(stderr, "\n");
}

std::string make_filelist(std::vector<std::string> files, std::string output_prefix) {
    std::string fname = output_prefix + "_filelist.txt";
    std::ofstream outfile(fname);
    for (size_t i = 0; i < files.size(); ++i) {
        outfile << std::filesystem::absolute(files[i]).string() << std::endl;
    }
    outfile.close();
    return fname;
}

void remove_temp_files(std::string filename) {
    std::vector<std::string> temp_files = {".fna", ".fna.dict", ".fna.occ", ".fna.parse_old", ".fna.last", ".fna.parse", "_filelist.txt"};
    for (auto &ext : temp_files) {
        std::filesystem::remove(std::filesystem::path(filename + ext));
    }
}


void parse_build_options(int argc, char** argv, BuildOptions* opts) {
    /* parses the arguments for the build sub-command, and returns a struct with arguments */

    static struct option long_options[] = {
        {"help",      no_argument, NULL,  'h'},
        {"input",   required_argument, NULL,  'i'},
        {"output",       required_argument, NULL,  'o'},
        {"revcomp",   no_argument, NULL,  'r'},
        {"minimum-genomes",   required_argument, NULL,  'k'},
        {"no-overlap",   no_argument, NULL,  's'},
        {"modulus", required_argument, NULL, 'm'},
        {"from-parse",   no_argument, NULL,  'p'},
        {"min-match-len",   required_argument, NULL,  'l'},
        {"max-freq",   required_argument, NULL,  'F'},
        {"arrays-out",   no_argument, NULL,  'A'},
        {"arrays-in",   required_argument, NULL,  'a'},
        {"keep-temp-files",   no_argument, NULL,  'K'},
        {"window",   required_argument, NULL,  'w'},
        {"rare",   required_argument, NULL,  'f'},
        {"binary",   no_argument, NULL,  'b'},
        {"merge",   no_argument, NULL,  'M'},
        {0, 0, 0,  0}
    };
    int c = 0;
    int long_index = 0;
    
    while ((c = getopt_long(argc, argv, "hi:F:o:w:sl:ra:AKk:p:m:f:bM", long_options, &long_index)) >= 0) {
        switch(c) {
            case 'h': mumemto_usage(); std::exit(0);
            case 'i': opts->input_list.assign(optarg); break;
            case 'o': opts->output_prefix.assign(optarg); break;
            case 'w': opts->pfp_w = std::atoi(optarg); break;
            case 'r': opts->use_rcomp = false; break;
            case 's': opts->overlap = false; break;
            case 'k': opts->num_distinct_docs = std::atoi(optarg); break;
            case 'm': opts->hash_mod = std::atoi(optarg); break;
            case 'p': opts->parse_prefix.assign(optarg); opts->from_parse_flag = true; break;
            case 'l': opts->min_match_len = std::atoi(optarg); break;
            case 'F': opts->max_mem_freq = std::atoi(optarg); break;
            case 'A': opts->arrays_out = true; break;
            case 'a': opts->arrays_in.assign(optarg); opts->arrays_in_flag = true; break;
            case 'K': opts->keep_temp = true; break;
            case 'f': opts->rare_freq = std::atoi(optarg); break;
            case 'b': opts->binary = true; break;
            case 'M': opts->merge = true; break;
            default: mumemto_usage(); std::exit(1);
        }
    }

    for (int i = optind; i < argc; ++i) {
        opts->files.push_back(argv[i]);
    }
}

int mumemto_usage() {
    /* prints out the usage information for the build method */
    std::fprintf(stderr, "\nmumemto - find maximal [unique | exact] matches using PFP.\n");
    std::fprintf(stderr, "Usage: mumemto [mum | mem] [options] [input_fasta [...]]\n\n");
    std::fprintf(stderr, "*** for all options, N = # of sequences ***\n");
    std::fprintf(stderr, "I/O options:\n");
    std::fprintf(stderr, "\t%-32sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-22s%-10spath to a file-list of genomes to use (overrides positional args)\n", "-i, --input", "[FILE]");
    std::fprintf(stderr, "\t%-22s%-10soutput prefix path\n", "-o, --output", "[arg]");
    std::fprintf(stderr, "\t%-32sinclude the reverse complement of the sequences (default: true)\n\n", "-r, --no-revcomp");
    std::fprintf(stderr, "\t%-32soutput binary format (multi-MUMs only)\n\n", "-b, --binary");

    std::fprintf(stderr, "\t%-32swrite LCP, BWT, and SA to file\n", "-A, --arrays-out");
    std::fprintf(stderr, "\t%-22s%-10scompute matches from precomputed LCP, BWT, SA (with shared PREFIX.bwt/sa/lcp)\n\n", "-a, --arrays-in", "[PREFIX]");
    std::fprintf(stderr, "\t%-32soutput extra metadata to enable merging multi-MUMs\n", "-M, --merge");

    
    std::fprintf(stderr, "Exact match parameters:\n");
    std::fprintf(stderr, "\t%-22s%-10sminimum MUM or MEM length (default: 20)\n\n", "-l, --min-match-len", "[INT]");
    std::fprintf(stderr, "\t%-22s%-10sfind matches in at least k sequences (k < 0 sets the sequences relative to N, i.e. matches must occur in at least N - |k| sequences)\n\t%-32s(default: match occurs in all sequences, i.e. strict multi-MUM/MEM)\n\n", "-k, --minimum-genomes", "[INT]", "");
    std::fprintf(stderr, "\t%-22s%-10smaximum number of occurences per sequence (set to 0 for no upper limit)\n\t%-32s(default: 1 [unique])\n\n", "-f, --per-seq-freq", "[INT]", "");
    std::fprintf(stderr, "\t%-22s%-10smaximum number of total occurences of match (if negative, then relative to N) \n\t%-32s(default: no upper limit)\n\n", "-F, --max-total-freq","[INT]", "");
    

    std::fprintf(stderr, "PFP options:\n");
    std::fprintf(stderr, "\t%-22s%-10swindow size used for pfp (default: 10)\n", "-w, --window", "[INT]");
    std::fprintf(stderr, "\t%-22s%-10shash-modulus used for pfp (default: 100)\n", "-m, --modulus", "[INT]");
    std::fprintf(stderr, "\t%-32suse pre-computed pf-parse\n", "-p, --from-parse");
    std::fprintf(stderr, "\t%-32skeep PFP files\n\n", "-K, --keep-temp-files");

    std::fprintf(stderr, "Overview:\n");
        std::fprintf(stderr, "\tBy default, Mumemto computes multi-MUMs. Exact match parameters can be additionally tuned in three main ways:\n");
        std::fprintf(stderr, "\t1) Choosing the number of sequences a match must appear in [-k]\n");
        std::fprintf(stderr, "\t2) Choosing the maximum number of occurences in each genome [-f]\n");
        std::fprintf(stderr, "\t3) Choosing the total maximum number of occurences [-F]\n");
    std::fprintf(stderr, "Examples:\n");
    std::fprintf(stderr, "\t - Find all strict multi-MUMs across a collection: mumemto [OPTIONS] [input_fasta [...]] (equivalently -k 0 -f 1 -F 0)\n");
    std::fprintf(stderr, "\t - Find partial multi-MUMs in all sequences but one: mumemto -k -1 [OPTIONS] [input_fasta [...]]\n");
    std::fprintf(stderr, "\t - Find multi-MEMs that appear at most 3 times in each sequence: mumemto -f 3 [OPTIONS] [input_fasta [...]]\n");
    std::fprintf(stderr, "\t - Find all MEMs that appear at most 100 times within a collection: mumemto -f 0 -k 2 -F 100 [OPTIONS] [input_fasta [...]]\n");
    return 0;
}

int mumemto_short_usage() {
    /* Prints the usage information for mumemto */
    std::fprintf(stderr, "\nmumemto - find maximal [unique | exact] matches using PFP.\n");
    std::fprintf(stderr, "Usage: mumemto [options] [input_fasta [...]]\n");
    std::fprintf(stderr, "\t%-18sprints detailed usage message\n", "-h, --help");
    return 0;
}

int main(int argc, char** argv) {
    /* main method for mumemto */
    std::fprintf(stderr, "\033[1m\033[31m\nmumemto version: %s\033[m\033[0m\n", PFPMUM_VERSION);
    if (argc > 1) {
        if (std::strcmp(argv[1], "mori") == 0)
        {
            std::fprintf(stdout, "\033[1m\033[31m\nDeath is inevitable.\033[m\033[0m\n");
            std::fprintf(stdout, SKULL.c_str());
            return 0;
        }
        return build_main(argc, argv);
    }
    return mumemto_short_usage();
}
