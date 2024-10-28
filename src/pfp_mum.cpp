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
#include <mum_finder.hpp>
#include <getopt.h>
#include <queue>
#include <read_arrays.hpp>

int build_main(int argc, char** argv, bool mum_mode) {
    /* main method for finding matches */
    if (argc == 1) return mumemto_build_usage();

    // grab the command-line options, and validate them
    BuildOptions build_opts;
    parse_build_options(argc, argv, &build_opts);
    // print_build_status_info(&build_opts);
    build_opts.validate();

    // determine output path for reference, and print all info
    build_opts.output_ref.assign(build_opts.output_prefix + ".fa");
    print_build_status_info(&build_opts, mum_mode);

    if (build_opts.input_list.length() == 0) {build_opts.input_list = make_filelist(build_opts.files, build_opts.output_prefix);}
    // Build the input reference file, and bitvector labeling the end for each doc
    STATUS_LOG("build_main", "building the reference file based on file-list");
    auto start = std::chrono::system_clock::now();
    
    RefBuilder ref_build(build_opts.input_list, build_opts.output_prefix, build_opts.use_rcomp);
    DONE_LOG((std::chrono::system_clock::now() - start));

    // Make sure that document numbers can be store in 2 bytes
    if (ref_build.num_docs >= MAXDOCS)
        FATAL_ERROR("An index cannot be build over %ld documents, "
                    "please reduce to a max of 65,535 docs.", ref_build.num_docs);

    if (1) {
        bool use_input_paths = build_opts.files.size() == ref_build.seq_lengths.size();
        int includes_rc = build_opts.use_rcomp ? 2 : 1;
        std::string lengths_fname = build_opts.output_prefix + ".lengths";
        std::ofstream outfile(lengths_fname);
        std::string doc_name;
        for (size_t i = 0; i < ref_build.seq_lengths.size() - 1; ++i) {
            doc_name = use_input_paths ? build_opts.files[i] : ("sequence_" + std::to_string(i + 1));
            outfile << doc_name << " " << ref_build.seq_lengths[i] / includes_rc << std::endl;
        }
        doc_name = use_input_paths ? build_opts.files[ref_build.seq_lengths.size() - 1] : ("sequence_" + std::to_string(ref_build.seq_lengths.size()));
        outfile << doc_name << " " << (ref_build.seq_lengths[ref_build.seq_lengths.size() - 1] - 1) / includes_rc << std::endl;
        outfile.close();
    }

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
    helper_bins.build_paths((path / "bin/").string());
    helper_bins.validate();

    if (build_opts.missing_genomes >= ref_build.num_docs - 1)
    {
        std::string match_type = mum_mode ? "MUMs" : "MEMs";
        std::string message = "Too few number of sequences, defaulting to multi-" + match_type + " in 2 or more sequences";
        FORCE_LOG("build_main", message.c_str());
        build_opts.missing_genomes = ref_build.num_docs - 2;
    }

    // just read from files if provided
    if (build_opts.arrays_in.length() > 0) {
        file_lcp input_lcp(build_opts.arrays_in, &ref_build);
        start = std::chrono::system_clock::now();
        if (mum_mode){
            STATUS_LOG("build_main", "finding multi-MUMs from pfp\n");
            mum_finder match_finder(build_opts.output_prefix, &ref_build, build_opts.min_match_len, build_opts.missing_genomes + 1, build_opts.overlap, build_opts.col_mum_mode);
            input_lcp.process(match_finder);
            match_finder.close();
        }
        else {
            STATUS_LOG("build_main", "finding multi-MEMs from pfp");
            mem_finder match_finder(build_opts.output_prefix, &ref_build, build_opts.min_match_len, ref_build.num_docs - build_opts.missing_genomes, build_opts.max_mem_freq);
            input_lcp.process(match_finder);
            match_finder.close();
        }

        input_lcp.close();

        if (mum_mode)
            FORCE_LOG("build_main", "finding multi-MUMs from pfp");
        else
            FORCE_LOG("build_main", "finding multi-MEMs from pfp");
        DONE_LOG((std::chrono::system_clock::now() - start));

        FORCE_LOG("build_main", "finished computing matches");
        return 0;
    }

    if (!build_opts.from_parse){
        // Parse the input text with BigBWT, and load it into pf object
        STATUS_LOG("build_main", "generating the prefix-free parse for given reference");
        start = std::chrono::system_clock::now();

        run_build_parse_cmd(&build_opts, &helper_bins);
        DONE_LOG((std::chrono::system_clock::now() - start));
    }
    STATUS_LOG("build_main", "building the parse and dictionary objects");
    start = std::chrono::system_clock::now();

    pf_parsing pf(build_opts.output_ref, build_opts.pfp_w);
    DONE_LOG((std::chrono::system_clock::now() - start));

    std::cerr << "\n";

    // Builds the BWT, SA, LCP, and MUMs and writes to a file
    
    start = std::chrono::system_clock::now();
    
    pfp_lcp lcp(pf, build_opts.output_prefix, &ref_build, build_opts.arrays_out, build_opts.rlbwt_out, build_opts.thresholds_out);
    size_t count = 0;

    if (build_opts.rare_freq > 1){
        STATUS_LOG("build_main", "finding rare multi-MEMs from pfp");
        rare_mem_finder match_finder(build_opts.output_prefix, &ref_build, build_opts.min_match_len, ref_build.num_docs - build_opts.missing_genomes, build_opts.rare_freq);
        count = lcp.process(match_finder);
        match_finder.close();
    }
    else if (mum_mode){
        STATUS_LOG("build_main", "finding multi-MUMs from pfp");
        mum_finder match_finder(build_opts.output_prefix, &ref_build, build_opts.min_match_len, build_opts.missing_genomes + 1, build_opts.overlap, build_opts.col_mum_mode);
        count = lcp.process(match_finder);
        match_finder.close();
    }
    else {
        STATUS_LOG("build_main", "finding multi-MEMs from pfp");
        mem_finder match_finder(build_opts.output_prefix, &ref_build, build_opts.min_match_len, ref_build.num_docs - build_opts.missing_genomes, build_opts.max_mem_freq);
        count = lcp.process(match_finder);
        match_finder.close();
    }

    lcp.close();
    
    auto sec = std::chrono::duration<double>((std::chrono::system_clock::now() - start)); std::fprintf(stderr, " done.  (%.3f sec)\n", sec.count());

    FORCE_LOG("build_main", "Found %d matches!", count);
    FORCE_LOG("build_main", "Run stats:");
    double n_r = static_cast<double>(pf.n) / lcp.num_run;
    std::cerr << "n/r = " << pf.n << " / " << lcp.num_run << " = " << std::fixed << std::setprecision(3) << std::round(n_r * 1000) / 1000 << std::endl;

    

    if (!build_opts.keep_temp && !build_opts.from_parse)
        remove_temp_files(build_opts.output_ref);
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

    // Check if the command failed
    size_t exit_code = WEXITSTATUS(pclose(pipe));
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

void print_build_status_info(BuildOptions* opts, bool mum_mode) {
    /* prints out the information being used in the current run */
    std::fprintf(stderr, "\nOverview of Parameters:\n");
    if (opts->input_list.length())
        std::fprintf(stderr, "\tInput file-list: %s\n", opts->input_list.data());
    else if (opts->files.size() > 5) {
        std::fprintf(stderr, "\tInput files: ");
        std::fprintf(stderr, "%s,", opts->files.at(0).data());
        std::fprintf(stderr, "%s,", opts->files.at(1).data());
        std::fprintf(stderr, " ... ,", opts->files.at(1).data());
        std::fprintf(stderr, "%s,", opts->files.at(opts->files.size() - 2).data());
        std::fprintf(stderr, "%s\n", opts->files.at(opts->files.size() - 1).data());
    }
    else {
        std::fprintf(stderr, "\tInput files: ");
        for (int i = 0; i < opts->files.size() - 1; i++)
            {
                std::fprintf(stderr, "%s,", opts->files.at(i).data());
            }
            std::fprintf(stderr, "%s\n", opts->files.at(opts->files.size() - 1).data());
    }
    std::string match_type = mum_mode ? "MUM" : "MEM";
    std::fprintf(stderr, "\tOutput ref path: %s\n", opts->output_ref.data());
    if (opts->arrays_in.length() > 0)
        std::fprintf(stderr, "\tUsing pre-computed LCP/BWT/SA arrays from files with prefix: %s\n", opts->output_prefix.data());
    else 
        std::fprintf(stderr, "\tPFP window size: %d\n", opts->pfp_w);
    if (opts->arrays_out) {std::fprintf(stderr, "\tWriting LCP, BWT and suffix arrays\n");}
    if (opts->rlbwt_out) {std::fprintf(stderr, "\tWriting RLBWT and run-sampled suffix arrays\n");}
    if (opts->thresholds_out) {std::fprintf(stderr, "\tWriting thresholds\n");}
    if (opts->col_mum_mode) {std::fprintf(stderr, "\tFinding col-MUMs, written in byte output\n");}
    else {std::fprintf(stderr, "\tFinding multi-MUMs written in human readable output\n");}
    std::fprintf(stderr, "\tMinimum %s length: %d\n", match_type.data(), opts->min_match_len);
    std::fprintf(stderr, "\tInclude reverse complement?: %d\n", opts->use_rcomp);
    if (opts->missing_genomes == 0)
        std::fprintf(stderr, "\tfinding multi-%ss present in all genomes\n", match_type.data(), opts->missing_genomes);
    else
        std::fprintf(stderr, "\tfinding multi-%ss present in N - %d genomes\n", match_type.data(), opts->missing_genomes);
    if (!mum_mode && opts->max_mem_freq >= 0)
        std::fprintf(stderr, "\t\t- excluding multi-MEMs that occur more than N + %d times\n", opts->max_mem_freq);
    if (!mum_mode && opts->rare_freq > 1)
        std::fprintf(stderr, "\t\t- only finding rare multi-MEMs, that occur at most %d times in each sequence\n", opts->rare_freq);
    if (opts->overlap && opts->missing_genomes > 0 && mum_mode)
        std::fprintf(stderr, "\t\t- including overlapping multi-MUMs\n");
}

std::string make_filelist(std::vector<std::string> files, std::string output_prefix) {
    std::string fname = output_prefix + "_filelist.txt";
    std::ofstream outfile(fname);
    for (size_t i = 0; i < files.size(); ++i) {
        outfile << files[i] << " " << i + 1 << std::endl;
    }
    outfile.close();
    return fname;
}

void remove_temp_files(std::string filename) {
    std::vector<std::string> temp_files = {".dict", ".occ", ".parse_old", ".last", ".parse"};
    for (auto &ext : temp_files) {
        std::filesystem::remove(std::filesystem::path(filename + ext));
    }
    std::filesystem::remove(std::filesystem::path(filename));
}


void parse_build_options(int argc, char** argv, BuildOptions* opts) {
    /* parses the arguments for the build sub-command, and returns a struct with arguments */

    static struct option long_options[] = {
        {"help",      no_argument, NULL,  'h'},
        {"input",   required_argument, NULL,  'i'},
        {"output",       required_argument, NULL,  'o'},
        {"revcomp",   no_argument, NULL,  'r'},
        {"missing-genomes",   required_argument, NULL,  'k'},
        {"no-overlap",   no_argument, NULL,  'p'},
        {"modulus", required_argument, NULL, 'm'},
        {"from-parse",   no_argument, NULL,  's'},
        {"min-match-len",   required_argument, NULL,  'l'},
        {"max-freq",   required_argument, NULL,  'F'},
        {"arrays-out",   no_argument, NULL,  'A'},
        {"arrays-in",   required_argument, NULL,  'a'},
        {"rlbwt-out",   no_argument, NULL,  'R'},
        {"thresholds-out",   no_argument, NULL,  't'},
        {"regular-mum",   no_argument, NULL,  'M'},
        {"keep-temp-files",   no_argument, NULL,  'K'},
        {"window",   required_argument, NULL,  'w'},
        {"rare",   required_argument, NULL,  'x'},
        {0, 0, 0,  0}
    };
    int c = 0;
    int long_index = 0;
    
    while ((c = getopt_long(argc, argv, "hi:F:o:w:sl:ra:AKk:p:m:x:RTM", long_options, &long_index)) >= 0) {
        switch(c) {
            case 'h': mumemto_build_usage(); std::exit(1);
            case 'i': opts->input_list.assign(optarg); break;
            case 'o': opts->output_prefix.assign(optarg); break;
            case 'w': opts->pfp_w = std::atoi(optarg); break;
            case 'r': opts->use_rcomp = true; break;
            case 'p': opts->overlap = false; break;
            case 'k': opts->missing_genomes = std::atoi(optarg); break;
            case 'm': opts->hash_mod = std::atoi(optarg); break;
            case 's': opts->from_parse = true; break;
            case 'l': opts->min_match_len = std::atoi(optarg); break;
            case 'F': opts->max_mem_freq = std::atoi(optarg); break;
            case 'R': opts->rlbwt_out = true; break;
            case 'T': opts->thresholds_out = true; break;
            case 'M': opts->col_mum_mode = false; break;
            case 'A': opts->arrays_out = true; break;
            case 'a': opts->arrays_in.assign(optarg); break;
            case 'K': opts->keep_temp = true; break;
            case 'x': opts->rare_freq = std::atoi(optarg); break;
            default: mumemto_build_usage(); std::exit(1);
        }
    }

    for (int i = optind; i < argc; ++i) {
        opts->files.push_back(argv[i]);
    }
}

int mumemto_build_usage() {
    /* prints out the usage information for the build method */
    std::fprintf(stderr, "\nmumemto - find maximal [unique | exact] matches using PFP.\n");
    std::fprintf(stderr, "Usage: mumemto [mum | mem] [options] [input_fasta [...]]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-28sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-18s%-10spath to a file-list of genomes to use (overrides positional args)\n", "-i, --input", "[FILE]");
    std::fprintf(stderr, "\t%-18s%-10soutput prefix path\n", "-o, --output", "[arg]");
    std::fprintf(stderr, "\t%-28sinclude the reverse-complement of sequence (default: false)\n\n", "-r, --revcomp");
    std::fprintf(stderr, "\t%-18s%-10sminimum MUM or MEM length (default: 20)\n\n", "-l, --min-match-len", "[INT]");
    std::fprintf(stderr, "\t%-28swrite LCP, BWT, and SA to file\n\n", "-A, --arrays-out");
    std::fprintf(stderr, "\t%-28swrite RLBWT and SA (run-sampled) to file\n\n", "-R, --rlbwt-out");
    std::fprintf(stderr, "\t%-28swrite threshold positions and values to file\n\n", "-T, --thresholds-out");
    std::fprintf(stderr, "\t%-18s%-10scompute matches from precomputed LCP, BWT, SA (with shared PREFIX.bwt/sa/lcp)\n\n", "-a, --arrays-in", "[PREFIX]");

    std::fprintf(stderr, "\t%-18s%-10sfind multi-MUMs or multi-MEMs in at least N - k genomes\n\t%-28s(default: 0, match must occur in all sequences, i.e. strict multi-MUM/MEM)\n\n", "-k, --missing-genomes", "[INT]", "");
    std::fprintf(stderr, "MUM mode options:\n");
    std::fprintf(stderr, "\t%-28soutput subset multi-MUMs that overlap shorter, more complete multi-MUMs (default: true w/ -k)\n", "-p, --no-overlap");
    std::fprintf(stderr, "\t%-28scompute regular multi-MUM output (default: col-MUMs in binary output)\n\n", "-M, --regular-mum");
    
    std::fprintf(stderr, "MEM mode options:\n");
    std::fprintf(stderr, "\t%-18s%-10smaximum number of occurences of MEM, beyond N \n\t%-28s(default: -1, no upper limit. Note: -k 0 -F 0 will result in strict multi-MUMs)\n", "-F, --max-freq","[INT]", "");
    std::fprintf(stderr, "\t%-18s%-10sfind rare multi-MEMs, occuring at most k times in each genome\n", "--rare", "[INT]");


    std::fprintf(stderr, "PFP options:\n");
    std::fprintf(stderr, "\t%-18s%-10swindow size used for pfp (default: 10)\n", "-w, --window", "[INT]");
    std::fprintf(stderr, "\t%-18s%-10shash-modulus used for pfp (default: 100)\n\n", "-m, --modulus", "[INT]");
    std::fprintf(stderr, "\t%-28suse pre-computed pf-parse\n\n", "-s, --from-parse");
    std::fprintf(stderr, "\t%-28skeep PFP files\n\n", "-K, --keep-temp-files");

    return 0;
}

int mumemto_usage() {
    /* Prints the usage information for mumemto */
    std::fprintf(stderr, "\nmumemto - find maximal [unique | exact] matches using PFP.\n");
    std::fprintf(stderr, "Usage: mumemto [mum | mem] [options] [input_fasta [...]]\n\n");

    std::fprintf(stderr, "\nmumemto has different modes to run:\n");
    std::fprintf(stderr, "\tmum\tcomputes maximal unique matches (MUMs) in the collection of sequences\n");
    std::fprintf(stderr, "\tmem\tcomputes maximal exact matches (MEMs) in the collection of sequences\n");
    return 0;
}

int main(int argc, char** argv) {
    /* main method for mumemto */
    std::fprintf(stderr, "\033[1m\033[31m\nmumemto version: %s\033[m\033[0m\n", PFPMUM_VERSION);
    bool mum_mode;
    if (argc > 1) {
        if (std::strcmp(argv[1], "mum") == 0)
            mum_mode = true;
        else if (std::strcmp(argv[1], "mem") == 0)
            mum_mode = false;
        else if (std::strcmp(argv[1], "mori") == 0)
        {
            std::fprintf(stdout, "\033[1m\033[31m\nDeath is inevitable.\033[m\033[0m\n");
            std::fprintf(stdout, SKULL.c_str());
            return 0;
        }
        else
        {
            std::fprintf(stderr, "\nOne of [mum | mem] mode selection required!\n");
            return mumemto_usage();
        }
        return build_main(argc-1, argv+1, mum_mode);
    }
    return mumemto_usage();
}
