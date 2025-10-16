/*
 * File: extract_mums.cpp
 * Description: Standalone executable for extracting MUMs from a MUM and length file.
 */

#include <iostream>
#include <string>
#include <getopt.h>
#include <filesystem>
#include <zlib.h>
#include <kseq.h>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>  // Include this for std::runtime_error
#include "parse_mums.hpp"


KSEQ_INIT(gzFile, gzread);

void print_usage() {
    std::fprintf(stderr, "\nextract_mums - extract MUMs from a MUM and length file\n");
    std::fprintf(stderr, "Usage: extract_mums [options] -m mum_file -o output_file\n\n");
    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-32sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-22s%-10spath to a mum file\n", "-m, --mums", "[FILE]");
    std::fprintf(stderr, "\t%-22s%-10spath to a length file (optional, uses associated length file if not provided)\n", "-l, --lengths", "[FILE]");
    std::fprintf(stderr, "\t%-22s%-10soutput path\n", "-o, --output", "[FILE]");
    std::fprintf(stderr, "\t%-32sdo not add terminator (#) to end of each MUM sequence\n", "-t, --no-terminator");
}

bool endsWith(const std::string& str, const std::string& suffix) {
    // Checks if the string ends the suffix
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

int is_file(std::string path) {
    /* Checks if the path is a valid file-path */
    return std::filesystem::is_regular_file(path);
}

int read_ref_file(const std::string& length_file_path, std::string& ref_seq) {
    std::ifstream length_file(length_file_path);
    std::string line;
    if (!std::getline(length_file, line)) {
        std::cerr << "Error: Length file is empty" << std::endl;
        return 1; // Error: Length file is empty
    }
    std::string fasta_path = line.substr(0, line.find_first_of(" \t\n\r"));
    if (!is_file(fasta_path)) {
        std::cerr << "Error: Invalid FASTA path in length file" << std::endl;
        return 1; // Error: Invalid FASTA path in length file
    }        

    // read in the fasta file using gzopen for both compressed and uncompressed files
    gzFile gzfp = gzopen(fasta_path.data(), "r");
    if (gzfp == 0) {
        std::cerr << "Error: Unable to open FASTA file" << std::endl;
        return 1; // Error opening file
    }

    kseq_t* seq = kseq_init(gzfp);
    ref_seq.clear();
    while (kseq_read(seq) >= 0) {
        ref_seq += seq->seq.s;
    }
    kseq_destroy(seq);
    gzclose(gzfp);

    if (ref_seq.empty()) {
        std::cerr << "Error: Empty input file found" << std::endl;
        return 1; // Error: Empty input file found
    }
    return 0; // Success
}

std::vector<std::pair<size_t, size_t>> parse_mums(const std::string& path) {
    std::vector<std::pair<size_t, size_t>> mum_list;
    if (endsWith(path, ".bumbl")) {
        mumsio::stream_bumbl_first(path, [&](uint32_t len, int64_t off0){
            mum_list.emplace_back(static_cast<size_t>(off0), static_cast<size_t>(len));
        });
    } else {
        mumsio::stream_mums_first(path, [&](uint32_t len, int64_t off0){
            mum_list.emplace_back(static_cast<size_t>(off0), static_cast<size_t>(len));
        });
    }
    return mum_list;
}

int extract_mums(std::string mum_file_path, std::string length_file_path, std::string output_path, bool add_terminator) {
    std::string ref_seq;
    int result = read_ref_file(length_file_path, ref_seq);
    if (result != 0) {
        return result;
    }
    std::vector<std::pair<size_t, size_t>> mum_list = parse_mums(mum_file_path);
    std::ofstream output_file(output_path);

    // extract the mums from the reference sequence
    std::string mum_seq;
    size_t count = 0;
    for (const auto& mum : mum_list) {
        output_file << ">mum_" << count << "\n";
        output_file << ref_seq.substr(mum.first, mum.second);
        if (add_terminator) {
            output_file << "#";
        }
        output_file << "\n";
        count++;
    }
    output_file.close();
    
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 1) {
        print_usage();
        return 1;
    }

    std::string mum_file_path;
    std::string output_path = "";
    bool add_terminator = true;  // default is true
    std::string length_file_path = "";

    static struct option long_options[] = {
        {"help",          no_argument,       NULL, 'h'},
        {"mums",          required_argument, NULL, 'm'},
        {"lengths",       required_argument, NULL, 'l'},
        {"output",        required_argument, NULL, 'o'},
        {"no-terminator", no_argument,       NULL, 't'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "hm:l:o:t", long_options, NULL)) >= 0) {
        switch(c) {
            case 'h': print_usage(); return 0;
            case 'm': mum_file_path.assign(optarg); break;
            case 'l': length_file_path.assign(optarg); break;
            case 'o': output_path.assign(optarg); break;
            case 't': add_terminator = false; break;
            default: print_usage(); return 1;
        }
    }

    // Check for positional argument if mum_file_path is still empty
    if (optind < argc && mum_file_path.empty()) {
        mum_file_path.assign(argv[optind]);
    }

    if (mum_file_path.empty()) {
        std::cerr << "Error: No mum file provided\n";
        return 1;
    }
    if (!endsWith(mum_file_path, ".mums") && !endsWith(mum_file_path, ".bumbl")) {
        mum_file_path += ".mums";
    }
    if (!is_file(mum_file_path)) {
        std::cerr << "Error: Invalid mum file provided: " << mum_file_path << std::endl;
        return 1;
    }
    if (length_file_path.empty()) {
        length_file_path = mum_file_path.substr(0, mum_file_path.find_last_of(".")) + ".lengths";
    }
    if (!is_file(length_file_path)) {
        std::cerr << "Error: Invalid lengths file: " << length_file_path << std::endl;
        std::cerr << "Provide a lengths file with the -l option or ensure the mum file has an associated lengths file" << std::endl;
        return 1;
    }

    if (output_path.empty()) {
        output_path = mum_file_path.substr(0, mum_file_path.find_last_of(".")) + "_mums.fa";
    }
    else if (!endsWith(output_path, ".fa")) {
        output_path += ".fa";
    }

    return extract_mums(mum_file_path, length_file_path, output_path, add_terminator);
}
