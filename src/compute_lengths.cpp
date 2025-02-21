/*
 * File: compute_lengths.cpp
 * Description: Standalone executable for computing sequence lengths from FASTA files
 *              without writing sequences to file.
 */

#include <iostream>
#include <string>
#include <getopt.h>
#include <filesystem>
#include <kseq.h>
#include <unistd.h>
#include <vector>
#include <fstream>


KSEQ_INIT(int, read);

void print_usage() {
    std::fprintf(stderr, "\ncompute_lengths - compute lengths of sequences in FASTA files\n");
    std::fprintf(stderr, "Usage: compute_lengths [options] [input_fastas [...]]\n\n");
    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-32sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-22s%-10spath to a file-list of genomes to use (overrides positional args)\n", "-i, --input", "[FILE]");
    std::fprintf(stderr, "\t%-22s%-10soutput prefix path\n", "-o, --output", "[PREFIX]");
}

bool endsWith(const std::string& str, const std::string& suffix) {
    // Checks if the string ends the suffix
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

int is_file(std::string path) {
    /* Checks if the path is a valid file-path */
    std::ifstream test_file(path.data());
    if (test_file.fail()) {return 0;}
    
    test_file.close();
    return 1;
}

void write_length_files(const std::vector<std::string>& files,
                        const std::vector<std::vector<size_t>>& multifasta_lengths,
                        const std::vector<std::vector<std::string>>& multifasta_names,
                        const std::string& output_prefix,
                        bool multi) {
    size_t total_file_length = 0;
    // Write main lengths file
    std::string lengths_fname = output_prefix + ".lengths";
    std::ofstream outfile(lengths_fname);

    if (multi) {
        for (size_t i = 0; i < multifasta_lengths.size(); ++i) {
            total_file_length = 0;
            for (auto n : multifasta_lengths[i]) {
                total_file_length += n;
            }
            outfile << std::filesystem::absolute(files[i]).string() << " * " << total_file_length << std::endl;
            for (auto idx = 0; idx < multifasta_lengths[i].size(); ++idx) {
                outfile << std::filesystem::absolute(files[i]).string() << " " << multifasta_names[i][idx] << " " << multifasta_lengths[i][idx] << std::endl;
            }
        }
    } else {
        for (size_t i = 0; i < multifasta_lengths.size(); ++i) {
            outfile << std::filesystem::absolute(files[i]).string() << " " << multifasta_lengths[i][0] << std::endl;
        }
    }
    outfile.close();
}

int compute_lengths(std::vector<std::string>& input_files, std::string output_prefix) {
    FILE* fp;
    kseq_t* seq;
    std::vector<std::vector<size_t>> multifasta_lengths;
    std::vector<std::vector<std::string>> multifasta_names;
    std::vector<size_t> temp_lengths;
    std::vector<std::string> temp_names;
    bool multi = false;

    size_t curr_id_seq_length = 0;
    for (const auto& input_file : input_files) {
        temp_lengths.clear();
        temp_names.clear();
        
        fp = fopen(input_file.data(), "r");
        if (fp == 0) {
            std::cerr << "Error opening file: " << input_file << std::endl;
            return 1;
        }

        seq = kseq_init(fileno(fp));
        while (kseq_read(seq) >= 0) {
            curr_id_seq_length += seq->seq.l;
            temp_lengths.push_back(seq->seq.l);
            temp_names.push_back(seq->name.s);
        }

        if (curr_id_seq_length == 0) {
            std::cerr << std::endl << "Error: Empty input file found: " << input_file << std::endl;
            kseq_destroy(seq);
            fclose(fp);
            return 1;
        }
        if (temp_lengths.size() > 1) {
            multi = true;
        }
        multifasta_lengths.push_back(temp_lengths);
        multifasta_names.push_back(temp_names);

        kseq_destroy(seq);
        fclose(fp);
    }

    // Write lengths files
    write_length_files(input_files, multifasta_lengths, multifasta_names, output_prefix, multi);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 1) {
        print_usage();
        return 1;
    }

    std::string input_list = "";
    std::string output_prefix = "output";
    bool use_revcomp = true;
    std::vector<std::string> files;

    static struct option long_options[] = {
        {"help",    no_argument,       NULL, 'h'},
        {"input",   required_argument, NULL, 'i'},
        {"output",  required_argument, NULL, 'o'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "hi:o:", long_options, NULL)) >= 0) {
        switch(c) {
            case 'h': print_usage(); return 0;
            case 'i': input_list.assign(optarg); break;
            case 'o': output_prefix.assign(optarg); break;
            default: print_usage(); return 1;
        }
    }

    for (int i = optind; i < argc; ++i) {
        files.push_back(argv[i]);
    }

    if (input_list.empty() && files.empty()) {
        std::cerr << "Error: No input files provided\n";
        return 1;
    }
    // If no input list provided, create one from files
    if (!input_list.empty()) {
        if (!files.empty()) {
            std::cerr << "Filelist provided, ignoring positional args\n";
            files.clear();
        }
        std::ifstream infile(input_list);
        if (!infile) {
            std::cerr << "Error opening file: " << input_list << std::endl;
            return 1;
        }
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::string file;
            iss >> file;
            if (!file.empty()) {
                files.push_back(file);
            }
        }
        infile.close();
    }

    for (const auto& input_file : files) {
        if (!is_file(input_file)) {
            std::cerr << "The following path in the input list is not valid: " << input_file << std::endl;
            return 1;
        }
        if (!endsWith(input_file, ".fa") && !endsWith(input_file, ".fasta") && !endsWith(input_file, ".fna")) {
            std::cerr << "The following input-file is not a FASTA file: " << input_file << std::endl;
            return 1;
        }
    }

    return compute_lengths(files, output_prefix);
}
