/*
 * File: compute_lengths.cpp
 * Description: Standalone executable for computing sequence lengths from FASTA files
 *              without writing sequences to file.
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

KSEQ_INIT(gzFile, gzread);

/* Complement Table from: https://github.com/lh3/seqtk/blob/master/seqtk.c */
char comp_tab[] = {
    0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
   16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
   32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
   48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
   64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
  'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
   64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
  'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

void rev_comp(std::string &seq) {
  int c0, c1;
  for (size_t i = 0; i < seq.length()>>1; ++i) { // reverse complement sequence
      c0 = comp_tab[(int)seq[i]];
      c1 = comp_tab[(int)seq[seq.length() - 1 - i]];
      seq[i] = c1;
      seq[seq.length() - 1 - i] = c0;
  }
  if (seq.length() & 1) // complement the remaining base
      seq[seq.length()>>1] = comp_tab[static_cast<int>(seq[seq.length()>>1])];
}

void print_usage() {
    std::fprintf(stderr, "\ncompute_lengths - compute lengths of sequences in FASTA files\n");
    std::fprintf(stderr, "Usage: compute_lengths [options] [input_fastas [...]]\n\n");
    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-32sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-22s%-10spath to a file-list of genomes to use (overrides positional args)\n", "-i, --input", "[FILE]");
    std::fprintf(stderr, "\t%-22s%-10soutput prefix path\n", "-o, --output", "[PREFIX]");
    std::fprintf(stderr, "\t%-22s%-10swrite concatenated FASTA files (fwd$rc$)\n", "-p, --print-fasta", "");
}

bool endsWith(const std::string& str, const std::string& suffix) {
    // Checks if the string ends the suffix
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

int is_file(std::string path) {
    /* Checks if the path is a valid file-path */
    return std::filesystem::is_regular_file(path);
}

void write_length_files(const std::vector<std::string>& files,
                        const std::vector<std::vector<size_t>>& multifasta_lengths,
                        const std::vector<std::vector<std::string>>& multifasta_names,
                        const std::string& output_prefix) {
    size_t total_file_length = 0;
    // Write main lengths file
    std::string lengths_fname = output_prefix + ".lengths";
    std::ofstream outfile(lengths_fname);

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
    outfile.close();
}


int compute_lengths(std::vector<std::string>& input_files, std::string output_prefix, bool print_fasta) {
    gzFile gzfp;
    kseq_t* seq;
    std::vector<std::vector<size_t>> multifasta_lengths;
    std::vector<std::vector<std::string>> multifasta_names;
    std::vector<size_t> temp_lengths;
    std::vector<std::string> temp_names;

    size_t curr_id_seq_length = 0;
    for (size_t i = 0; i < input_files.size(); ++i) {
        const auto& input_file = input_files[i];
        temp_lengths.clear();
        temp_names.clear();
        
        // If printing FASTA, prepare output filename
        std::string output_filename;
        if (print_fasta) {
            output_filename = output_prefix + "_file" + std::to_string(i + 1) + ".fna";
        }
        
        // Use gzopen for both compressed and uncompressed files
        gzfp = gzopen(input_file.data(), "r");
        if (gzfp == 0) {
            std::cerr << "Error opening file: " << input_file << std::endl;
            return 1;
        }
        seq = kseq_init(gzfp);
        
        // For FASTA output, we need to concatenate sequences
        std::string concatenated_seq;
        if (print_fasta) {
            while (kseq_read(seq) >= 0) {
                // Convert to uppercase and append to concatenated sequence
                std::string sequence = seq->seq.s;
                for (size_t k = 0; k < sequence.length(); ++k) {
                    sequence[k] = static_cast<char>(std::toupper(sequence[k]));
                }
                concatenated_seq += sequence;
                
                // Also collect length and name info
                curr_id_seq_length += sequence.length();
                temp_lengths.push_back(sequence.length());
                temp_names.push_back(seq->name.s);
            }
        } else {
            // Just collect length and name info
            while (kseq_read(seq) >= 0) {
                curr_id_seq_length += seq->seq.l;
                temp_lengths.push_back(seq->seq.l);
                temp_names.push_back(seq->name.s);
            }
        }

        if (curr_id_seq_length == 0) {
            std::cerr << std::endl << "Error: Empty input file found: " << input_file << std::endl;
            kseq_destroy(seq);
            gzclose(gzfp);
            return 1;
        }
        multifasta_lengths.push_back(temp_lengths);
        multifasta_names.push_back(temp_names);

        // Write FASTA file if requested
        if (print_fasta) {
            std::ofstream outfile(output_filename);

            // Write single sequence: fwd$rc$
            outfile << ">file" << (i + 1) << std::endl;
            outfile << concatenated_seq << "$";
            // Create reverse complement
            rev_comp(concatenated_seq);
            outfile << concatenated_seq << "$" << std::endl;
            outfile.close();
        }

        kseq_destroy(seq);
        gzclose(gzfp);
    }

    // Write lengths files
    write_length_files(input_files, multifasta_lengths, multifasta_names, output_prefix);
    
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 1) {
        print_usage();
        return 1;
    }

    std::string input_list = "";
    std::string output_prefix = "output";
    bool print_fasta = false;
    std::vector<std::string> files;

    static struct option long_options[] = {
        {"help",    no_argument,       NULL, 'h'},
        {"input",   required_argument, NULL, 'i'},
        {"output",  required_argument, NULL, 'o'},
        {"print-fasta", no_argument,   NULL, 'p'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "hi:o:p", long_options, NULL)) >= 0) {
        switch(c) {
            case 'h': print_usage(); return 0;
            case 'i': input_list.assign(optarg); break;
            case 'o': output_prefix.assign(optarg); break;
            case 'p': print_fasta = true; break;
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
        if (!endsWith(input_file, ".fa") && !endsWith(input_file, ".fasta") && !endsWith(input_file, ".fna") &&
            !endsWith(input_file, ".fa.gz") && !endsWith(input_file, ".fasta.gz") && !endsWith(input_file, ".fna.gz")) {
            std::cerr << "The following input-file is not a FASTA file: " << input_file << std::endl;
            return 1;
        }
    }

    return compute_lengths(files, output_prefix, print_fasta);
}
