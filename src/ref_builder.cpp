/*
 * File: ref_builder.cpp
 * Description: Implementation of RefBuilder class which builds
 *              input reference file and determine the number of
 *              bp in each class.
 * Note: this is based on the ref_builder.cpp in the SPUMONI repo (written by Omar Ahmed).
 * Author: Vikram Shivakumar
 * Date: August 31, 2022
 */

#include <ref_builder.hpp>
#include <pfp_mum.hpp> 
#include <newscan.hpp>
#include <zlib.h> 
#include <kseq.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <numeric>
#include <sdsl/bit_vectors.hpp>
#include <unordered_set>
#include <filesystem>

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

RefBuilder::RefBuilder(std::string input_data, std::string output_prefix, 
                        bool use_rcomp): input_file(input_data), use_revcomp(use_rcomp), output_prefix(output_prefix) {
    /* Constructor of RefBuilder - builds input reference and determines size of each class */

    // Verify every file in the filelist is valid
    std::string line = "";
    size_t curr_id = 0, member_num = 0;

    std::ifstream input_fd (input_file.data(), std::ifstream::in);

    while (std::getline(input_fd, line)) {
        auto word_list = split(line, ' ');
        if (word_list.size() == 0) 
            continue;
        // Make sure the filelist has at least 2 columns (name and doc_id)
        // ASSERT((word_list.size() >= 2), "Input file-list does not have expected structure.");
        // if (word_list.size() == 1)
        //     word_list.push_back(std::to_string(curr_id + 1));

        if (!is_file(word_list[0])) {
            FATAL_ERROR("The following path in the input list is not valid: %s", word_list[0].data());}
        if (!endsWith(word_list[0], ".fa") && !endsWith(word_list[0], ".fasta") && !endsWith(word_list[0], ".fna") &&
            !endsWith(word_list[0], ".fa.gz") && !endsWith(word_list[0], ".fasta.gz") && !endsWith(word_list[0], ".fna.gz")) {
            FATAL_ERROR("The following input-file is not a FASTA file: %s", word_list[0].data());}
        input_files.push_back(word_list[0]);

        // Make sure second column is valid integer, and starts at 1
        // if (!is_integer(word_list[1])) 
        //     FATAL_ERROR("A document ID in the file_list is not an integer: %s", word_list[1].data());  
        // if (member_num == 0 && std::stoi(word_list[1]) != 1) 
        //     FATAL_ERROR("The first ID in file_list must be 1");
            
        // if (std::stoi(word_list[1]) ==  static_cast<int>(curr_id) || std::stoi(word_list[1]) == static_cast<int>(curr_id+1)) {
        //     if (std::stoi(word_list[1]) == static_cast<int>(curr_id+1)) 
        //         {
        //             curr_id+=1;
        //         }
        //     document_ids.push_back(curr_id);
        // } else
        //     FATAL_ERROR("The IDs in the file_list must be staying constant or increasing by 1.");
        // member_num += 1;
    }

    // Remove duplicates while preserving order
    std::vector<std::string> unique_files;
    std::unordered_set<std::string> seen;
    for (const auto& file : input_files) {
        if (seen.insert(file).second) {
            unique_files.push_back(file);
        }
    }
    input_files = std::move(unique_files);
    
    // Make sure we have parsed each line, and it has multiple groups
    // ASSERT((document_ids.size() == input_files.size()), "Issue with file-list parsing occurred.");
    if (input_files.size() <= 1) {
        FATAL_ERROR("Multiple FASTA inputs required. Perhaps split a multi-FASTA into multiple files?");}

    this->num_docs = input_files.size();
}

RefBuilder::RefBuilder(std::string output_prefix, bool use_rcomp): use_revcomp(use_rcomp), output_prefix(output_prefix) {
    /* Alternative constructor for running from the lengths file */
    from_parse = true;
    std::string lengths_fname = output_prefix + ".lengths";
    if (!is_file(lengths_fname)) {
        FATAL_ERROR("Lengths file required for using intermediate PFP files. File should match parse prefix: %s", lengths_fname.data());}
    std::ifstream lengths_fd(lengths_fname.data(), std::ifstream::in);
    std::string line;
    size_t cur_length = 0;
    while (std::getline(lengths_fd, line)) {
        auto word_list = split(line, ' ');
        if (word_list.size() == 0)
            continue;
        
        if (word_list.size() == 2) {
            cur_length = std::stoull(word_list[1]) + 1;
        } else if ((word_list.size() == 3) && (word_list[1] == "*")) {
            cur_length = std::stoull(word_list[2]) + 1;
        }
        else {continue;}
        input_files.push_back(word_list[0]);
        
        if (use_revcomp) {
            cur_length *= 2;
        }
        seq_lengths.push_back(cur_length);
    }

    this->num_docs = input_files.size();
}

void RefBuilder::build_bv() {
    // Assuming seq_lengths has been built correctly
    size_t total_input_length = 0;
    for (auto length: seq_lengths) {
        total_input_length += length;
    }

    this->total_length = total_input_length;
    // sanity check
    ASSERT((this->num_docs == seq_lengths.size()), "Issue with file-list parsing occurred.");
    
    // Build bitvector/rank support marking the end of each document
    doc_ends = sdsl::bit_vector(total_input_length, 0);
    size_t curr_sum = 0;

    for (size_t i = 0; i < seq_lengths.size(); i++) {
        curr_sum += seq_lengths[i];
        doc_ends[curr_sum-1] = 1;
    }
    doc_ends_rank = sdsl::rank_support_v<1> (&doc_ends); 
}

void RefBuilder::write_lengths_file(bool multi) {
    // assuming multifasta_lengths is built correctly
    std::string lengths_fname = output_prefix + ".lengths";
    std::ofstream outfile(lengths_fname);
    if (multi) {
        size_t total_file_length = 0;
        for (size_t i = 0; i < multifasta_lengths.size(); ++i) {
            total_file_length = 0;
            for (auto n : multifasta_lengths[i]) {
                total_file_length += n;
            }
            outfile << std::filesystem::absolute(input_files[i]).string() << " * " << total_file_length << std::endl;
            for (auto idx = 0; idx < multifasta_lengths[i].size(); ++idx) {
                outfile << std::filesystem::absolute(input_files[i]).string() << " " << multifasta_names[i][idx] << " " << multifasta_lengths[i][idx] << std::endl;
            }
        }
    } else {
        for (size_t i = 0; i < multifasta_lengths.size(); ++i) {
            outfile << std::filesystem::absolute(input_files[i]).string() << " " << multifasta_lengths[i][0] << std::endl;
        }
    }
    outfile.close();
}

int RefBuilder::build_input_file(size_t w = 10, size_t p = 100, bool probing = false, bool keep_seqs = false) {
    bool multi = false;
    if (!from_parse) {
        // Declare needed parameters for reading/writing
        output_ref = this->output_prefix + ".fna";
        pfparser parser(output_ref, w, p, probing);
        gzFile gzfp; kseq_t* seq;
        std::vector<std::string> seq_vec;
        std::vector<size_t> temp_lengths;
        std::vector<std::string> temp_names;
        
        // Start working on building the reference file by reading each file ...
        size_t curr_id = 1;
        size_t curr_id_seq_length = 0;
        for (auto iter = input_files.begin(); iter != input_files.end(); ++iter) {
            // Use gzopen for both compressed and uncompressed files
            gzfp = gzopen((*iter).data(), "r");
            if (gzfp == 0) {std::exit(1);}
            seq = kseq_init(gzfp);
            size_t iter_index = static_cast<size_t>(iter-input_files.begin());

            while (kseq_read(seq)>=0) {
                // Get forward seq, and write to file
                for (size_t i = 0; i < seq->seq.l; ++i) {
                    seq->seq.s[i] = static_cast<char>(std::toupper(seq->seq.s[i]));
                }
                seq_vec.push_back(seq->seq.s);
                curr_id_seq_length += seq->seq.l;
                temp_lengths.push_back(seq->seq.l);
                temp_names.push_back(seq->name.s);
            }
            if (temp_lengths.size() > 1) {
                multi = true;
            }
            multifasta_lengths.push_back(temp_lengths);
            multifasta_names.push_back(temp_names);
            temp_lengths.clear();
            temp_names.clear();

            kseq_destroy(seq);
            gzclose(gzfp);
            if (curr_id_seq_length == 0) {
                std::cerr << std::endl << "Empty input file found: " + *iter << std::endl;
                return 1;
            }

            // for the terminating $
            curr_id_seq_length += 1;
            if (keep_seqs) {
                std::cout << "Keeping sequences" << std::endl;
                this->text.reserve(this->text.size() + curr_id_seq_length);
                for (auto i = 0; i < seq_vec.size(); ++i) {
                    this->text.insert(this->text.end(), seq_vec.at(i).begin(), seq_vec.at(i).end());
                }
                this->text.push_back('$');
            }
            else {
                std::cout << "Not keeping sequences" << std::endl;
                for (auto i = 0; i < seq_vec.size(); ++i) {
                    parser.process_string(seq_vec.at(i));
                }
                parser.process_string("$");
            }
            // Get reverse complement, and print it
            // Based on seqtk reverse complement code, that does it 
            // in place. (https://github.com/lh3/seqtk/blob/master/seqtk.c)
            if (use_revcomp) {
                for (auto i = seq_vec.size(); i-- != 0; ) {
                    rev_comp(seq_vec.at(i));
                    curr_id_seq_length += seq_vec.at(i).length();
                }
                // for the terminating $
                curr_id_seq_length += 1;
                if (keep_seqs) {
                    this->text.reserve(this->text.size() + curr_id_seq_length);
                    for (auto i = seq_vec.size(); i-- != 0; ) {
                        this->text.insert(this->text.end(), seq_vec.at(i).begin(), seq_vec.at(i).end());
                    }
                    this->text.push_back('$');
                }
                else {
                    for (auto i = seq_vec.size(); i-- != 0; ) {
                        parser.process_string(seq_vec.at(i));
                    }
                    parser.process_string("$");
                }
            }

            seq_lengths.push_back(curr_id_seq_length);
            curr_id += 1; curr_id_seq_length = 0;
            seq_vec.clear();
        }

        // Write final phrase to file, sort dictionary, and remap final parse file
        if (!keep_seqs) {
            parser.finish_parse();
        }
        // Write out lengths file
        this->write_lengths_file(multi);
    }

    // build the bv
    this->build_bv();

    exit(0);

    return 0;
}


