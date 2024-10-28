/*
 * File: ref_builder.cpp
 * Description: Implementation of RefBuilder class which builds
 *              input reference file and determine the number of
 *              bp in each class.
 * Note: this is based on the ref_builder.cpp in the SPUMONI repo.
 * Author: Omar Ahmed
 * Date: August 31, 2022
 */

#include <ref_builder.hpp>
#include <pfp_mum.hpp> 
#include <zlib.h> 
#include <kseq.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <numeric>
#include <algorithm>
#include <sdsl/bit_vectors.hpp>

KSEQ_INIT(int, read);

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

char process(char c) {
    c = static_cast<char>(std::toupper(c));
    if (c == 'N') {c = 'A';}
    return c;
}

RefBuilder::RefBuilder(std::string input_data, std::string output_prefix, 
                        bool use_rcomp): input_file(input_data), use_revcomp(use_rcomp) {
    /* Constructor of RefBuilder - builds input reference and determines size of each class */

    // Verify every file in the filelist is valid
    std::string line = "";
    size_t curr_id = 0, member_num = 0;

    std::ifstream input_fd (input_file.data(), std::ifstream::in);
    std::vector<std::string> input_files;
    std::vector<size_t> document_ids;

    while (std::getline(input_fd, line)) {
        auto word_list = split(line, ' ');

        // Make sure the filelist has at least 2 columns (name and doc_id)
        ASSERT((word_list.size() >= 2), "Input file-list does not have expected structure.");

        if (!is_file(word_list[0])) {
            FATAL_ERROR("The following path in the input list is not valid: %s", word_list[0].data());}
        if (!endsWith(word_list[0], ".fa") && !endsWith(word_list[0], ".fasta") && !endsWith(word_list[0], ".fna")) {
            FATAL_ERROR("The following input-file is not a FASTA file: %s", word_list[0].data());}
        input_files.push_back(word_list[0]);

        // Make sure second column is valid integer, and starts at 1
        if (!is_integer(word_list[1])) 
            FATAL_ERROR("A document ID in the file_list is not an integer: %s", word_list[1].data());  
        if (member_num == 0 && std::stoi(word_list[1]) != 1) 
            FATAL_ERROR("The first ID in file_list must be 1");
            
        if (std::stoi(word_list[1]) ==  static_cast<int>(curr_id) || std::stoi(word_list[1]) == static_cast<int>(curr_id+1)) {
            if (std::stoi(word_list[1]) == static_cast<int>(curr_id+1)) {curr_id+=1;}
            document_ids.push_back(static_cast<size_t>(std::stoi(word_list[1])));
        } else
            FATAL_ERROR("The IDs in the file_list must be staying constant or increasing by 1.");
        member_num += 1;
    }
    
    // Make sure we have parsed each line, and it has multiple groups
    ASSERT((document_ids.size() == input_files.size()), "Issue with file-list parsing occurred.");
    if (document_ids.back() == 1) {
        FATAL_ERROR("If you only have one class ID, you should not build a document array.");}

    // Declare needed parameters for reading/writing
    output_ref = output_prefix + ".fa";
    std::ofstream output_fd (output_ref.data(), std::ofstream::out);
    FILE* fp; kseq_t* seq;
    // std::vector<size_t> seq_lengths;

    // Start working on building the reference file by reading each file ...
    curr_id = 1;
    size_t curr_id_seq_length = 0;
    for (auto iter = input_files.begin(); iter != input_files.end(); ++iter) {

        fp = fopen((*iter).data(), "r"); 
        if(fp == 0) {std::exit(1);}

        seq = kseq_init(fileno(fp));
        size_t iter_index = static_cast<size_t>(iter-input_files.begin());

        while (kseq_read(seq)>=0) {
            // Get forward seq, and write to file
			for (size_t i = 0; i < seq->seq.l; ++i) {
				seq->seq.s[i] = process(seq->seq.s[i]);
            }
           
            // Added dollar sign as separator, and added 1 to length
            output_fd << '>' << seq->name.s << '\n' << seq->seq.s << '\n';
            curr_id_seq_length += seq->seq.l;

            // Get reverse complement, and print it
            // Based on seqtk reverse complement code, that does it 
            // in place. (https://github.com/lh3/seqtk/blob/master/seqtk.c)
            if (use_revcomp) {
                int c0, c1;
                for (size_t i = 0; i < seq->seq.l>>1; ++i) { // reverse complement sequence
                    c0 = comp_tab[(int)seq->seq.s[i]];
                    c1 = comp_tab[(int)seq->seq.s[seq->seq.l - 1 - i]];
                    seq->seq.s[i] = c1;
                    seq->seq.s[seq->seq.l - 1 - i] = c0;
                }
                if (seq->seq.l & 1) // complement the remaining base
                    seq->seq.s[seq->seq.l>>1] = comp_tab[static_cast<int>(seq->seq.s[seq->seq.l>>1])];

                // Added dollar sign as separator, and added 1 to length
                output_fd << '>' << seq->name.s << "_rev_comp" << '\n' << seq->seq.s << '\n';
                curr_id_seq_length += seq->seq.l;
            }
        }
        kseq_destroy(seq);
        fclose(fp);

        // Check if we are transitioning to a new group
        if (iter_index < document_ids.size()-1 && document_ids[iter_index] != document_ids[iter_index+1]){
            seq_lengths.push_back(curr_id_seq_length);
            curr_id += 1; curr_id_seq_length = 0;
        // If it is the last file, output current sequence length
        } else if (iter_index == document_ids.size()-1){
            seq_lengths.push_back(curr_id_seq_length);
            curr_id_seq_length = 0;
        }
    }
    output_fd.close();

    // Add 1 to last document for $ and find total length
    size_t total_input_length = 0;
    seq_lengths[seq_lengths.size()-1] += 1; // for $
    for (auto length: seq_lengths) {
        total_input_length += length;
    }

    this->total_length = total_input_length;
    this->num_docs = seq_lengths.size();
    
    // Build bitvector/rank support marking the end of each document
    doc_ends = sdsl::bit_vector(total_input_length, 0);
    size_t curr_sum = 0;

    for (size_t i = 0; i < seq_lengths.size(); i++) {
        curr_sum += seq_lengths[i];
        doc_ends[curr_sum-1] = 1;
    }
    doc_ends_rank = sdsl::rank_support_v<1> (&doc_ends); 
}