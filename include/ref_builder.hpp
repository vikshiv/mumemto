/*
 * File: ref_builder.hpp
 * Description: Originally written by Omar Ahmed in docprofiles: https://github.com/oma219/docprofiles/tree/main
 *              Header for ref_builder.cpp that includes 
 *              definition of the RefBuilder class.
 * Date: August 31, 2022
 */

#ifndef REF_BUILD_H
#define REF_BUILD_H

#include <string>
#include <sdsl/bit_vectors.hpp>

class RefBuilder {
public:
    std::string input_file = "";
    std::string output_ref = "";
    bool use_revcomp = false;

    sdsl::bit_vector doc_ends;
    sdsl::rank_support_v<1> doc_ends_rank;
    std::vector<size_t> seq_lengths;

    size_t num_docs = 0;
    size_t total_length = 0;
    
    RefBuilder(std::string input_data, std::string output_prefix, bool use_rcomp);
    RefBuilder(std::string output_prefix, bool use_rcomp);
    int build_input_file(size_t w, size_t p, bool probing, bool keep_seqs);
    
    // only used for direct_gsacak
    std::vector<uint8_t> text;

private:
    std::vector<std::string> input_files;
    std::string output_prefix;
    bool from_parse = false;
    std::vector<std::vector<size_t>> multifasta_lengths;
    std::vector<std::vector<std::string>> multifasta_names;
    

    void build_bv();
    void write_lengths_file();
}; // end of RefBuilder class


#endif /* end of REF_BUILD_H */