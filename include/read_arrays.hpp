/* pfp_lcp - lcp from prefix free parsing 
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file pfp_lcp.hpp
   \brief pfp_lcp.hpp define and build the lcp from the prefix-free parsing.
            Adapted by Vikram Shivakumar to compute Maximal Unique/Exact Matches (MUM/MEM) between sequences (12/20/2023)
   \author Massimiliano Rossi
   \date 01/07/2020
*/


#ifndef _LCP_FILE_HH
#define _LCP_FILE_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <ref_builder.hpp>

class file_lcp{
public:
    RefBuilder* ref_build;

    file_lcp(std::string filename, RefBuilder* ref_build) : 
                ref_build(ref_build),
                SA(filename + std::string(".sa"), std::ios::binary),
                LCP(filename + std::string(".lcp"), std::ios::binary),
                BWT(filename + std::string(".bwt"), std::ios::binary)
    {
        // Opening input files
        std::string outfile;
        
        if (!SA.is_open()) {
            std::fprintf(stderr, "\nError: opening SA file\n\n");
            std::exit(1);
        }
        if (!BWT.is_open()) {
            std::fprintf(stderr, "\nError: opening BWT file\n\n");
            std::exit(1);
        }
        if (!LCP.is_open()) {
            std::fprintf(stderr, "\nError: opening LCP file\n\n");
            std::exit(1);
        }
    }

    void close() {
        // Close output files
        LCP.close();
        BWT.close();
        SA.close();
        
    }


    template <class T>
    size_t process(T &match_finder) {
        size_t count;
        for (auto j = 0; j < ref_build->total_length; j++)
        {    
            if (j % (ref_build->total_length / PBWIDTH) == 0)
                        printProgress((double) j / ref_build->total_length);
            // Start of MUM computation code
            uint8_t bwt_i = get_BWT_offset();
            size_t sa_i = get_SA_offset();
            size_t lcp_i = get_LCP_offset();
            size_t doc_i = ref_build->doc_ends_rank(sa_i);

            count += match_finder.update(j, bwt_i, doc_i, sa_i, lcp_i, 0);
            // End of MUM computation code
        }
        printProgress(1.0);
        return count;
    }

private:
    std::ifstream LCP;
    std::ifstream BWT;
    std::ifstream SA;

    // read in SA 5 bytes as a low byte and 4 high bytes
    // uint64_t SA_offset;
    uint8_t low_SA;
    uint32_t high_SA;
    uint8_t low_LCP;
    uint32_t high_LCP;

    uint8_t bwt_char;

    inline uint64_t get_SA_offset() {
        // SA.read(reinterpret_cast<char *>(&SA_offset), 5);
        SA.read(reinterpret_cast<char *>(&low_SA), sizeof(low_SA));
        SA.read(reinterpret_cast<char *>(&high_SA), sizeof(high_SA));
        return ((uint64_t)high_SA) << 8 | (uint64_t)low_SA;
    }
    inline uint8_t get_BWT_offset() {
        BWT.read(reinterpret_cast<char *>(&bwt_char), sizeof(bwt_char));
        return bwt_char;
    }
    inline uint64_t get_LCP_offset() {
        LCP.read(reinterpret_cast<char *>(&low_LCP), sizeof(low_LCP));
        LCP.read(reinterpret_cast<char *>(&high_LCP), sizeof(high_LCP));
        return ((uint64_t)high_LCP) << 8 | (uint64_t)low_LCP;
    }
    
};

#endif /* end of include guard: _LCP_FILE_HH */