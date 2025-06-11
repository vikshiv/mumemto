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
   \date 01/07/2020
*/


#ifndef _GSACAK_FILE_HH
#define _GSACAK_FILE_HH

#include <common.hpp>

#include <kseq.h>

#include <ref_builder.hpp>
extern "C" {
    #include<gsacak.h>
}

KSEQ_INIT(int, read);

class gsacak_lcp{
public:
    RefBuilder* ref_build;
    std::vector<uint8_t> text;
    std::vector<uint_t> sa;
    std::vector<int_t> lcp;
    std::vector<uint8_t> bwt;

    gsacak_lcp(std::string filename, RefBuilder* ref_build, bool write_arrays = false) : 
                ref_build(ref_build)
    {
        text.reserve(ref_build->total_length);
        readFasta(filename + ".fna");
        text.push_back(1);
        text.push_back(0);

        sa.resize(text.size());
        lcp.resize(text.size());

        gsacak(&text[0], &sa[0], &lcp[0], nullptr, text.size());
        
        bwt.resize(text.size());
        for (size_t i = 0; i < text.size(); ++i) {
            bwt[i] = text[(sa[i] == 0 ? sa.size() : sa[i]) - 1];
        }

        if (write_arrays) {
            // Opening output files
            FILE *sa_file;
            std::string outfile = filename + ".sa";
            if ((sa_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

            FILE *lcp_file;
            outfile = filename + ".lcp";
            if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

            for (size_t i = 0; i < sa.size(); ++i) {
                if (fwrite(&sa[i], SSABYTES, 1, sa_file) != 1)
                    error("SA write error");
                if (fwrite(&lcp[i], THRBYTES, 1, lcp_file) != 1)
                    error("LCP write error");
            }
            fclose(sa_file);
            fclose(lcp_file);

            outfile = filename + std::string(".bwt");
            write_file(outfile.c_str(), bwt);
        }    
    }

    template <class T>
    size_t process(T &match_finder) {
        size_t count = 0;
        size_t doc_i;
        for (size_t j = 0; j < text.size(); j++)
        {    
            if (j % (text.size() / PBWIDTH) == 0){
                printProgress((double) j / text.size());
            }
                        
            // Start of MUM computation code
            doc_i = ref_build->doc_ends_rank(sa[j]);
            count += match_finder.update(j, bwt[j], doc_i, sa[j], lcp[j]);
            // End of MUM computation code
        }
        printProgress(1.0);
        return count;
    }

private:
    void readFasta(std::string filename) {
        FILE* fp; kseq_t* seq;
        fp = fopen(filename.data(), "r"); 
        if(fp == 0) {std::exit(1);}
        seq = kseq_init(fileno(fp));
        while (kseq_read(seq)>=0) {
            this->text.insert(this->text.end(), seq->seq.s, seq->seq.s + seq->seq.l);
        }
        kseq_destroy(seq);
        fclose(fp);
    }
};

#endif /* end of include guard: _GSACAK_FILE_HH */