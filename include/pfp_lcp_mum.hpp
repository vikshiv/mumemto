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


#ifndef _LCP_PFP_HH
#define _LCP_PFP_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C"
{
#include <gsacak.h>
}

#include <ref_builder.hpp>
#include <pfp.hpp>

// Adapted from https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
#define PBSTR "============================================================"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    if (val == 0) {std::fprintf(stderr, "\n\033[33m%3d%%\033[m [%.*s%*s]", val, lpad, PBSTR, rpad, ""); std::fflush(stderr);}
    else if (val == 100) {std::fprintf(stderr, "\r\033[33m%3d%%\033[m [%.*s%*s]", val, lpad, PBSTR, rpad, "");}
    else {std::fprintf(stderr, "\r\033[33m%3d%%\033[m [%.*s>%*s]", val, lpad, PBSTR, rpad-1, ""); std::fflush(stderr);}
    
}

class pfp_lcp{
public:

    pf_parsing& pf;
    std::vector<size_t> min_s; // Value of the minimum lcp_T in each run of BWT_T
    std::vector<size_t> pos_s;    // Position of the minimum lcp_T in each run of BWT_T

    uint8_t head;
    size_t length = 0; // Length of the current run of BWT_T
    RefBuilder* ref_build;
    bool write_arrays;

    size_t num_run;

    pfp_lcp(pf_parsing &pfp_, std::string filename, RefBuilder* ref_build, bool write_arrays) : 
                pf(pfp_),
                min_s(1, pf.n),
                pos_s(1,0),
                head(0),
                num_run(0),
                ref_build(ref_build),
                write_arrays(write_arrays)
    {
        // Opening output files
        if (write_arrays) {
            std::string outfile = filename + std::string(".lcp");
            if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

            outfile = filename + std::string(".sa");
            if ((sa_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

            outfile = filename + std::string(".bwt");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
        }
        assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict);
    }

    void close() {
        // Close output files
        if (write_arrays) {
            fclose(sa_file);
            fclose(bwt_file);
            fclose(lcp_file);
        }
        
    }

    template <class T>
    size_t process(T &match_finder) {
        phrase_suffix_t curr;
        phrase_suffix_t prev;
        size_t count = 0;
        size_t pb_inc = pf.n / PBWIDTH;
        inc(curr);
        while (curr.i < pf.dict.saD.size())
        {

            if(is_valid(curr)){
                // Compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);  // Store the list of all phrase ids with the same suffix.

                phrase_suffix_t next = curr;

                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {

                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_suffix.push_back(next);
                    }

                }

                // Hard case
                int_t lcp_suffix = compute_lcp_suffix(curr,prev);
                
                // really ugly change, but it let's us use int_vector for ilist instead of vector (i.e we can't do pointer arithmetic)
                // typedef std::pair<size_t, std::pair<size_t, uint8_t>> pq_t;
                // type -> ( (idx_start, val_start) , ( (idx_end, val_end), bwt_char) ), replacing the pointer with (idx, val)
                // occ.first.first = idx_start, occ.first.second = val, occ.second.first = idx_end, occ.second.second = bwt
                typedef std::pair<std::pair<size_t, size_t>, std::pair<size_t, uint8_t>> pq_t;

                // using lambda to compare elements.
                auto cmp = [](const pq_t lhs, const pq_t rhs) {
                    return lhs.first.second > rhs.first.second;
                };

                std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
                for (auto s: same_suffix)
                {
                    size_t begin = pf.pars.select_ilist_s(s.phrase + 1);
                    size_t end = pf.pars.select_ilist_s(s.phrase + 2);
                    // pq.push({&pf.pars.ilist[begin], {&pf.pars.ilist[end], s.bwt_char}});
                    pq.push({{begin, pf.pars.ilist[begin]}, {end, s.bwt_char}});
                }

                size_t prev_occ;
                bool first = true;
                while (!pq.empty())
                {
                    if (j % pb_inc == 0)
                        printProgress((double) j / pf.n);
                    auto curr_occ = pq.top();
                    pq.pop();
                    auto curr_occ_first = curr_occ.first.second;
                    if (!first)
                    {
                        // Compute the minimum s_lcpP of the the current and previous occurrence of the phrase in BWT_P
                        lcp_suffix = curr.suffix_length + min_s_lcp_T(curr_occ_first, prev_occ);
                    }
                    first = false;
                    // Update min_s
                    if (write_arrays)
                        print_lcp(lcp_suffix, j);
                    update_ssa(curr, curr_occ_first);
                    if (write_arrays) {
                        print_sa(); 
                    }
                    update_bwt(curr_occ.second.second, 1);

                    // Start of MUM computation code
                    uint8_t curr_bwt_ch = curr_occ.second.second;
                    size_t sa_i = ssa;
                    size_t doc_i = ref_build->doc_ends_rank(ssa);
                    // lcp is in lcp_suffix

                    count += match_finder.update(j, curr_bwt_ch, doc_i, sa_i, lcp_suffix);
                    // End of MUM computation code

                    
                    // Update prevs
                    prev_occ = curr_occ_first;

                    // Update pq
                    curr_occ.first.first++;
                    if (curr_occ.first.first != curr_occ.second.first) {
                        curr_occ.first.second = pf.pars.ilist[curr_occ.first.first];
                        pq.push(curr_occ);
                    }

                    j += 1;
                }

                prev = same_suffix.back();
                
                curr = next;
                
            }
            else
            {
                inc(curr);
            }
        }
        // print last BWT char and SA sample
        if (write_arrays) {
            print_sa();
            print_bwt();
        }
        printProgress(1.0);
        return count;
    }



private:
    typedef struct
    {
        size_t i = 0; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t phrase = 0;
        size_t suffix_length = 0;
        int_da sn = 0;
        uint8_t bwt_char = 0;
    } phrase_suffix_t;

    size_t j = 0;

    size_t ssa = 0;
    size_t esa = 0;

    FILE *lcp_file;

    FILE *bwt_file;

    FILE *sa_file;
    // FILE *esa_file;

    inline bool inc(phrase_suffix_t& s)
    {
        s.i++;
        if (s.i >= pf.dict.saD.size())
            return false;
        s.sn = pf.dict.saD[s.i];
        s.phrase = pf.dict.rank_b_d(s.sn);
        // s.phrase = pf.dict.daD[s.i] + 1; // + 1 because daD is 0-based
        assert(!is_valid(s) || (s.phrase > 0 && s.phrase < pf.pars.ilist.size()));
        s.suffix_length = pf.dict.select_b_d(pf.dict.rank_b_d(s.sn + 1) + 1) - s.sn - 1;
        if(is_valid(s))
            s.bwt_char = (s.sn == pf.w ? 0 : pf.dict.d[s.sn - 1]);
        return true;
    }

    inline bool is_valid(phrase_suffix_t& s)
    {
        // avoid the extra w # at the beginning of the text
        if (s.sn < pf.w)
            return false;
        // Check if the suffix has length at least w and is not the complete phrase.
        if (pf.dict.b_d[s.sn] != 0 || s.suffix_length < pf.w)
            return false;
        
        return true;
    }
    
    inline int_t min_s_lcp_T(size_t left, size_t right)
    {
        // assume left < right
        if (left > right)
            std::swap(left, right);

        assert(pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] >= pf.w);

        return (pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] - pf.w);
    }

    inline int_t compute_lcp_suffix(phrase_suffix_t& curr, phrase_suffix_t& prev)
    {
        int_t lcp_suffix = 0;

        if (j > 0)
        {
            // Compute phrase boundary lcp
            lcp_suffix = pf.dict.lcpD[curr.i];
            for (size_t k = prev.i + 1; k < curr.i; ++k)
            {
                lcp_suffix = std::min(lcp_suffix, pf.dict.lcpD[k]);
            }

            if (lcp_suffix >= curr.suffix_length && curr.suffix_length == prev.suffix_length)
            {
                // Compute the minimum s_lcpP of the phrases following the two phrases
                // we take the first occurrence of the phrase in BWT_P
                size_t left = pf.pars.ilist[pf.pars.select_ilist_s(curr.phrase + 1)]; //size_t left = first_P_BWT_P[phrase];
                // and the last occurrence of the previous phrase in BWT_P
                size_t right = pf.pars.ilist[pf.pars.select_ilist_s(prev.phrase + 2) - 1]; //last_P_BWT_P[prev_phrase];
                
                lcp_suffix += min_s_lcp_T(left,right);
            }
        }

        return lcp_suffix;
    }

    inline void print_lcp(int_t val, size_t pos)
    {
        size_t tmp_val = val;
        if (fwrite(&tmp_val, THRBYTES, 1, lcp_file) != 1)
            error("LCP write error 1");
    }

    // We can put here the check if we want to store the LCP or stream it out
    inline void new_min_s(int_t val, size_t pos)
    {
        min_s.push_back(val);
        pos_s.push_back(j);
    }

    inline void update_ssa(phrase_suffix_t &curr, size_t pos)
    {
        ssa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
        assert(ssa < (pf.n - pf.w + 1ULL));
    }

    inline void update_esa(phrase_suffix_t &curr, size_t pos)
    {
        esa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
        assert(esa < (pf.n - pf.w + 1ULL));
    }

    inline void print_sa()
    {
        if (j < (pf.n - pf.w + 1ULL))
        {
            // std::cout << ssa << std::endl;
            size_t pos = j;
            // if (fwrite(&pos, SSABYTES, 1, sa_file) != 1)
            //     error("SA write error 1");
            if (fwrite(&ssa, SSABYTES, 1, sa_file) != 1)
                error("SA write error");
        }

        // if (j > 0)
        // {
        //     size_t pos = j - 1;
        //     if (fwrite(&pos, SSABYTES, 1, esa_file) != 1)
        //         error("SA write error 1");
        //     if (fwrite(&esa, SSABYTES, 1, esa_file) != 1)
        //         error("SA write error 2");
        // }
    }

    inline void print_bwt()
    {   
        for (size_t i = 0; i < length; ++i)
        {
            if (fputc(head, bwt_file) == EOF)
                error("BWT write error 1");
        }
    }

    inline void update_bwt(uint8_t next_char, size_t length_)
    {
        if (head != next_char)
        {
            if (write_arrays)
                print_bwt();

            head = next_char;

            length = 0;

            num_run++;
        }
        length += length_;

    }

};

#endif /* end of include guard: _LCP_PFP_HH */
