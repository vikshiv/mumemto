/* mumento - finding maximal matches with PFP
    Copyright (C) 2024 Vikram Shivakumar

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
   \file mum_finder.hpp
   \brief mum_finder.hpp compute Maximal Unique Matches (MUMs) between sequences 
   \author Vikram Shivakumar
   \date 12/20/2023
*/


#ifndef _MUM_HH
#define _MUM_HH

#include <common.hpp>
#include <ref_builder.hpp>

#include <unordered_map>
#include <deque>
#include <boost/circular_buffer.hpp>

// build a struct for unique-counter-like types
struct unique_counter {
    std::vector<std::vector<int>> windows;
    std::vector<int> total_unique;
    boost::circular_buffer<size_t> sliding_window;
    size_t num_docs;
    unique_counter(size_t unique, size_t num_docs, size_t topk) : 
        sliding_window(num_docs), 
        total_unique(topk), 
        windows(topk, std::vector<int>(unique,0))
    {
        this->num_docs = num_docs;
    }
    void add(size_t d)
    {
        if(windows[0][d])
        {
            windows[0][d]++;
        }
        else
        {
            windows[0][d] = 1;
            total_unique[0]++;
        }
        for(int i = 1; i < total_unique.size(); i++){
            if (sliding_window.size() < i)
                break;
            d = sliding_window[sliding_window.size() - i];
            if(windows[i][d])
            {
                windows[i][d]++;
            }
            else
            {
                windows[i][d] = 1;
                total_unique[i]++;
            }
        }
    }
    void remove(size_t d)
    {
        // assert(windows[d] > 0);
        for(int i=0; i<windows.size(); i++){
            if(windows[i][d] == 1)
            {
                windows[i][d]--;
                total_unique[i]--;
            }
            else
            {
                windows[i][d]--;
            }
        }
    }
};

// struct to track the minimum value in a sliding window
// Alternate implementation with the linear time RMQ window algorithm,
// but naively a sliding window and min function
struct pq_window {
    // boost::circular_buffer<std::pair<size_t, size_t>> pq;
    // std::deque<std::pair<size_t, size_t>> pq;
    boost::circular_buffer<size_t> sliding_window;
    size_t num_docs;
    size_t left_lcp;
    pq_window(size_t n) : sliding_window(n)//,  pq(n)
    {
        num_docs = n;
    }
    void update(size_t lcp){
        sliding_window.push_back(lcp);
        left_lcp = sliding_window.front();

        // int j needed as arg for below
        // while(!pq.empty() && pq.front().second <= (j + 1 - num_docs))
        //     pq.pop_front();
        // while(!pq.empty() && pq.back().first > lcp)
        //     pq.pop_back();
        // pq.push_back(std::pair<size_t, size_t>(lcp, j));
    }
    size_t min(int i)
    {
        return *std::min_element(std::next(sliding_window.begin()), sliding_window.end() - i);
        // return pq.front().first;
    }
}; 

class mum_finder{
public:

    size_t MIN_MUM_LENGTH;
    size_t num_docs = 0;
    std::vector<size_t> doc_offsets;
    std::vector<size_t> doc_lens;
    size_t topk;
    bool overlap_mum;
    bool col_mums;
    bool revcomp;

    mum_finder(std::string filename, RefBuilder* ref_build, size_t min_mum_len, size_t topk, bool overlap, bool col_mum_mode) : 
                MIN_MUM_LENGTH(min_mum_len),
                num_docs(ref_build->num_docs),
                revcomp(ref_build->use_revcomp),
                doc_lens(ref_build->seq_lengths),
                topk(topk),
                overlap_mum(overlap),
                col_mums(col_mum_mode),
                sa_window(num_docs),
                lcp_pq(num_docs),
                window_docs(num_docs, num_docs, topk),
                window_bwt(6, num_docs, topk),
                doc_offsets(num_docs, 0)
    {
        // get cumulative offset
        size_t curr_sum = 0;
        for (size_t i = 0; i < num_docs - 1; i++) {
            curr_sum += doc_lens[i];
            doc_offsets[i + 1] = curr_sum;
        }
        if (revcomp) {
            for (auto i = 0; i < doc_lens.size(); i++) {
                doc_lens[i] = doc_lens[i] / 2;
            }
        }
        
        // Opening output file
        std::string outfile = filename + ((col_mum_mode) ? std::string(".fa.col_mums") : std::string(".mums"));
        mum_file.open(outfile);
        if (col_mum_mode) {
            size_t docs = num_docs;
            mum_file.write(reinterpret_cast<const char*>(&docs), BWTBYTES);
        }
    }

    void close()
    {
        mum_file.close();
    }

    // main update function, takes in the current streamed value of each array and write mum if found
    inline size_t update(size_t j, uint8_t bwt_c, size_t doc, size_t sa_entry, size_t lcp)
    {
        right_lcp = lcp;
        size_t count = 0;
        bool valid_window = sa_window.size() == num_docs;
        if(valid_window)
        {
            mum_idxs = is_mum();
            if (mum_idxs.size() > 0) {
                if (col_mums) {
                    count = write_mum(j, mum_idxs);
                }
                else {
                    count = write_mum(mum_idxs);
                }
            }

        }
        
        update_lcp_window(right_lcp, valid_window);
        update_sa_window(sa_entry, valid_window);
        update_bwt_window(bwt_c, valid_window);
        update_doc_window(doc, valid_window);
        return count;
    }

private:    

    std::ofstream mum_file;

    // Helper functions and variables to compute MUMs
    //
    // variables needed to track mums
    size_t left_lcp = 0; // lcp of window and suffix preceding window
    size_t right_lcp = 0; // lcp of window and suffix succeeding window

    
    // PQ window to track and compute the minimum LCP value in a window
    pq_window lcp_pq;
    
    // circular buffer for the sliding SA offset window
    boost::circular_buffer<size_t> sa_window;

    // unique counter to track number of unique doc ids in a window
    unique_counter window_docs;

    // track unique bwt chars in the sliding window
    unique_counter window_bwt;

    // map characters to ids, useful to count frequencies in an array
    std::unordered_map<uint8_t,int> nucMap = {
        {'A', 1},
        {'C', 2},
        {'G', 3},
        {'T', 4},
        {'N', 0}, // all default to this too
        {0, 5} //null char
    };

    // stores the current RMQ (avoiding recomputing)
    size_t mum_length = 0;

    // store MUMs to write out
    std::vector<std::pair<int, int>> mum_idxs;

    // min number of iterations to avoid checking for a mum
    // int skip = 0;
    // int total_skips = 0;
    // skip = num_docs - 1;    <- comes after writing mum

    // update the 4 sliding windows with helper functions below

    inline void update_bwt_window(uint8_t bwt_c, bool valid_window)
    {
        int bwt_idx = nucMap[bwt_c];
        window_bwt.add(bwt_idx);
        if(valid_window)
        {
            window_bwt.remove(window_bwt.sliding_window.front());
        }
        window_bwt.sliding_window.push_back(bwt_idx);
    }

    inline void update_doc_window(size_t doc, bool valid_window)
    {
        window_docs.add(doc);
        if(valid_window)
        {
            window_docs.remove(window_docs.sliding_window.front());
        }
        window_docs.sliding_window.push_back(doc);
        // if(num_docs - window_docs.size() > skip)
        //     skip = num_docs - window_docs.size();
    }

    inline void update_sa_window(size_t sa_entry, bool valid_window)
    {
        // slide over sa window
        sa_window.push_back(sa_entry);
    }

    inline void update_lcp_window(size_t lcp, bool valid_window)
    {
        // update pq
        lcp_pq.update(lcp);
    }

    // find the minimum LCP value in the window
    inline size_t rmq_of_window(int i)
    {
        return lcp_pq.min(i);
    }

    // Determine if the windows consitute a MUM by checking the properties
    inline std::vector<std::pair<int, int>> is_mum()
    {
        // Check each condition: (check the fast conditions first, then compute RMQ if needed)

        // pair of window size and mum_length
        std::vector<std::pair<int, int>> mum_subset;
        for(int i = 0; i < topk; i++)
        {
            // Check that every doc appears once
            if(window_docs.total_unique[i] != (num_docs - i))
                continue;
            // Check BWT chars in that range are not all identical (i.e. can be left extended by 1, not maximal)
            if(window_bwt.total_unique[i] == 1)
                continue;
            // check RMQ LCP of window > min_mum
            mum_length = rmq_of_window(i);
            // long enough MUM
            if(mum_length < MIN_MUM_LENGTH)
                continue;
            // Suffix preceding and succeding window don't share long enough prefix (mum is not unique!)
            size_t this_right_lcp;
            if (i == 0)
                this_right_lcp = right_lcp;
            else
                this_right_lcp = lcp_pq.sliding_window[lcp_pq.sliding_window.size() - i];
            if(lcp_pq.left_lcp >= mum_length || this_right_lcp >= mum_length)
                continue;

            mum_subset.push_back(std::pair<int, int>(i, mum_length));

            if (!overlap_mum)
                return mum_subset;
        }
        return mum_subset;
    }

    // write the mum to file
    inline size_t write_mum(size_t const j, std::vector<std::pair<int, int>> const &idxs)
    {
        size_t count = 0;
        int idx;
        size_t mum_length;
        // TODO
        // size_t mum_start = j - num_docs;
        for (auto data : idxs){
            mum_length = data.second;
            mum_file.write(reinterpret_cast<const char*>(&mum_length), BWTBYTES);
            mum_file.write(reinterpret_cast<const char*>(&j), SSABYTES);
            // mum_file.write(reinterpret_cast<const char*>(&sa_window[0]), SSABYTES);
            // mum_file.write(reinterpret_cast<const char*>(&sa_window[num_docs - 1]), SSABYTES);
            count += 1;
        }
        return count;
    }

    // write the mum to file
    inline size_t write_mum(std::vector<std::pair<int, int>> const &idxs)
    {
        std::vector<size_t> offsets(num_docs);
        std::vector<char> strand(num_docs);
        size_t count = 0;
        int doc;
        int idx;
        int mum_length;
        for (auto data : idxs){
            idx = data.first;
            mum_length = data.second;
            std::fill(offsets.begin(), offsets.end(), -1);
            std::fill(strand.begin(), strand.end(), 0);
            for (int i = 0; i < num_docs - idx; i++)
            {
                doc = window_docs.sliding_window[i];
                offsets[doc] = sa_window[i] - doc_offsets[doc];
                if (revcomp && offsets[doc] >= doc_lens[doc]) {
                    strand[doc] = '-';
                    offsets[doc] = offsets[doc] - doc_lens[doc];
                }
                else 
                    strand[doc] = '+';
            }
            // temporarory fix: only write MUMs where 1st genome is + stranded
            if (strand[0] == '-')
                continue;

            mum_file << std::to_string(mum_length) << '\t';
            for (int i = 0; i < num_docs - 1; i++)
            {
                if (offsets[i] == -1) 
                    mum_file << ',';
                else
                    mum_file << offsets[i] << ',';
            }
            if (offsets[num_docs - 1] == -1) 
                mum_file << '\t';
            else
                mum_file << offsets[num_docs - 1] << '\t';

            for (int i = 0; i < num_docs - 1; i++) {
                if (offsets[i] == -1) 
                    mum_file << ',';
                else
                    mum_file << strand[i] << ',';
            }   
            if (offsets[num_docs - 1] == -1) 
                mum_file << std::endl;
            else
                mum_file << strand[num_docs - 1] << std::endl;
            
            count += 1;
        }
        return count;
    }
};

#endif /* end of include guard: _MUM_HH */
