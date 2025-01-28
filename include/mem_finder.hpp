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
   \file mem_finder.hpp
   \brief mem_finder.hpp compute Maximal Exact Matches (MEMs) between sequences 
   \author Vikram Shivakumar
   \date 2/21/2024
*/


#ifndef _MEM_HH
#define _MEM_HH

#include <common.hpp>
#include <ref_builder.hpp>
#include <unordered_set>

#include <deque>

class mem_finder{
public:

    size_t min_mem_length;
    size_t num_distinct;
    int max_freq;
    int max_doc_freq;
    bool no_max_freq;
    size_t num_docs = 0;
    std::vector<size_t> doc_offsets;
    std::vector<size_t> doc_lens;
    bool revcomp;
    bool mummode;

    mem_finder(std::string filename, RefBuilder& ref_build, size_t min_mem_len, size_t num_distinct, int max_doc_freq, int max_total_freq): 
                min_mem_length(min_mem_len),
                num_docs(ref_build.num_docs),
                revcomp(ref_build.use_revcomp),
                doc_lens(ref_build.seq_lengths),
                num_distinct(num_distinct),
                doc_offsets(num_docs, 0),
                max_freq(max_total_freq),
                max_doc_freq(max_doc_freq)
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
        this->mummode = (max_doc_freq == 1);
        // Initialize stack
        init_stack();
        
        // Set parameters and limits
        // this->max_freq = num_docs + max_freq;
        this->no_max_freq = max_freq == 0;

        // Opening output file
        std::string outfile = filename + (mummode ? std::string(".mums") : std::string(".mems"));
        mem_file.open(outfile);
    }

    virtual void close()
    {
        mem_file.close();
    }

    // main update function, takes in the current streamed value of each array and write mem if found
    inline size_t update(size_t j, uint8_t bwt_c, size_t doc, size_t sa_entry, size_t lcp)
    {   
        // bwt last change checker
        size_t count = update_mems(j, lcp);
        if (bwt_buffer.size() == 0 || bwt_buffer.back() != bwt_c)
            last_bwt_change = j;
        update_buffers(j, bwt_c, sa_entry, lcp, doc);
        return count;
    }
protected:
    std::ofstream mem_file;

    // Helper functions and variables to compute MEMs
    
    // Speed up checking BWT property by storing the last change position
    size_t last_bwt_change = 0;

    size_t buffer_start = 0;
    std::deque<size_t> sa_buffer;
    std::deque<uint8_t> bwt_buffer;
    std::deque<size_t> da_buffer;

    inline bool check_bwt_range(size_t start, size_t end) 
    {
        // std::cout << last_bwt_change <<", " << start << ", " << end << std::endl;
        return last_bwt_change <= start;
        // size_t iterations = end - start;
        // size_t idx = 0;
        // std::deque<uint8_t>::iterator it = bwt_buffer.begin() + (start - buffer_start);
        // uint8_t cur_char = *it;
        // while (idx < iterations) {
        //     it++;
        //     if (*it != cur_char)
        //     {
        //         // std::cout <<"good range: "<< last_bwt_change <<", " << start << ", " << end << std::endl;
        //         return false;
        //     }
        //     idx++;
        // }
        // // std::cout <<"bad range: "<< last_bwt_change <<", " << start << ", " << end << std::endl;
        // return true;
    }

    inline size_t write_mem(size_t length, size_t start, size_t end)
    {
        std::string pos = "";
        std::string strand = "";
        std::string docs = "";
        
        size_t curpos;
        size_t curdoc;
        char curstrand;

        // size_t min_offset = -1;
        // char min_offset_strand;

        for (size_t i = start; i < end; i++)
        {
            curdoc = da_buffer.at(i - buffer_start);
            curpos = sa_buffer.at(i - buffer_start) - doc_offsets[curdoc];
            if (revcomp && curpos >= doc_lens[curdoc]) {
                curstrand = '-';
                curpos = doc_lens[curdoc] + doc_lens[curdoc] - curpos - length;
            }
            else
                curstrand = '+';
            
            // if (revcomp && sa_buffer.at(i - buffer_start) < min_offset) {
            //     min_offset = sa_buffer.at(i - buffer_start);
            //     min_offset_strand = curstrand;
            // }

            pos += std::to_string(curpos) + ",";
            docs += std::to_string(curdoc) + ",";
            strand += curstrand;
            strand += ",";
        }
        curdoc = da_buffer.at(end - buffer_start);
        curpos = sa_buffer.at(end - buffer_start) - doc_offsets[curdoc];
        if (revcomp && curpos >= doc_lens[curdoc]) {
            curstrand = '-';
            curpos = doc_lens[curdoc] + doc_lens[curdoc] - curpos - length;
        }
        else
            curstrand = '+';
        pos += std::to_string(curpos);
        docs += std::to_string(curdoc);
        strand += curstrand;

        // if (revcomp && sa_buffer.at(end - buffer_start) < min_offset) {
        //     min_offset = sa_buffer.at(end - buffer_start);
        //     min_offset_strand = curstrand;
        // }
        // if (min_offset_strand == '+')
            mem_file << std::to_string(length) << '\t' << pos << '\t' << docs << '\t' << strand << std::endl;
        return 1;
    }

    inline size_t write_mum(size_t length, size_t start, size_t end)
    {
        std::vector<size_t> offsets(num_docs, -1);
        std::vector<char> strand(num_docs, 0);
        size_t curpos;
        size_t curdoc;
        char curstrand;
        std::string pos_string = "";
        std::string strand_string = "";
        for (size_t i = start; i <= end; i++)
        {
            curdoc = da_buffer.at(i - buffer_start);
            curpos = sa_buffer.at(i - buffer_start) - doc_offsets[curdoc];
            if (revcomp && curpos >= doc_lens[curdoc]) {
                curstrand = '-';
                curpos = doc_lens[curdoc] + doc_lens[curdoc] - curpos - length - 1;
            }
            else
                curstrand = '+';

            offsets[curdoc] = curpos;    
            strand[curdoc] = curstrand;
        }
        // temporarory fix: only write MUMs where 1st genome is + stranded
        int i = 0;
        while (i < num_docs - 1)
        {
            if (strand[i] != 0)
                break;
            i++;
        }
        if (strand[i] == '-')
            return 0;

        for (int i = 0; i < num_docs - 1; i++)
        {
            if (offsets[i] == -1) 
            {
                pos_string += ",";
                strand_string += ",";
            }
            else
            {
                pos_string += std::to_string(offsets[i]) + ",";
                strand_string += strand[i];
                strand_string += ",";
            }
        }
        if (offsets[num_docs - 1] != -1) 
        {
            pos_string += std::to_string(offsets[num_docs - 1]);
            strand_string += strand[num_docs - 1];
        }
        mem_file << std::to_string(length) << '\t' << pos_string << '\t' << strand_string << std::endl;
        return 1;
    }

    inline bool check_doc_range(size_t start, size_t end) 
    {
        std::unordered_map<size_t, size_t> seen;
        size_t unique = 0;
        size_t iterations = end - start + 1;
        size_t idx = 0;
        std::deque<size_t>::iterator it = da_buffer.begin() + (start - buffer_start);
        size_t cur_doc = *it;
        while (idx < iterations) {
            if (!seen.count(cur_doc)) {
                unique++;
                seen[cur_doc] = 1;
                // if (max_doc_freq == 0 && unique >= num_distinct)
                //     return true;
            } else {
                // seen[cur_doc]++;
                if (max_doc_freq && (++seen[cur_doc]) > max_doc_freq)
                    return false;
            }
            it++;
            idx++;
            cur_doc = *it;
        }
        return unique >= num_distinct;
    }
        

private:    

    std::vector<std::pair<size_t, size_t>> current_mems; // list of pairs, (start idx in SA, length of mem)
    
    inline size_t update_mems(size_t j, size_t lcp)
    {
        // three cases for LCP, increase, decrease, or stagnant (nothing changes)
        // j = idx in SA
        size_t count = 0;
        size_t start = j - 1;
        std::pair<size_t, size_t> interval;
        while (lcp < current_mems.back().second) {
            interval = current_mems.back();
            current_mems.pop_back();

            // check conditions of MEM/MUM
            if (interval.second >= min_mem_length && 
                j - interval.first >= num_distinct && 
                (no_max_freq || j - interval.first <= max_freq) &&
                !check_bwt_range(interval.first, j-1) && 
                check_doc_range(interval.first, j-1)) 
                {
                    if (mummode)
                        count += write_mum(interval.second, interval.first, j - 1);
                    else
                        count += write_mem(interval.second, interval.first, j - 1);
                }
            start = interval.first;
        }

        if (lcp > current_mems.back().second) {
            if (lcp >= min_mem_length)
                current_mems.push_back(std::make_pair(start, lcp));
        }

        return count;
    }

    inline void update_buffers(size_t j, uint8_t bwt_c, size_t sa_pos, size_t lcp, size_t docid) {
        if (current_mems.size() <= 1) { // empty stack, only the null interval. < 1 should never happen, but included nonetheless
            if (sa_buffer.size() > 0) {
                sa_buffer.clear();
                bwt_buffer.clear();
                da_buffer.clear();
            }
            buffer_start = j;
        }
        else if (current_mems[1].first > buffer_start) {
            size_t to_remove = current_mems[1].first - buffer_start;
            buffer_start = current_mems[1].first;
            sa_buffer.erase(sa_buffer.begin(), sa_buffer.begin() + to_remove);
            bwt_buffer.erase(bwt_buffer.begin(), bwt_buffer.begin() + to_remove);
            da_buffer.erase(da_buffer.begin(), da_buffer.begin() + to_remove);
        }
        bwt_buffer.push_back(bwt_c);
        sa_buffer.push_back(sa_pos);
        da_buffer.push_back(docid);
    }

    virtual void init_stack()
    {
        current_mems.push_back(std::make_pair(0, 0));
    }

};

class mergable_mem_finder : public mem_finder {
public:
    mergable_mem_finder(std::string filename, RefBuilder& ref_build, size_t min_mem_len, size_t num_distinct, int max_doc_freq, int max_total_freq):
        mem_finder(filename, ref_build, min_mem_len, num_distinct, max_doc_freq, max_total_freq),
        candidate_thresh(doc_lens[0], 0),
        filename(filename)
    {
        init_stack();
    }

    inline size_t update(size_t j, uint8_t bwt_c, size_t doc, size_t sa_entry, size_t lcp)
    {   
        // bwt last change checker
        size_t count = update_mems(j, lcp);
        if (bwt_buffer.size() == 0 || bwt_buffer.back() != bwt_c)
            last_bwt_change = j;
        update_buffers(j, bwt_c, sa_entry, lcp, doc);
        prev_lcp = lcp;
        return count;
    }

    void close() override
    {
        mem_file.close();

        std::string thresh_file = filename + ".thresh";
        std::ofstream thresh_out(thresh_file, std::ios::binary);
        thresh_out.write(reinterpret_cast<const char*>(candidate_thresh.data()), candidate_thresh.size() * sizeof(uint16_t));
        thresh_out.close();
    }
    
protected:
    size_t prev_lcp = 0;

private:
    // Override the data structure for current_mems to include prev_lcp
    std::vector<std::pair<std::pair<size_t, size_t>, size_t>> current_mems;

    // data structure to hold meta data for merging
    std::vector<uint16_t> candidate_thresh;
    size_t MAX_THRESH = static_cast<size_t>(UINT16_MAX);
    std::string filename; // filename for output metadata

    // Override update_mems to handle prev_lcp and next_best
    inline size_t update_mems(size_t j, size_t lcp)
    {
        // three cases for LCP, increase, decrease, or stagnant (nothing changes)
        // j = idx in SA
        size_t count = 0;
        size_t start = j - 1;
        size_t prev = 0;
        size_t next_best = 0;
        size_t start_offset = 0;
        std::pair<size_t, size_t> interval;
        while (lcp < current_mems.back().first.second) {
            interval = current_mems.back().first;
            prev = current_mems.back().second;
            current_mems.pop_back();

            // check conditions of MEM/MUM
            if (interval.second >= min_mem_length && 
                j - interval.first >= num_distinct && 
                (no_max_freq || j - interval.first <= max_freq) &&
                // !check_bwt_range(interval.first, j-1) && 
                check_doc_range(interval.first, j-1)) 
                {
                    next_best = std::min(std::max(prev, lcp), MAX_THRESH); // cap at max uint16
                    for (size_t i = interval.first; i < j-1; i++) {
                        if (da_buffer[i - buffer_start] == 0) {
                            start_offset = sa_buffer[i - buffer_start] - doc_offsets[0];
                            if (revcomp && start_offset >= doc_lens[0]) {
                                start_offset = doc_lens[0] + doc_lens[0] - start_offset - interval.second - 1;
                            }
                            break;
                        }
                    }
                    candidate_thresh[start_offset] = next_best;
                    if (!check_bwt_range(interval.first, j-1))
                        count += write_mum(interval.second, interval.first, j - 1);
                }
            start = interval.first;
            prev_lcp = prev;
        }

        if (lcp > current_mems.back().first.second) {
            if (lcp >= min_mem_length)
                current_mems.push_back(std::make_pair(std::make_pair(start, lcp), prev_lcp));
        }
        return count;
    }

    inline void update_buffers(size_t j, uint8_t bwt_c, size_t sa_pos, size_t lcp, size_t docid) {
        if (current_mems.size() <= 1) { // empty stack, only the null interval. < 1 should never happen, but included nonetheless
            if (sa_buffer.size() > 0) {
                sa_buffer.clear();
                bwt_buffer.clear();
                da_buffer.clear();
            }
            buffer_start = j;
        }
        else if (current_mems[1].first.first > buffer_start) {
            size_t to_remove = current_mems[1].first.first - buffer_start;
            buffer_start = current_mems[1].first.first;
            sa_buffer.erase(sa_buffer.begin(), sa_buffer.begin() + to_remove);
            bwt_buffer.erase(bwt_buffer.begin(), bwt_buffer.begin() + to_remove);
            da_buffer.erase(da_buffer.begin(), da_buffer.begin() + to_remove);
        }
        bwt_buffer.push_back(bwt_c);
        sa_buffer.push_back(sa_pos);
        da_buffer.push_back(docid);
    }

    void init_stack() override
    {
        current_mems.push_back(std::make_pair(std::make_pair(0, 0), 0));
    }
};

#endif /* end of include guard: _MEM_HH */
