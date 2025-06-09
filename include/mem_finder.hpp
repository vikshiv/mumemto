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
#include <filesystem>

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
    bool binary;
    bool merge;
    bool anchor_merge;
    std::string filename;

    mem_finder(std::string filename, RefBuilder& ref_build, size_t min_mem_len, size_t num_distinct, int max_doc_freq, int max_total_freq, bool binary, bool merge, bool anchor_merge): 
                min_mem_length(min_mem_len),
                num_docs(ref_build.num_docs),
                revcomp(ref_build.use_revcomp),
                doc_lens(ref_build.seq_lengths),
                num_distinct(num_distinct),
                doc_offsets(num_docs, 0),
                max_freq(max_total_freq),
                max_doc_freq(max_doc_freq),
                binary(binary),
                filename(filename),
                merge(merge),
                candidate_thresh(0), 
                anchor_merge(anchor_merge)
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
        if (merge) {
            candidate_thresh.resize(doc_lens[0] * 2, 0);
        }

        this->mummode = (max_doc_freq == 1);
        // Initialize stack
        init_stack();
        
        // Set parameters and limits
        // this->max_freq = num_docs + max_freq;
        this->no_max_freq = max_freq == 0;
        // Opening output file
        if (!mummode)
            mem_file.open(filename + std::string(".mems"));
        else if (binary)
        {
            bums_lengths.open(filename + std::string(".bumbl.lengths"), std::ios::binary);
            bums_starts.open(filename + std::string(".bumbl.starts"), std::ios::binary);
            bums_strands.open(filename + std::string(".bumbl.strands"), std::ios::binary);
        }
            
        else
            mem_file.open(filename + std::string(".mums"));
    }

    virtual void close()
    {
        if (binary)
            write_bums();
        else
            mem_file.close();
        if (anchor_merge) {
            std::string thresh_file = filename + ".athresh";
            std::ofstream thresh_out(thresh_file, std::ios::binary);
            thresh_out.write(reinterpret_cast<const char*>(candidate_thresh.data()), doc_lens[0] * sizeof(uint16_t));
            thresh_out.close();
        }
        else if (merge) {
            // write merging metadata
            // write thresholds in terms of positions on MUMs, rather than first genome coords
            size_t total_mum_length = 0;
            for (size_t i = 0; i < mum_positions.size(); i++) {
                total_mum_length += mum_positions[i].second + 1;
            }
            size_t offset = 0;
            std::vector<uint16_t> mum_based_thresh(total_mum_length, 0);
            std::vector<uint16_t> mum_based_thresh_rev(total_mum_length, 0);
            size_t revpos;
            for (size_t i = 0; i < mum_positions.size(); i++) {
                // revpos = doc_lens[0] - mum_positions[i].first - mum_positions[i].second - 1;
                // curpos = doc_lens[curdoc] + doc_lens[curdoc] - curpos - length - 1;
                // 
                revpos = doc_lens[0] + doc_lens[0] - (mum_positions[i].first) - mum_positions[i].second - 1;
                for (size_t j = 0; j < mum_positions[i].second; j++) {
                    if (candidate_thresh[mum_positions[i].first + j] < mum_positions[i].second - j)
                        mum_based_thresh[offset] = candidate_thresh[mum_positions[i].first + j];
                    if (candidate_thresh.at(revpos + j) < mum_positions[i].second - j)
                        mum_based_thresh_rev[offset] = candidate_thresh.at(revpos + j);
                    offset++;
                }
                mum_based_thresh[offset] = 0;
                mum_based_thresh_rev[offset] = 0;
                offset++;
            }
            // write to file
            std::string thresh_file = filename + ".thresh";
            std::ofstream thresh_out(thresh_file, std::ios::binary);
            thresh_out.write(reinterpret_cast<const char*>(mum_based_thresh.data()), mum_based_thresh.size() * sizeof(uint16_t));
            thresh_out.close();

            std::string thresh_file_rev = filename + ".thresh_rev";
            std::ofstream thresh_out_rev(thresh_file_rev, std::ios::binary);
            thresh_out_rev.write(reinterpret_cast<const char*>(mum_based_thresh_rev.data()), mum_based_thresh_rev.size() * sizeof(uint16_t));
            thresh_out_rev.close();
        }
    }

    // main update function, takes in the current streamed value of each array and write mem if found
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

protected:
    std::ofstream mem_file;
    std::ofstream bums_lengths;
    std::ofstream bums_starts;
    std::ofstream bums_strands;
    std::vector<std::vector<char>> bums_strands_vec;
    
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
                curpos = doc_lens[curdoc] + doc_lens[curdoc] - curpos - length - 1;
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
    // current stack of MEMs, ((start idx in SA, length of mem), prev_lcp)
    std::vector<std::pair<std::pair<size_t, size_t>, size_t>> current_mems; 
    // positions of MUMs in the first genome (offset, length), for merging
    std::vector<std::pair<size_t, size_t>> mum_positions;
    // data structure to hold meta data for merging
    std::vector<uint16_t> candidate_thresh;
    size_t MAX_THRESH = static_cast<size_t>(UINT16_MAX);

    // stores the LCP value preceding the current MEM interval
    size_t prev_lcp = 0;

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
                    if (merge) {
                        // update candidate threshold
                        next_best = std::min(std::max(prev, lcp), MAX_THRESH); // cap at max uint16
                        for (size_t i = interval.first; i <= j-1; i++) {
                            if (da_buffer[i - buffer_start] == 0) {
                                start_offset = sa_buffer[i - buffer_start] - doc_offsets[0];
                                candidate_thresh[start_offset] = next_best;
                                break;
                            }
                        }
                    }
                    
                    if (!check_bwt_range(interval.first, j-1)) {
                        if (mummode)
                            count += write_mum(interval.second, interval.first, j - 1);
                        else
                            count += write_mem(interval.second, interval.first, j - 1);
                    }
                    
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
                if (curpos + length >= (doc_lens[curdoc] + doc_lens[curdoc]))
                    return 0;
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

        // store the offset in the first genome to pull later
        if (merge) {
            mum_positions.push_back(std::make_pair(offsets[0], length));
        }
        
        if (binary) {
            uint16_t length_uint16 = static_cast<uint16_t>(length);
            bums_lengths.write(reinterpret_cast<const char*>(&length_uint16), sizeof(length_uint16));
            std::vector<uint64_t> converted(offsets.begin(), offsets.end());
            bums_starts.write(reinterpret_cast<const char*>(converted.data()), sizeof(uint64_t) * offsets.size());
            bums_strands_vec.push_back(strand);
        }
        else {
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
        }
        return 1;
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

    inline uint16_t get_flags() {
        uint16_t flags = 0;
        if (num_distinct < num_docs)
            flags |= 1 << 0; 
        return flags;
    }

    inline void write_bums() {
        size_t num_mems = bums_strands_vec.size();
        size_t num_seqs = bums_strands_vec[0].size();
        // Prepare packed data buffer
        std::vector<uint8_t> buffer((num_mems * num_seqs + 7) / 8, 0);
        // Pack bits into bytes
        for (size_t i = 0; i < num_mems; ++i) {
            for (size_t j = 0; j < num_seqs; ++j) {
                size_t bitIndex = i * num_seqs + j;
                if (bums_strands_vec[i][j] == '+') {
                    buffer[bitIndex / 8] |= (1 << (7 - (bitIndex % 8))); // Set bit for '+'
                }
            }
        }
        // Write packed data to file
        bums_strands.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
        bums_strands.close();
        bums_lengths.close();
        bums_starts.close();
        mem_file.open(filename + std::string(".bumbl"), std::ios::binary);
        // Combine intermediate files into a single output file
        if (mem_file.is_open()) {
            uint16_t flags = get_flags();
            mem_file.write(reinterpret_cast<const char*>(&flags), sizeof(flags));
            mem_file.write(reinterpret_cast<const char*>(&num_seqs), sizeof(uint64_t));
            mem_file.write(reinterpret_cast<const char*>(&num_mems), sizeof(uint64_t));

            // Read and write the contents of the intermediate files
            std::ifstream lengths_file(filename + std::string(".bumbl.lengths"), std::ios::binary);
            std::ifstream starts_file(filename + std::string(".bumbl.starts"), std::ios::binary);
            std::ifstream strands_file(filename + std::string(".bumbl.strands"), std::ios::binary);

            mem_file << lengths_file.rdbuf();
            mem_file << starts_file.rdbuf();
            mem_file << strands_file.rdbuf();

            mem_file.close();
        }

        // Optionally, remove the intermediate files
        std::filesystem::remove(filename + std::string(".bumbl.lengths"));
        std::filesystem::remove(filename + std::string(".bumbl.starts"));
        std::filesystem::remove(filename + std::string(".bumbl.strands"));
    }

    virtual void init_stack()
    {
        current_mems.push_back(std::make_pair(std::make_pair(0, 0), 0));
    }
};

#endif /* end of include guard: _MEM_HH */
