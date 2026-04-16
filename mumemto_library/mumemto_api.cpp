#include "mumemto_api.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <deque>
#include <new>
#include <queue>
#include <stdexcept>
#include <unordered_map>

#include <sdsl/bit_vectors.hpp>

#include <pfp.hpp>

#include <ref_builder.hpp>
#include <pfp_lcp_mum.hpp>
#include <direct_gsacak.hpp>
#include <mem_finder.hpp>

namespace mumemto {

MumResult mumemto_mum(
    const std::vector<std::vector<std::string>>& sequences,
    std::uint32_t l,
    bool r,
    std::int32_t k,
    bool g);

MemResult mumemto_mem(
    const std::vector<std::vector<std::string>>& sequences,
    std::uint32_t l,
    bool r,
    std::int32_t k,
    std::int32_t F,
    std::int32_t f,
    bool g);


class mem_finder_library{
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
        std::vector<mumsio::Mum> mums;
        std::vector<mumsio::Mem> mems;
    
        mem_finder_library(
            RefBuilder& ref_build,
            bool use_revcomp,
            size_t min_mem_len,
            size_t num_distinct_,
            int max_doc_freq_,
            int max_total_freq_)
            : min_mem_length(min_mem_len),
              num_docs(ref_build.num_docs),
              doc_lens(ref_build.seq_lengths),
              num_distinct(num_distinct_),
              doc_offsets(num_docs, 0),
              max_freq(max_total_freq_),
              max_doc_freq(max_doc_freq_),
              revcomp(use_revcomp)
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
    
            mummode = (max_doc_freq == 1);
            // Initialize stack
            init_stack();
            
            // Set parameters and limits
            // this->max_freq = num_docs + max_freq;
            this->no_max_freq = max_freq == 0;
            
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
            mumsio::Mem mem;
            mem.length = static_cast<uint32_t>(length);
            const size_t occurrences = end - start + 1;
            mem.offsets.reserve(occurrences);
            mem.seq_ids.reserve(occurrences);
            mem.strands.reserve(occurrences);

            for (size_t i = start; i <= end; i++) {
                const size_t curdoc = da_buffer.at(i - buffer_start);
                size_t curpos = sa_buffer.at(i - buffer_start) - doc_offsets[curdoc];

                bool curstrand = true; // '+' in the original output
                if (revcomp && curpos >= doc_lens[curdoc]) {
                    curstrand = false; // '-' strand
                    if (i == end)
                        curpos = doc_lens[curdoc] + doc_lens[curdoc] - curpos - length;
                    else
                        curpos = doc_lens[curdoc] + doc_lens[curdoc] - curpos - length - 1;
                }

                mem.offsets.push_back(static_cast<int64_t>(curpos));
                mem.seq_ids.push_back(curdoc);
                mem.strands.push_back(curstrand ? 1u : 0u);
            }

            mems.push_back(std::move(mem));
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
    
        // stores the LCP value preceding the current MEM interval
        size_t prev_lcp = 0;
    
        inline size_t update_mems(size_t j, size_t lcp)
        {
            // three cases for LCP, increase, decrease, or stagnant (nothing changes)
            // j = idx in SA
            size_t count = 0;
            size_t start = j - 1;
            size_t prev = 0;
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
            std::vector<int64_t> offsets(num_docs, -1);
            std::vector<char> strand(num_docs, 0);
            size_t curpos;
            size_t curdoc;
            char curstrand;

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

                offsets[curdoc] = static_cast<int64_t>(curpos);
                strand[curdoc] = curstrand;
            }

            // temporarory fix: only write MUMs where 1st genome is + stranded
            int i = 0;
            while (i < static_cast<int>(num_docs) - 1)
            {
                if (strand[i] != 0)
                    break;
                i++;
            }
            if (strand[i] == '-')
                return 0;

            mumsio::Mum m;
            m.length = static_cast<uint32_t>(length);
            m.offsets = std::move(offsets);
            m.strands.assign(num_docs, 0u);
            for (size_t doc = 0; doc < num_docs; doc++) {
                m.strands[doc] = (strand[doc] == '+') ? 1u : 0u;
            }
            mums.push_back(std::move(m));
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
    
        inline void init_stack()
        {
            current_mems.push_back(std::make_pair(std::make_pair(0, 0), 0));
        }
    };

    static inline std::vector<std::vector<size_t>> compute_record_lengths(
        const std::vector<std::vector<std::string>>& sequences) {
        std::vector<std::vector<size_t>> lengths;
        lengths.reserve(sequences.size());
        for (const auto& doc : sequences) {
            std::vector<size_t> doc_lens;
            doc_lens.reserve(doc.size());
            for (const auto& rec : doc) {
                doc_lens.push_back(rec.size());
            }
            lengths.push_back(std::move(doc_lens));
        }
        return lengths;
    }

    // --- Public library API implementations ---

    MumResult mumemto_mum(
        const std::vector<std::vector<std::string>>& sequences,
        std::uint32_t min_match_len,
        bool use_revcomp,
        size_t num_distinct,
        bool use_gsacak) {
        if (sequences.empty()) {
            return MumResult{};
        }

        const size_t num_docs = sequences.size();
        
        if (num_distinct == 0) {
            num_distinct = num_docs;
        }

        auto lengths = compute_record_lengths(sequences);
        RefBuilder ref_build = RefBuilder(lengths, use_revcomp);
        ref_build.build_input_file_lib(sequences, use_gsacak);

        // MUM mode: at most 1 occurrence per doc.
        mem_finder_library match_finder(ref_build, use_revcomp, min_match_len, num_distinct, 1, 0);
        if (use_gsacak) {
            mumemto_set_progress_enabled(false);
            gsacak_lcp gsacak("", &ref_build, false);
            gsacak.process(match_finder);
            mumemto_set_progress_enabled(true);
        } else {
            mumemto_set_progress_enabled(false);
            pf_parsing pf;
            pf.build_from_ref_builder(ref_build, 10);
            pfp_lcp lcp(pf, "", &ref_build, false);
            lcp.process(match_finder);
            mumemto_set_progress_enabled(true);
        }
        return MumResult{std::move(match_finder.mums), std::move(lengths)};
    }

    MemResult mumemto_mem(
        const std::vector<std::vector<std::string>>& sequences,
        std::uint32_t min_match_len,
        bool use_revcomp,
        size_t num_distinct,
        size_t max_total_freq,
        size_t max_doc_freq,
        bool use_gsacak) {
        if (sequences.empty()) {
            return MemResult{};
        }
        if (max_doc_freq <= 1) {
            throw std::invalid_argument("per-sequence MEM frequency f must be > 1 (use mumemto_mum instead)");
        }

        const size_t num_docs = sequences.size();

        if (num_distinct == 0) {
            num_distinct = num_docs;
        }

        auto lengths = compute_record_lengths(sequences);
        RefBuilder ref_build = RefBuilder(lengths, use_revcomp);
        ref_build.build_input_file_lib(sequences, use_gsacak);

        mem_finder_library match_finder(ref_build, use_revcomp, min_match_len, num_distinct, max_doc_freq, max_total_freq);
        if (use_gsacak) {
            mumemto_set_progress_enabled(false);
            gsacak_lcp gsacak("", &ref_build, false);
            gsacak.process(match_finder);
            mumemto_set_progress_enabled(true);
        } else {
            mumemto_set_progress_enabled(false);
            pf_parsing pf;
            pf.build_from_ref_builder(ref_build, 10);
            pfp_lcp lcp(pf, "", &ref_build, false);
            lcp.process(match_finder);
            mumemto_set_progress_enabled(true);
        }

        return MemResult{std::move(match_finder.mems), std::move(lengths)};
    }

} // namespace mumemto

// -----------------------------------------------------------------------------
// C ABI implementation
// -----------------------------------------------------------------------------

// These must be defined in the global namespace to match the forward
// declarations in the public C header (avoid anonymous-namespace type mismatch).

struct mumemto_mum_result {
    size_t num_docs = 0;
    std::vector<size_t> doc_record_offsets; // size num_docs+1
    std::vector<size_t> record_lengths;     // flattened

    std::vector<mumsio::Mum> mums;            // native result storage
};

struct mumemto_mem_result {
    size_t num_docs = 0;
    std::vector<size_t> doc_record_offsets; // size num_docs+1
    std::vector<size_t> record_lengths;     // flattened

    std::vector<mumsio::Mem> mems; // native result storage
};

namespace {

thread_local std::string g_last_error;

static inline void set_last_error(const char* msg) {
    g_last_error = (msg ? msg : "unknown error");
}

static inline void set_last_error(const std::string& msg) {
    g_last_error = msg;
}

static inline std::vector<std::vector<std::string>> build_sequences_from_docs(
    const mumemto_doc_view* docs,
    size_t num_docs) {
    std::vector<std::vector<std::string>> sequences;
    sequences.reserve(num_docs);
    for (size_t d = 0; d < num_docs; ++d) {
        std::vector<std::string> recs;
        recs.reserve(docs[d].num_records);
        for (size_t r = 0; r < docs[d].num_records; ++r) {
            const char* s = docs[d].records ? docs[d].records[r] : nullptr;
            recs.emplace_back(s ? s : "");
        }
        sequences.emplace_back(std::move(recs));
    }
    return sequences;
}

static inline void flatten_lengths(
    const std::vector<std::vector<size_t>>& lengths,
    std::vector<size_t>& out_doc_offsets,
    std::vector<size_t>& out_flat) {
    out_doc_offsets.clear();
    out_flat.clear();

    out_doc_offsets.reserve(lengths.size() + 1);
    out_doc_offsets.push_back(0);
    size_t total = 0;
    for (const auto& doc : lengths) {
        total += doc.size();
        out_doc_offsets.push_back(total);
    }
    out_flat.reserve(total);
    for (const auto& doc : lengths) {
        out_flat.insert(out_flat.end(), doc.begin(), doc.end());
    }
}

} // namespace

extern "C" {

const char* mumemto_last_error(void) {
    return g_last_error.c_str();
}

int mumemto_mum(
    const mumemto_doc_view* docs,
    size_t num_docs,
    uint32_t min_match_len,
    uint8_t use_revcomp,
    size_t num_distinct,
    uint8_t use_gsacak,
    mumemto_mum_result** out_result) {
    if (!out_result) {
        set_last_error("out_result must be non-null");
        return 1;
    }
    *out_result = nullptr;
    try {
        if (!docs && num_docs != 0) {
            set_last_error("docs must be non-null when num_docs != 0");
            return 2;
        }

        auto sequences = build_sequences_from_docs(docs, num_docs);
        auto cpp = mumemto::mumemto_mum(
            sequences,
            min_match_len,
            use_revcomp != 0,
            num_distinct,
            use_gsacak != 0);

        auto* r = new mumemto_mum_result();
        r->num_docs = sequences.size();
        flatten_lengths(cpp.lengths, r->doc_record_offsets, r->record_lengths);

        r->mums = std::move(cpp.matches);

        *out_result = r;
        return 0;
    } catch (const std::exception& e) {
        set_last_error(e.what());
        return 3;
    } catch (...) {
        set_last_error("unknown exception");
        return 4;
    }
}

int mumemto_mem(
    const mumemto_doc_view* docs,
    size_t num_docs,
    uint32_t min_match_len,
    uint8_t use_revcomp,
    size_t num_distinct,
    size_t max_total_freq,
    size_t max_doc_freq,
    uint8_t use_gsacak,
    mumemto_mem_result** out_result) {
    if (!out_result) {
        set_last_error("out_result must be non-null");
        return 1;
    }
    *out_result = nullptr;
    try {
        if (!docs && num_docs != 0) {
            set_last_error("docs must be non-null when num_docs != 0");
            return 2;
        }

        auto sequences = build_sequences_from_docs(docs, num_docs);
        auto cpp = mumemto::mumemto_mem(
            sequences,
            min_match_len,
            use_revcomp != 0,
            num_distinct,
            max_total_freq,
            max_doc_freq,
            use_gsacak != 0);

        auto* r = new mumemto_mem_result();
        r->num_docs = sequences.size();
        flatten_lengths(cpp.lengths, r->doc_record_offsets, r->record_lengths);

        r->mems = std::move(cpp.matches);

        *out_result = r;
        return 0;
    } catch (const std::exception& e) {
        set_last_error(e.what());
        return 3;
    } catch (...) {
        set_last_error("unknown exception");
        return 4;
    }
}

size_t num_docs(const mumemto_mum_result* r) {
    return r ? r->num_docs : 0;
}
const size_t* doc_record_offsets(const mumemto_mum_result* r) {
    return r ? r->doc_record_offsets.data() : nullptr;
}
const size_t* record_lengths(const mumemto_mum_result* r) {
    return r ? r->record_lengths.data() : nullptr;
}

size_t num_docs_mem(const mumemto_mem_result* r) {
    return r ? r->num_docs : 0;
}
const size_t* doc_record_offsets_mem(const mumemto_mem_result* r) {
    return r ? r->doc_record_offsets.data() : nullptr;
}
const size_t* record_lengths_mem(const mumemto_mem_result* r) {
    return r ? r->record_lengths.data() : nullptr;
}

size_t num_mums(const mumemto_mum_result* r) {
    return r ? r->mums.size() : 0;
}

mumemto_mum_match_view mum_at(const mumemto_mum_result* r, size_t idx) {
    mumemto_mum_match_view v{};
    if (!r) return v;
    if (idx >= r->mums.size()) return v;
    v.length = r->mums[idx].length;
    v.offsets = r->mums[idx].offsets.data();
    v.strands = r->mums[idx].strands.data();
    return v;
}

size_t num_mems(const mumemto_mem_result* r) {
    return r ? r->mems.size() : 0;
}

mumemto_mem_match_view mem_at(const mumemto_mem_result* r, size_t idx) {
    mumemto_mem_match_view v{};
    if (!r) return v;
    if (idx >= r->mems.size()) return v;
    v.length = r->mems[idx].length;
    v.occurrences = r->mems[idx].offsets.size();
    v.offsets = r->mems[idx].offsets.data();
    v.seq_ids = r->mems[idx].seq_ids.data();
    v.strands = r->mems[idx].strands.data();
    return v;
}

void mum_free(mumemto_mum_result* r) {
    delete r;
}
void mem_free(mumemto_mem_result* r) {
    delete r;
}

} // extern "C"