#pragma once

// Header-only C++ convenience wrapper around the C API (mumemto_api.hpp).
// Provides RAII so C++ callers don't manually call *_free.

#include "mumemto_api.hpp"

#include <cstdint>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace mumemto_cxx {

struct MumResultDeleter {
    void operator()(mumemto_mum_result* p) const noexcept { mum_free(p); }
};
struct MemResultDeleter {
    void operator()(mumemto_mem_result* p) const noexcept { mem_free(p); }
};

class MumResult {
public:
    explicit MumResult(mumemto_mum_result* p) : handle_(p) {}
    MumResult(const MumResult&) = delete;
    MumResult& operator=(const MumResult&) = delete;
    MumResult(MumResult&&) noexcept = default;
    MumResult& operator=(MumResult&&) noexcept = default;

    size_t num_docs() const { return ::num_docs(handle_.get()); }
    size_t num_matches() const { return num_mums(handle_.get()); }

    struct Match {
        uint32_t length = 0;
        const int64_t* offsets = nullptr;   // size num_docs()
        const uint8_t* strands = nullptr;   // size num_docs(), 1 => '+', 0 => '-'
        size_t num_docs = 0;
    };

    Match match_at(size_t idx) const {
        const auto view = mum_at(handle_.get(), idx);
        Match m;
        m.length = view.length;
        m.offsets = view.offsets;
        m.strands = view.strands;
        m.num_docs = num_docs();
        return m;
    }

    // Flattened record lengths for all docs:
    //   doc_record_offsets()[doc] .. doc_record_offsets()[doc+1]-1 indexes record_lengths()
    const size_t* doc_record_offsets() const { return ::doc_record_offsets(handle_.get()); }
    const size_t* record_lengths() const { return ::record_lengths(handle_.get()); }

    // Zero-copy text serialization compatible with `mumsio::serialize_mum` output.
    void write_mums(const std::string& path) const {
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Failed to open output file: " + path);

        const size_t nd = num_docs();
        const size_t nm = num_matches();
        for (size_t i = 0; i < nm; ++i) {
            const auto v = mum_at(handle_.get(), i);
            out << v.length << "\t";
            for (size_t d = 0; d < nd; ++d) {
                out << v.offsets[d];
                if (d + 1 != nd) out << ",";
            }
            out << "\t";
            for (size_t d = 0; d < nd; ++d) {
                out << (v.strands[d] ? "+" : "-");
                if (d + 1 != nd) out << ",";
            }
            out << "\n";
        }
    }

    // Zero-copy .bumbl writer (same layout as `mumsio::write_bumbl`).
    // If `partial` is false, we auto-detect partials (-1 offsets) without allocating.
    void write_bumbl(const std::string& path, bool partial = false, bool coll_blocks = false) const {
        const uint64_t nd = static_cast<uint64_t>(num_docs());
        const uint64_t nm = static_cast<uint64_t>(num_matches());
        if (nd == 0 || nm == 0) {
            std::ofstream out(path, std::ios::binary);
            if (!out) throw std::runtime_error("Failed to open output file: " + path);
            // Write a minimal header (flags + zeros) to avoid surprises.
            uint16_t flags = 0;
            flags |= static_cast<uint16_t>(1u << 15); // length32 always set
            out.write(reinterpret_cast<const char*>(&flags), sizeof(flags));
            out.write(reinterpret_cast<const char*>(&nd), sizeof(nd));
            out.write(reinterpret_cast<const char*>(&nm), sizeof(nm));
            return;
        }

        bool any_partial = partial;
        if (!any_partial) {
            for (uint64_t i = 0; i < nm && !any_partial; ++i) {
                const auto v = mum_at(handle_.get(), static_cast<size_t>(i));
                for (uint64_t d = 0; d < nd; ++d) {
                    if (v.offsets[d] == -1) { any_partial = true; break; }
                }
            }
        }

        std::ofstream out(path, std::ios::binary);
        if (!out) throw std::runtime_error("Failed to open output file: " + path);

        uint16_t flags = 0;
        flags |= static_cast<uint16_t>(1u << 15); // length32 always set
        if (any_partial)  flags |= static_cast<uint16_t>(1u << 13);
        if (coll_blocks)  flags |= static_cast<uint16_t>(1u << 14);
        out.write(reinterpret_cast<const char*>(&flags), sizeof(flags));
        out.write(reinterpret_cast<const char*>(&nd), sizeof(nd));
        out.write(reinterpret_cast<const char*>(&nm), sizeof(nm));

        // lengths (uint32_t, nm entries)
        for (uint64_t i = 0; i < nm; ++i) {
            const auto v = mum_at(handle_.get(), static_cast<size_t>(i));
            const uint32_t len = v.length;
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        }

        // starts (int64_t, nm*nd entries)
        for (uint64_t i = 0; i < nm; ++i) {
            const auto v = mum_at(handle_.get(), static_cast<size_t>(i));
            for (uint64_t d = 0; d < nd; ++d) {
                const int64_t start = v.offsets[d];
                out.write(reinterpret_cast<const char*>(&start), sizeof(start));
            }
        }

        // packed strands (MSB-first)
        uint8_t byte = 0;
        uint8_t bit_in_byte = 0; // 0..7, where 0 corresponds to MSB
        for (uint64_t i = 0; i < nm; ++i) {
            const auto v = mum_at(handle_.get(), static_cast<size_t>(i));
            for (uint64_t d = 0; d < nd; ++d) {
                const uint8_t b = v.strands[d] ? 1u : 0u;
                if (b) {
                    byte |= static_cast<uint8_t>(1u << (7u - bit_in_byte));
                }
                bit_in_byte++;
                if (bit_in_byte == 8) {
                    out.write(reinterpret_cast<const char*>(&byte), 1);
                    byte = 0;
                    bit_in_byte = 0;
                }
            }
        }
        if (bit_in_byte != 0) {
            out.write(reinterpret_cast<const char*>(&byte), 1);
        }

        if (coll_blocks) {
            // Blocks are not produced by this library path; write 0 blocks to keep
            // the format consistent with `parse_bumbl`.
            const uint64_t num_blocks = 0;
            out.write(reinterpret_cast<const char*>(&num_blocks), sizeof(num_blocks));
        }
    }

private:
    std::unique_ptr<mumemto_mum_result, MumResultDeleter> handle_;
};

class MemResult {
public:
    explicit MemResult(mumemto_mem_result* p) : handle_(p) {}
    MemResult(const MemResult&) = delete;
    MemResult& operator=(const MemResult&) = delete;
    MemResult(MemResult&&) noexcept = default;
    MemResult& operator=(MemResult&&) noexcept = default;

    size_t num_docs() const { return num_docs_mem(handle_.get()); }
    size_t num_matches() const { return num_mems(handle_.get()); }

    struct Match {
        uint32_t length = 0;
        size_t occurrences = 0;
        const int64_t* offsets = nullptr;  // size occurrences
        const size_t* seq_ids = nullptr;   // size occurrences
        const uint8_t* strands = nullptr;  // size occurrences, 1 => '+', 0 => '-'
    };

    Match match_at(size_t idx) const {
        const auto view = mem_at(handle_.get(), idx);
        Match m;
        m.length = view.length;
        m.occurrences = view.occurrences;
        m.offsets = view.offsets;
        m.seq_ids = view.seq_ids;
        m.strands = view.strands;
        return m;
    }

    const size_t* doc_record_offsets() const { return ::doc_record_offsets_mem(handle_.get()); }
    const size_t* record_lengths() const { return ::record_lengths_mem(handle_.get()); }

    // Zero-copy text serialization compatible with `mumsio::serialize_mem` output.
    void write_mems(const std::string& path) const {
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Failed to open output file: " + path);
        const size_t nm = num_matches();
        for (size_t i = 0; i < nm; ++i) {
            const auto v = mem_at(handle_.get(), i);
            out << v.length << "\t";
            for (size_t j = 0; j < v.occurrences; ++j) {
                out << v.offsets[j];
                if (j + 1 != v.occurrences) out << ",";
            }
            out << "\t";
            for (size_t j = 0; j < v.occurrences; ++j) {
                out << v.seq_ids[j];
                if (j + 1 != v.occurrences) out << ",";
            }
            out << "\t";
            for (size_t j = 0; j < v.occurrences; ++j) {
                out << (v.strands[j] ? "+" : "-");
                if (j + 1 != v.occurrences) out << ",";
            }
            out << "\n";
        }
    }

private:
    std::unique_ptr<mumemto_mem_result, MemResultDeleter> handle_;
};

namespace detail {
inline uint8_t as_u8(bool b) { return b ? 1u : 0u; }

inline std::vector<std::vector<const char*>> build_record_ptrs(
    const std::vector<std::vector<std::string>>& sequences) {
    std::vector<std::vector<const char*>> out;
    out.reserve(sequences.size());
    for (const auto& doc : sequences) {
        std::vector<const char*> recs;
        recs.reserve(doc.size());
        for (const auto& s : doc) recs.push_back(s.c_str());
        out.push_back(std::move(recs));
    }
    return out;
}

inline std::vector<mumemto_doc_view> build_docs_from_strings(
    const std::vector<std::vector<std::string>>& sequences,
    std::vector<std::vector<const char*>>& record_ptrs) {
    record_ptrs = build_record_ptrs(sequences);
    std::vector<mumemto_doc_view> docs;
    docs.reserve(sequences.size());
    for (size_t d = 0; d < sequences.size(); ++d) {
        mumemto_doc_view dv;
        dv.records = record_ptrs[d].data();
        dv.num_records = record_ptrs[d].size();
        docs.push_back(dv);
    }
    return docs;
}
} // namespace detail

inline MumResult mum(
    const std::vector<std::vector<std::string>>& sequences,
    uint32_t min_match_len = 20,
    bool use_revcomp = true,
    size_t num_distinct = 0,
    bool use_gsacak = false) {
    std::vector<std::vector<const char*>> record_ptrs;
    auto docs = detail::build_docs_from_strings(sequences, record_ptrs);

    mumemto_mum_result* out = nullptr;
    const int rc = mumemto_mum(
        docs.data(),
        docs.size(),
        min_match_len,
        detail::as_u8(use_revcomp),
        num_distinct,
        detail::as_u8(use_gsacak),
        &out);
    if (rc != 0) throw std::runtime_error(mumemto_last_error());
    return MumResult(out);
}

inline MemResult mem(
    const std::vector<std::vector<std::string>>& sequences,
    uint32_t min_match_len = 20,
    bool use_revcomp = true,
    size_t num_distinct = 0,
    size_t max_total_freq = 0,
    size_t max_doc_freq = 1,
    bool use_gsacak = false) {
    std::vector<std::vector<const char*>> record_ptrs;
    auto docs = detail::build_docs_from_strings(sequences, record_ptrs);

    mumemto_mem_result* out = nullptr;
    const int rc = mumemto_mem(
        docs.data(),
        docs.size(),
        min_match_len,
        detail::as_u8(use_revcomp),
        num_distinct,
        max_total_freq,
        max_doc_freq,
        detail::as_u8(use_gsacak),
        &out);
    if (rc != 0) throw std::runtime_error(mumemto_last_error());
    return MemResult(out);
}

} // namespace mumemto_cxx

