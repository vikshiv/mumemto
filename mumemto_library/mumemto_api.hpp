#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

#include <mumsio.hpp>

namespace mumemto {

struct MumResult {
    std::vector<mumsio::Mum> matches;
    // Per input "FASTA record" lengths (matches RefBuilder.multifasta_lengths),
    // stored in the same nested structure: lengths[doc][record].
    std::vector<std::vector<size_t>> lengths;
};

struct MemResult {
    std::vector<mumsio::Mem> matches;
    // Per input "FASTA record" lengths (matches RefBuilder.multifasta_lengths),
    // stored in the same nested structure: lengths[doc][record].
    std::vector<std::vector<size_t>> lengths;
};

MumResult mumemto_mum(
    std::vector<std::vector<std::string>>& sequences,
    std::uint32_t min_match_len = 20,
    bool use_revcomp = true,
    size_t num_distinct = 0,
    bool use_gsacak = false);

MemResult mumemto_mem(
    std::vector<std::vector<std::string>>& sequences,
    std::uint32_t min_match_len = 20,
    bool use_revcomp = true,
    size_t num_distinct = 0,
    size_t max_total_freq = 0,
    size_t max_doc_freq = 1,
    bool use_gsacak = false);

} // namespace mumemto

/* C ABI: pure C declaration in mumemto.h. */
#include "mumemto.h"

