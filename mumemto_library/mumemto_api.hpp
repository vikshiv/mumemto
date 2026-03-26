#pragma once

#include <cstdint>
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
    const std::vector<std::vector<std::string>>& sequences,
    std::uint32_t l = 20,
    bool r = true,
    std::int32_t k = 0,
    bool g = false);

MemResult mumemto_mem(
    const std::vector<std::vector<std::string>>& sequences,
    std::uint32_t l = 20,
    bool r = true,
    std::int32_t k = 0,
    std::int32_t F = 0,
    std::int32_t f = 1,
    bool g = false);

} // namespace mumemto
