#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

#include <mumsio.hpp>

// The core library is built with -fvisibility=hidden (see CMakeModules/*),
// so we must explicitly export the native C++ API symbols for consumers
// (C++ callers and the pybind11 extension).
#if defined(_WIN32)
  #if defined(MUMEMTO_BUILD_DLL)
    #define MUMEMTO_CPP_API __declspec(dllexport)
  #else
    #define MUMEMTO_CPP_API __declspec(dllimport)
  #endif
#else
  #if defined(__GNUC__) || defined(__clang__)
    #define MUMEMTO_CPP_API __attribute__((visibility("default")))
  #else
    #define MUMEMTO_CPP_API
  #endif
#endif

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

MUMEMTO_CPP_API MumResult mumemto_mum(
    const std::vector<std::vector<std::string>>& sequences,
    std::uint32_t min_match_len = 20,
    bool use_revcomp = true,
    size_t num_distinct = 0,
    bool use_gsacak = false);

MUMEMTO_CPP_API MemResult mumemto_mem(
    const std::vector<std::vector<std::string>>& sequences,
    std::uint32_t min_match_len = 20,
    bool use_revcomp = true,
    size_t num_distinct = 0,
    size_t max_total_freq = 0,
    size_t max_doc_freq = 2,
    bool use_gsacak = false);

} // namespace mumemto

/* C ABI: pure C declaration in mumemto.h. */
#include "mumemto.h"

