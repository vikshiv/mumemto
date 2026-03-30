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

// -----------------------------------------------------------------------------
// C ABI wrappers
// -----------------------------------------------------------------------------
//
// These functions provide an `extern "C"` stable ABI for consuming `mumemto`
// from C (and other languages via FFI). The API:
// - Uses only C types in exported signatures.
// - Returns opaque result handles; results must be freed by the caller.
// - For match payloads, returned pointers reference data owned by the opaque
//   handle (no extra copying of offsets/seq_ids/strands).
//
// The input `docs` corresponds to the C++ `sequences` nested structure:
// docs[doc].records[rec] == sequences[doc][rec]

#if defined(_WIN32)
  #if defined(MUMEMTO_BUILD_DLL)
    #define MUMEMTO_C_API __declspec(dllexport)
  #else
    #define MUMEMTO_C_API __declspec(dllimport)
  #endif
#else
  #if defined(__GNUC__) || defined(__clang__)
    #define MUMEMTO_C_API __attribute__((visibility("default")))
  #else
    #define MUMEMTO_C_API
  #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mumemto_doc_view {
    const char* const* records;
    size_t num_records;
} mumemto_doc_view;

typedef struct mumemto_mum_result mumemto_mum_result;
typedef struct mumemto_mem_result mumemto_mem_result;

typedef struct mumemto_mum_match_view {
    uint32_t length;
    const int64_t* offsets;  // length == num_docs
    const uint8_t* strands;  // length == num_docs, 1 => '+', 0 => '-'
} mumemto_mum_match_view;

typedef struct mumemto_mem_match_view {
    uint32_t length;
    size_t occurrences;
    const int64_t* offsets;  // length == occurrences
    const size_t* seq_ids;   // length == occurrences
    const uint8_t* strands;  // length == occurrences, 1 => '+', 0 => '-'
} mumemto_mem_match_view;

// Error reporting: on failure, functions return non-zero and set an internal
// thread-local error string retrievable via this function.
MUMEMTO_C_API const char* mumemto_last_error(void);

// Build + run MUM/MEM search. On success returns 0 and sets *out_result.
MUMEMTO_C_API int mumemto_mum(
    const mumemto_doc_view* docs,
    size_t num_docs,
    uint32_t min_match_len,
    uint8_t use_revcomp,
    size_t num_distinct,
    uint8_t use_gsacak,
    mumemto_mum_result** out_result);

MUMEMTO_C_API int mumemto_mem(
    const mumemto_doc_view* docs,
    size_t num_docs,
    uint32_t min_match_len,
    uint8_t use_revcomp,
    size_t num_distinct,
    size_t max_total_freq,
    size_t max_doc_freq,
    uint8_t use_gsacak,
    mumemto_mem_result** out_result);

// MUM accessors
MUMEMTO_C_API size_t num_docs(const mumemto_mum_result* r);
MUMEMTO_C_API const size_t* doc_record_offsets(const mumemto_mum_result* r); // size num_docs(r)+1
MUMEMTO_C_API const size_t* record_lengths(const mumemto_mum_result* r);     // size doc_record_offsets(r)[num_docs(r)]
MUMEMTO_C_API size_t num_mums(const mumemto_mum_result* r);
MUMEMTO_C_API mumemto_mum_match_view mum_at(const mumemto_mum_result* r, size_t idx);
MUMEMTO_C_API void mum_free(mumemto_mum_result* r);

// MEM accessors
MUMEMTO_C_API size_t num_docs_mem(const mumemto_mem_result* r);
MUMEMTO_C_API const size_t* doc_record_offsets_mem(const mumemto_mem_result* r); // size num_docs_mem(r)+1
MUMEMTO_C_API const size_t* record_lengths_mem(const mumemto_mem_result* r);     // size doc_record_offsets_mem(r)[num_docs_mem(r)]
MUMEMTO_C_API size_t num_mems(const mumemto_mem_result* r);
MUMEMTO_C_API mumemto_mem_match_view mem_at(const mumemto_mem_result* r, size_t idx);
MUMEMTO_C_API void mem_free(mumemto_mem_result* r);

#ifdef __cplusplus
} // extern "C"
#endif

