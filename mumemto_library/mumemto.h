/* mumemto.h — C API for libmumemto (stable ABI for C / FFI).
 *
 * Link: -lmumemto (and use -Wl,-rpath,'$ORIGIN' or equivalent when loading from a non-system path).
 *
 * Example compile (after install to $PREFIX):
 *   cc -std=c11 -I"$PREFIX/include" test.c -L"$PREFIX/lib" -Wl,-rpath,"$PREFIX/lib" -lmumemto -o test
 */

#ifndef MUMEMTO_H
#define MUMEMTO_H

#include <stddef.h>
#include <stdint.h>

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
    const int64_t* offsets; /* length == num_docs for that result */
    const uint8_t* strands; /* 1 => '+', 0 => '-' */
} mumemto_mum_match_view;

typedef struct mumemto_mem_match_view {
    uint32_t length;
    size_t occurrences;
    const int64_t* offsets;
    const size_t* seq_ids;
    const uint8_t* strands;
} mumemto_mem_match_view;

/* Thread-local error string after a non-zero return from mumemto_mum / mumemto_mem. */
MUMEMTO_C_API const char* mumemto_last_error(void);

/* Build and run MUM search. Returns 0 on success and sets *out_result (caller must mum_free). */
MUMEMTO_C_API int mumemto_mum(
    const mumemto_doc_view* docs,
    size_t num_docs,
    uint32_t min_match_len,
    uint8_t use_revcomp,
    size_t num_distinct,
    uint8_t use_gsacak,
    mumemto_mum_result** out_result);

/* MEM search. Returns 0 on success; caller must mem_free. */
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

/* MUM result accessors */
MUMEMTO_C_API size_t num_docs(const mumemto_mum_result* r);
MUMEMTO_C_API const size_t* doc_record_offsets(const mumemto_mum_result* r); /* num_docs(r)+1 */
MUMEMTO_C_API const size_t* record_lengths(const mumemto_mum_result* r);
MUMEMTO_C_API size_t num_mums(const mumemto_mum_result* r);
MUMEMTO_C_API mumemto_mum_match_view mum_at(const mumemto_mum_result* r, size_t idx);
MUMEMTO_C_API void mum_free(mumemto_mum_result* r);

/* MEM result accessors */
MUMEMTO_C_API size_t num_docs_mem(const mumemto_mem_result* r);
MUMEMTO_C_API const size_t* doc_record_offsets_mem(const mumemto_mem_result* r);
MUMEMTO_C_API const size_t* record_lengths_mem(const mumemto_mem_result* r);
MUMEMTO_C_API size_t num_mems(const mumemto_mem_result* r);
MUMEMTO_C_API mumemto_mem_match_view mem_at(const mumemto_mem_result* r, size_t idx);
MUMEMTO_C_API void mem_free(mumemto_mem_result* r);

#ifdef __cplusplus
}
#endif

#endif /* MUMEMTO_H */
