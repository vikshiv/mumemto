# Mumemto library quick start (Python / C / C++)

Mumemto ships a shared library (`libmumemto.so`), headers, and Python bindings. This page is a **quick-start** for using the library APIs (not the CLI).

If you only want to run Mumemto from the command line, see the main `[README.md](../README.md)` and the [Mumemto wiki](https://github.com/vikshiv/mumemto/wiki).

---

## Installation

### Conda installation (recommended)

Conda installs:

- the shared library + headers (`mumemto.h`, `mumemto_api.hpp`, `mumemto.hpp`)
- the CMake package (`find_package(Mumemto)`)
- the Python package + extension module

```sh
conda create -n mumemto_env python
conda activate mumemto_env
conda install -c conda-forge bioconda::mumemto
```

**Linking C/C++ against conda `mumemto`:** install **`cxx-compiler`** in the same env (`mamba install cxx-compiler`) and **compile with conda’s toolchain** (`c++` / `cc` from the activated environment), not the system `g++`.

### Build from source (local install tree)

This installs the library, headers, and a CMake config into `build/install/` by default:

```sh
git clone https://github.com/vikshiv/mumemto
cd mumemto
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cmake --install build
```

> [!TIP]
> If you want to consume the library with `find_package(Mumemto)`, point `CMAKE_PREFIX_PATH` at the install root (e.g. `.../mumemto/build/install`).

---

## Python quick start

After installing via conda/pip:

```python
import mumemto

# sequences[doc][record]
res = mumemto.mem([["ACGTACGT"], ["TACGTAAA"]], min_match_len=3, max_doc_freq=2)

print(res.num_docs(), res.num_matches())
length, offsets, seq_ids, strands = res.match_at(0)

mum_res = mumemto.mum([["ACGTACGT"], ["TACGTAAA"]], min_match_len=3)
length, offsets, strands = mum_res.match_at(0)

first = res[0]  # same tuple as match_at(0)
for length, offsets, seq_ids, strands in len(res):
    ...
for length, offsets, strands in mum_res:
    ...
```

> [!NOTE]
> The arrays returned by `match_at()` are **views** into result storage. If you want to modify them, make a copy (e.g. `offsets = offsets.copy()`).

---

## Which API should I use?

### Python bindings

Use this if you want to call Mumemto from Python. The Python API is designed around returning **NumPy views** into result buffers, and Python inputs are `list[list[str]]`.

### C ABI (`mumemto.h`) — stable ABI for C / FFI

Use this if you want:

- C (or another language via FFI)
- a small stable surface (`extern "C"`)

You pass a list of assemblies, which can each contain multiple sequences represented as  pointers to NUL-terminated C strings, then query results through accessor functions. Results are freed via `mum_free` / `mem_free`.

### Native C++ API (`mumemto_api.hpp`) — best for C++ callers

Use this if you want:

- direct C++ calls and C++ data structures
- no C-ABI translation layer

This is the recommended interface for C++ programs that link against `libmumemto`.

### C++ convenience wrapper (`mumemto.hpp`) — RAII over the C ABI

Use this if you want:

- a header-only convenience layer that returns lightweight “view” structs
- RAII ownership of the C ABI result handles

This wraps the C ABI and is convenient for quick prototypes and tools.

---

## C quick start (C ABI)

Minimal example (in-memory strings):

```c
#include <stdio.h>
#include "mumemto.h"

int main(void) {
  const char* doc0[] = {"ACGTACGT"};
  const char* doc1[] = {"TACGTAAA"};
  mumemto_doc_view docs[2];
  docs[0].records = doc0; docs[0].num_records = 1;
  docs[1].records = doc1; docs[1].num_records = 1;

  mumemto_mum_result* res = NULL;
  if (mumemto_mum(docs, 2, /*min_match_len=*/3, /*use_revcomp=*/1u,
                 /*num_distinct=*/0, /*use_gsacak=*/0u, &res) != 0) {
    fprintf(stderr, "mumemto_mum failed: %s\n", mumemto_last_error());
    return 1;
  }

  for (size_t i = 0; i < num_mums(res); ++i) {
    mumemto_mum_match_view m = mum_at(res, i);
    printf("len=%u off0=%lld off1=%lld\n",
           (unsigned)m.length, (long long)m.offsets[0], (long long)m.offsets[1]);
  }

  mum_free(res);
  return 0;
}
```

Compile/link (if installed into `$PREFIX`, Linux):

```sh
cc -std=c11 -I"$PREFIX/include" example.c -L"$PREFIX/lib" -Wl,-rpath,"$PREFIX/lib" -lmumemto -lstdc++ -o example
./example
```

(`libmumemto` is C++; with `cc`, add `-lstdc++`. With conda, use the env’s `cc` / `c++` as above.)

---

## C++ quick start (native C++ API)

Minimal example:

```cpp
#include <string>
#include <vector>
#include "mumemto_api.hpp"

int main() {
  std::vector<std::vector<std::string>> sequences{
      {"ACGTACGT"},
      {"TACGTAAA"},
  };

  auto r = mumemto::mumemto_mem(
      sequences,
      /*min_match_len=*/3,
      /*use_revcomp=*/true,
      /*num_distinct=*/0,
      /*max_total_freq=*/0,
      /*max_doc_freq=*/2,
      /*use_gsacak=*/false);

  // r.matches is a vector of mumsio::Mem
  return (int)r.matches.size();
}
```

### C++ convenience wrapper (`mumemto.hpp`)

The same inputs as above, using **`mumemto_cxx::mum`** / **`mumemto_cxx::mem`**: RAII over the C ABI (no manual memory management needed).

```cpp
#include <vector>
#include "mumemto.hpp"

int main() {
  std::vector<std::vector<std::string>> sequences{
      {"ACGTACGT"},
      {"TACGTAAA"},
  };

  auto r = mumemto_cxx::mem(
      sequences,
      /*min_match_len=*/3,
      /*use_revcomp=*/true,
      /*num_distinct=*/0,
      /*max_total_freq=*/0,
      /*max_doc_freq=*/2,
      /*use_gsacak=*/false);

  return static_cast<int>(r.num_matches());
}
```

Use **`mumemto_cxx::mum`** for MUM mode.

### CMake integration (`find_package(Mumemto)`)

```cmake
find_package(Mumemto CONFIG REQUIRED)
add_executable(example example.cpp)
target_link_libraries(example PRIVATE Mumemto::mumemto)
target_compile_features(example PRIVATE cxx_std_17)
```

---

## AI disclosure

Portions of the C ABI and python bindings, along with this README, were written with using AI; the maintainers have reviewed the technical content and run manually (non-AI-written) tests for correctness.

