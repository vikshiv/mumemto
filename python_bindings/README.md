## mumemto Python bindings

This directory builds a `pybind11` extension that binds the **C ABI** of `mumemto`
and exposes **zero-copy NumPy views** into the result payload (offsets/strands/ids).

### Build prerequisites

- A built+installed Mumemto CMake package (`find_package(Mumemto CONFIG)` must work)
- Python + pip

From the repo root, one straightforward way is:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cmake --install build
```

Then build the wheel:

```bash
python -m pip install -U pip
python -m pip install -U build
python -m pip wheel ./python_bindings -w dist \
  -Ccmake.args="-DCMAKE_PREFIX_PATH=$(pwd)/build/install"
```

Or for an editable install:

```bash
python -m pip install -e ./python_bindings \
  -Ccmake.args="-DCMAKE_PREFIX_PATH=$(pwd)/build/install"
```

### Usage

```python
import mumemto

seqs = [["ACGTACGT"], ["TACGTAAA"]]
r = mumemto.mum(seqs, min_match_len=3)
print(r.num_docs(), r.num_matches())
m0 = r.match_at(0)
print(m0["length"], m0["offsets"], m0["strands"])
```

Notes:
- `match_at()` returns a dict where `offsets` / `seq_ids` / `strands` are **NumPy arrays
  backed by the C++ result memory** (no copy). Keep the `MumResult`/`MemResult` alive.

