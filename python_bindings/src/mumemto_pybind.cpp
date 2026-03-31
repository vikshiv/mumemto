#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "mumemto_api.hpp"

namespace py = pybind11;

namespace {

struct MumHandleDeleter {
  void operator()(mumemto_mum_result* p) const noexcept { mum_free(p); }
};
struct MemHandleDeleter {
  void operator()(mumemto_mem_result* p) const noexcept { mem_free(p); }
};

template <typename T>
py::array view_1d(py::handle base, const T* ptr, size_t n) {
  if (!ptr) {
    return py::array(py::dtype::of<T>(), {static_cast<py::ssize_t>(0)});
  }
  return py::array(
      py::dtype::of<T>(),
      {static_cast<py::ssize_t>(n)},
      {static_cast<py::ssize_t>(sizeof(T))},
      ptr,
      base);
}

inline py::array view_bool_1d_from_u8(py::handle base, const uint8_t* ptr, size_t n) {
  if (!ptr) {
    return py::array(py::dtype("bool"), {static_cast<py::ssize_t>(0)});
  }
  // NumPy bool_ is stored as 1 byte; interpret the existing 0/1 bytes as bool.
  return py::array(
      py::dtype("bool"),
      {static_cast<py::ssize_t>(n)},
      {static_cast<py::ssize_t>(sizeof(uint8_t))},
      ptr,
      base);
}

inline uint8_t as_u8(bool b) { return b ? 1u : 0u; }

class MumResult {
public:
  explicit MumResult(mumemto_mum_result* h) : h_(h) {}

  size_t num_docs() const { return ::num_docs(h_.get()); }
  size_t num_matches() const { return ::num_mums(h_.get()); }

  py::array doc_record_offsets() const {
    const size_t nd = num_docs();
    return view_1d<size_t>(py::cast(this), ::doc_record_offsets(h_.get()), nd + 1);
  }

  py::array record_lengths() const {
    const size_t nd = num_docs();
    const size_t* off = ::doc_record_offsets(h_.get());
    const size_t total = (off && nd > 0) ? off[nd] : 0;
    return view_1d<size_t>(py::cast(this), ::record_lengths(h_.get()), total);
  }

  py::tuple match_at(size_t idx) const {
    const auto v = ::mum_at(h_.get(), idx);
    // (length, offsets[num_docs], strands[num_docs])
    return py::make_tuple(
        v.length,
        view_1d<int64_t>(py::cast(this), v.offsets, num_docs()),
        view_bool_1d_from_u8(py::cast(this), v.strands, num_docs()));
  }

private:
  std::unique_ptr<mumemto_mum_result, MumHandleDeleter> h_;
};

class MemResult {
public:
  explicit MemResult(mumemto_mem_result* h) : h_(h) {}

  size_t num_docs() const { return ::num_docs_mem(h_.get()); }
  size_t num_matches() const { return ::num_mems(h_.get()); }

  py::array doc_record_offsets() const {
    const size_t nd = num_docs();
    return view_1d<size_t>(py::cast(this), ::doc_record_offsets_mem(h_.get()), nd + 1);
  }

  py::array record_lengths() const {
    const size_t nd = num_docs();
    const size_t* off = ::doc_record_offsets_mem(h_.get());
    const size_t total = (off && nd > 0) ? off[nd] : 0;
    return view_1d<size_t>(py::cast(this), ::record_lengths_mem(h_.get()), total);
  }

  py::tuple match_at(size_t idx) const {
    const auto v = ::mem_at(h_.get(), idx);
    // (length, offsets[occ], seq_ids[occ], strands[occ])
    return py::make_tuple(
        v.length,
        view_1d<int64_t>(py::cast(this), v.offsets, v.occurrences),
        view_1d<size_t>(py::cast(this), v.seq_ids, v.occurrences),
        view_bool_1d_from_u8(py::cast(this), v.strands, v.occurrences));
  }

private:
  std::unique_ptr<mumemto_mem_result, MemHandleDeleter> h_;
};

MumResult mum_run(const std::vector<std::vector<std::string>>& sequences,
                  uint32_t min_match_len,
                  bool use_revcomp,
                  size_t num_distinct,
                  bool use_gsacak) {
  std::vector<std::vector<const char*>> record_ptrs;
  record_ptrs.reserve(sequences.size());
  for (const auto& doc : sequences) {
    std::vector<const char*> recs;
    recs.reserve(doc.size());
    for (const auto& s : doc) recs.push_back(s.c_str());
    record_ptrs.push_back(std::move(recs));
  }

  std::vector<mumemto_doc_view> docs;
  docs.reserve(sequences.size());
  for (size_t d = 0; d < sequences.size(); ++d) {
    mumemto_doc_view dv;
    dv.records = record_ptrs[d].data();
    dv.num_records = record_ptrs[d].size();
    docs.push_back(dv);
  }

  mumemto_mum_result* out = nullptr;
  const int rc = ::mumemto_mum(docs.data(),
                               docs.size(),
                               min_match_len,
                               as_u8(use_revcomp),
                               num_distinct,
                               as_u8(use_gsacak),
                               &out);
  if (rc != 0) throw std::runtime_error(::mumemto_last_error());
  return MumResult(out);
}

MemResult mem_run(const std::vector<std::vector<std::string>>& sequences,
                  uint32_t min_match_len,
                  bool use_revcomp,
                  size_t num_distinct,
                  size_t max_total_freq,
                  size_t max_doc_freq,
                  bool use_gsacak) {
  std::vector<std::vector<const char*>> record_ptrs;
  record_ptrs.reserve(sequences.size());
  for (const auto& doc : sequences) {
    std::vector<const char*> recs;
    recs.reserve(doc.size());
    for (const auto& s : doc) recs.push_back(s.c_str());
    record_ptrs.push_back(std::move(recs));
  }

  std::vector<mumemto_doc_view> docs;
  docs.reserve(sequences.size());
  for (size_t d = 0; d < sequences.size(); ++d) {
    mumemto_doc_view dv;
    dv.records = record_ptrs[d].data();
    dv.num_records = record_ptrs[d].size();
    docs.push_back(dv);
  }

  mumemto_mem_result* out = nullptr;
  const int rc = ::mumemto_mem(docs.data(),
                               docs.size(),
                               min_match_len,
                               as_u8(use_revcomp),
                               num_distinct,
                               max_total_freq,
                               max_doc_freq,
                               as_u8(use_gsacak),
                               &out);
  if (rc != 0) throw std::runtime_error(::mumemto_last_error());
  return MemResult(out);
}

} // namespace

PYBIND11_MODULE(_mumemto_core, m) {
  m.doc() = "mumemto Python bindings (simple API, backed by libmumemto.so)";

  py::class_<MumResult>(m, "MumResult")
      .def("num_docs", &MumResult::num_docs)
      .def("num_matches", &MumResult::num_matches)
      .def("__len__", &MumResult::num_matches)
      .def("__getitem__",
           [](const MumResult& self, py::ssize_t i) {
             const py::ssize_t n = static_cast<py::ssize_t>(self.num_matches());
             if (i < 0) i += n;
             if (i < 0 || i >= n) throw py::index_error("MumResult index out of range");
             return self.match_at(static_cast<size_t>(i));
           },
           py::arg("i"))
      .def("doc_record_offsets", &MumResult::doc_record_offsets)
      .def("record_lengths", &MumResult::record_lengths)
      .def("match_at", &MumResult::match_at, py::arg("idx"),
           "Return (length, offsets_ndarray, strands_ndarray) as zero-copy views.");

  py::class_<MemResult>(m, "MemResult")
      .def("num_docs", &MemResult::num_docs)
      .def("num_matches", &MemResult::num_matches)
      .def("__len__", &MemResult::num_matches)
      .def("__getitem__",
           [](const MemResult& self, py::ssize_t i) {
             const py::ssize_t n = static_cast<py::ssize_t>(self.num_matches());
             if (i < 0) i += n;
             if (i < 0 || i >= n) throw py::index_error("MemResult index out of range");
             return self.match_at(static_cast<size_t>(i));
           },
           py::arg("i"))
      .def("doc_record_offsets", &MemResult::doc_record_offsets)
      .def("record_lengths", &MemResult::record_lengths)
      .def("match_at", &MemResult::match_at, py::arg("idx"),
           "Return (length, offsets_ndarray, seq_ids_ndarray, strands_ndarray) as zero-copy views.");

  m.def("mum",
        &mum_run,
        py::arg("sequences"),
        py::arg("min_match_len") = 20,
        py::arg("use_revcomp") = true,
        py::arg("num_distinct") = static_cast<size_t>(0),
        py::arg("use_gsacak") = false);

  m.def("mem",
        &mem_run,
        py::arg("sequences"),
        py::arg("min_match_len") = 20,
        py::arg("use_revcomp") = true,
        py::arg("num_distinct") = static_cast<size_t>(0),
        py::arg("max_total_freq") = static_cast<size_t>(0),
        py::arg("max_doc_freq") = static_cast<size_t>(2),
        py::arg("use_gsacak") = false);
}

