#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "mumemto_api.hpp"

namespace py = pybind11;

namespace {

template <typename T>
py::array view_1d(py::handle base, const T* ptr, size_t n) {
  if (!ptr) {
    return py::array(py::dtype::of<T>(), py::array::ShapeContainer{static_cast<py::ssize_t>(0)});
  }
  py::array a = py::array(
      py::dtype::of<T>(),
      {static_cast<py::ssize_t>(n)},
      {static_cast<py::ssize_t>(sizeof(T))},
      ptr,
      base);
  a.attr("setflags")(false);
  return a;
}

inline py::array view_bool_1d_from_u8(py::handle base, const uint8_t* ptr, size_t n) {
  if (!ptr) {
    return py::array(py::dtype("bool"), py::array::ShapeContainer{static_cast<py::ssize_t>(0)});
  }
  // NumPy bool_ is stored as 1 byte; interpret the existing 0/1 bytes as bool.
  py::array a = py::array(
      py::dtype("bool"),
      {static_cast<py::ssize_t>(n)},
      {static_cast<py::ssize_t>(sizeof(uint8_t))},
      ptr,
      base);
  a.attr("setflags")(false);
  return a;
}

class MumResult {
public:
  explicit MumResult(mumemto::MumResult&& r) : r_(std::move(r)) {}

  size_t num_docs() const { return r_.lengths.size(); }
  size_t num_matches() const { return r_.matches.size(); }

  py::tuple match_at(size_t idx) const {
    if (idx >= r_.matches.size()) {
      throw py::index_error("MumResult index out of range");
    }
    const auto& m = r_.matches[idx];
    return py::make_tuple(
        m.length,
        view_1d<int64_t>(py::cast(this), m.offsets.data(), m.offsets.size()),
        view_bool_1d_from_u8(py::cast(this), m.strands.data(), m.strands.size()));
  }

private:
  mumemto::MumResult r_;
};

class MemResult {
public:
  explicit MemResult(mumemto::MemResult&& r) : r_(std::move(r)) {}

  size_t num_docs() const { return r_.lengths.size(); }
  size_t num_matches() const { return r_.matches.size(); }

  py::tuple match_at(size_t idx) const {
    if (idx >= r_.matches.size()) {
      throw py::index_error("MemResult index out of range");
    }
    const auto& m = r_.matches[idx];
    return py::make_tuple(
        m.length,
        view_1d<int64_t>(py::cast(this), m.offsets.data(), m.offsets.size()),
        view_1d<size_t>(py::cast(this), m.seq_ids.data(), m.seq_ids.size()),
        view_bool_1d_from_u8(py::cast(this), m.strands.data(), m.strands.size()));
  }

private:
  mumemto::MemResult r_;
};

MumResult mum_run(std::vector<std::vector<std::string>> sequences,
                  uint32_t min_match_len,
                  bool use_revcomp,
                  size_t num_distinct,
                  bool use_gsacak) {
  // `sequences` is already the pybind11 materialized buffer; forward to native C++ API without another copy.
  auto r = mumemto::mumemto_mum(sequences, min_match_len, use_revcomp, num_distinct, use_gsacak);
  return MumResult(std::move(r));
}

MemResult mem_run(std::vector<std::vector<std::string>> sequences,
                  uint32_t min_match_len,
                  bool use_revcomp,
                  size_t num_distinct,
                  size_t max_total_freq,
                  size_t max_doc_freq,
                  bool use_gsacak) {
  auto r = mumemto::mumemto_mem(sequences,
                                min_match_len,
                                use_revcomp,
                                num_distinct,
                                max_total_freq,
                                max_doc_freq,
                                use_gsacak);
  return MemResult(std::move(r));
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

