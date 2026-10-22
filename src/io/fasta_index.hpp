/*
 * Bindings for the random-access FASTA index (faidx). Mirrors
 * genogrove/io/fasta_index.hpp.
 *
 * fetch() releases the GIL and shares one BGZF cursor per handle, so it is
 * serialized by a per-handle mutex taken only after the GIL is released. The
 * metadata getters are read-only lookups and stay unlocked.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // std::vector<std::string> -> Python list

#include <cstddef>
#include <mutex>
#include <string>
#include <vector>

#include <genogrove/io/fasta_index.hpp>

namespace py = pybind11;
namespace gio = genogrove::io;

namespace pygg {

struct fasta_index_handle : gio::fasta_index {
    using gio::fasta_index::fasta_index;
    mutable std::mutex fetch_mutex;
};

}  // namespace pygg

inline void bind_fasta_index(py::module_& m) {
    using handle_t = pygg::fasta_index_handle;

    py::class_<handle_t>(m, "FastaIndex", R"pbdoc(
        Random-access reader for a FASTA file, backed by an .fai index.

        Opening builds the index if it is missing (this writes a sibling
        ``.fai`` file, so the FASTA's directory must be writable on first open).
        Region coordinates are 0-based half-open ``[start, end)`` — to fetch the
        bases of a GenomicCoordinate (0-based closed ``[start, end]``), call
        ``idx.fetch(index, gc.start, gc.end + 1)``.

        One handle may be shared between threads: fetch() calls on it are
        serialized by an internal lock (the GIL is released for the htslib
        work); the metadata getters are lock-free reads. For parallel fetch
        throughput open one handle per thread.

        Non-copyable; the underlying htslib handle is closed when the object is
        garbage-collected.
    )pbdoc")
        .def(py::init([](const std::string& path) {
                 // fai_load() is not thread-safe: serialize construction.
                 static std::mutex fai_load_mutex;
                 std::lock_guard<std::mutex> lock(fai_load_mutex);
                 return std::make_unique<handle_t>(path);
             }),
             py::arg("path"),
             py::call_guard<py::gil_scoped_release>(),
             "Open a FASTA file and load (or create) its .fai index.")
        .def("fetch",
             [](const handle_t& h, const std::string& name, std::size_t start,
                std::size_t end) {
                 std::lock_guard<std::mutex> lock(h.fetch_mutex);
                 return h.fetch(name, start, end);
             },
             py::arg("name"), py::arg("start"), py::arg("end"),
             py::call_guard<py::gil_scoped_release>(),
             R"pbdoc(
                 fetch(name, start, end) -> str

                 Bases of sequence `name` over the 0-based half-open region
                 [start, end). Raises IndexError if `name` is unknown or the
                 region is invalid (start >= end, or beyond htslib's limit).
             )pbdoc")
        .def("fetch",
             [](const handle_t& h, const std::string& name) {
                 std::lock_guard<std::mutex> lock(h.fetch_mutex);
                 return h.fetch(name);
             },
             py::arg("name"),
             py::call_guard<py::gil_scoped_release>(),
             R"pbdoc(
                 fetch(name) -> str

                 The entire sequence named `name`. Raises IndexError if unknown.
             )pbdoc")
        .def("sequence_count", &gio::fasta_index::sequence_count,
             "Number of sequences in the index.")
        .def("sequence_name", &gio::fasta_index::sequence_name, py::arg("index"),
             "Name of the i-th sequence (0-based). Raises IndexError if out of range.")
        .def("sequence_length", &gio::fasta_index::sequence_length,
             py::arg("name"),
             "Length in bases of sequence `name`. Raises IndexError if unknown.")
        .def("has_sequence", &gio::fasta_index::has_sequence, py::arg("name"),
             "Whether `name` is present in the index.")
        .def("names",
             [](const handle_t& idx) {
                 std::vector<std::string> out;
                 out.reserve(idx.sequence_count());
                 for (std::size_t i = 0; i < idx.sequence_count(); ++i) {
                     out.push_back(idx.sequence_name(i));
                 }
                 return out;
             },
             "List of all sequence names, in index order.")
        .def("__len__", &gio::fasta_index::sequence_count)
        .def("__contains__", &gio::fasta_index::has_sequence, py::arg("name"))
        .def("__repr__", [](const handle_t& idx) {
            return "FastaIndex(n_sequences=" + std::to_string(idx.sequence_count()) + ")";
        });
}