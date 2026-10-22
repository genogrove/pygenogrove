/*
 * Bindings for the random-access FASTA index (faidx). Mirrors
 * genogrove/io/fasta_index.hpp — a thin wrapper over htslib's faidx API.
 *
 * FastaIndex(path) opens a FASTA file and loads (or creates) its .fai index,
 * then serves region/whole-sequence fetches and per-sequence metadata. It
 * pairs with GenomicCoordinate: fetch a feature's bases with
 * idx.fetch(coord_index, gc.start, gc.end + 1) — fetch() is 0-based half-open
 * [start, end), GenomicCoordinate is 0-based closed [start, end].
 *
 * Every method releases the GIL (faidx is pure htslib I/O), so a handle
 * shared between Python threads would otherwise race on the one BGZF cursor
 * inside the faidx_t. The bound type is therefore a small wrapper pairing the
 * index with a mutex; each method takes the lock only after dropping the GIL
 * and drops it again before the GIL is reacquired, so lock order is always
 * GIL -> (release) -> mutex and there is no inversion to deadlock on.
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

struct fasta_index_handle {
    gio::fasta_index idx;
    mutable std::mutex mu;

    explicit fasta_index_handle(const std::string& path) : idx(path) {}

    template <typename F>
    auto locked(F&& f) const {
        std::lock_guard<std::mutex> lock(mu);
        return f(idx);
    }
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

        Thread-safe: one handle may be shared between threads; calls on it are
        serialized by an internal lock (the GIL is released for the htslib
        work). For parallel throughput open one handle per thread.

        Non-copyable; the underlying htslib handle is closed when the object is
        garbage-collected.
    )pbdoc")
        .def(py::init([](const std::string& path) {
                 // std::string -> std::filesystem::path; throws std::runtime_error
                 // (-> RuntimeError) if the file can't be opened or indexed.
                 //
                 // htslib's fai_load() is not thread-safe, so concurrent opens
                 // race and abort. Serialize construction with a
                 // process-wide lock; the GIL stays released (below) so a long
                 // multi-GB index build still won't freeze unrelated threads.
                 static std::mutex fai_load_mutex;
                 std::lock_guard<std::mutex> lock(fai_load_mutex);
                 return std::make_unique<handle_t>(path);
             }),
             py::arg("path"),
             // Building a missing .fai for a multi-GB genome is pure htslib I/O.
             py::call_guard<py::gil_scoped_release>(),
             "Open a FASTA file and load (or create) its .fai index.")
        .def("fetch",
             [](const handle_t& h, const std::string& name, std::size_t start,
                std::size_t end) {
                 return h.locked([&](const gio::fasta_index& idx) {
                     return idx.fetch(name, start, end);
                 });
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
                 return h.locked(
                     [&](const gio::fasta_index& idx) { return idx.fetch(name); });
             },
             py::arg("name"),
             py::call_guard<py::gil_scoped_release>(),
             R"pbdoc(
                 fetch(name) -> str

                 The entire sequence named `name`. Raises IndexError if unknown.
             )pbdoc")
        .def("sequence_count",
             [](const handle_t& h) {
                 return h.locked(
                     [](const gio::fasta_index& idx) { return idx.sequence_count(); });
             },
             py::call_guard<py::gil_scoped_release>(),
             "Number of sequences in the index.")
        .def("sequence_name",
             [](const handle_t& h, std::size_t index) {
                 return h.locked([&](const gio::fasta_index& idx) {
                     return idx.sequence_name(index);
                 });
             },
             py::arg("index"),
             py::call_guard<py::gil_scoped_release>(),
             "Name of the i-th sequence (0-based). Raises IndexError if out of range.")
        .def("sequence_length",
             [](const handle_t& h, const std::string& name) {
                 return h.locked([&](const gio::fasta_index& idx) {
                     return idx.sequence_length(name);
                 });
             },
             py::arg("name"),
             py::call_guard<py::gil_scoped_release>(),
             "Length in bases of sequence `name`. Raises IndexError if unknown.")
        .def("has_sequence",
             [](const handle_t& h, const std::string& name) {
                 return h.locked([&](const gio::fasta_index& idx) {
                     return idx.has_sequence(name);
                 });
             },
             py::arg("name"),
             py::call_guard<py::gil_scoped_release>(),
             "Whether `name` is present in the index.")
        .def("names",
             [](const handle_t& h) {
                 return h.locked([](const gio::fasta_index& idx) {
                     std::vector<std::string> out;
                     out.reserve(idx.sequence_count());
                     for (std::size_t i = 0; i < idx.sequence_count(); ++i) {
                         out.push_back(idx.sequence_name(i));
                     }
                     return out;
                 });
             },
             py::call_guard<py::gil_scoped_release>(),
             "List of all sequence names, in index order.")
        .def("__len__",
             [](const handle_t& h) {
                 return h.locked(
                     [](const gio::fasta_index& idx) { return idx.sequence_count(); });
             },
             py::call_guard<py::gil_scoped_release>())
        .def("__contains__",
             [](const handle_t& h, const std::string& name) {
                 return h.locked([&](const gio::fasta_index& idx) {
                     return idx.has_sequence(name);
                 });
             },
             py::arg("name"),
             py::call_guard<py::gil_scoped_release>())
        .def("__repr__",
             [](const handle_t& h) {
                 return "FastaIndex(n_sequences=" +
                        std::to_string(h.locked([](const gio::fasta_index& idx) {
                            return idx.sequence_count();
                        })) +
                        ")";
             },
             py::call_guard<py::gil_scoped_release>());
}