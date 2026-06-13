/*
 * Bindings for the random-access FASTA index (faidx). Mirrors
 * genogrove/io/fasta_index.hpp — a thin wrapper over htslib's faidx API.
 *
 * FastaIndex(path) opens a FASTA file and loads (or creates) its .fai index,
 * then serves region/whole-sequence fetches and per-sequence metadata. It
 * pairs with GenomicCoordinate: fetch a feature's bases with
 * idx.fetch(coord_index, gc.start, gc.end + 1) — fetch() is 0-based half-open
 * [start, end), GenomicCoordinate is 0-based closed [start, end].
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // std::vector<std::string> -> Python list

#include <cstddef>
#include <string>
#include <vector>

#include <genogrove/io/fasta_index.hpp>

namespace py = pybind11;
namespace gio = genogrove::io;

inline void bind_fasta_index(py::module_& m) {
    py::class_<gio::fasta_index>(m, "FastaIndex", R"pbdoc(
        Random-access reader for a FASTA file, backed by an .fai index.

        Opening builds the index if it is missing (this writes a sibling
        ``.fai`` file, so the FASTA's directory must be writable on first open).
        Region coordinates are 0-based half-open ``[start, end)`` — to fetch the
        bases of a GenomicCoordinate (0-based closed ``[start, end]``), call
        ``idx.fetch(index, gc.start, gc.end + 1)``.

        Non-copyable; the underlying htslib handle is closed when the object is
        garbage-collected.
    )pbdoc")
        .def(py::init([](const std::string& path) {
                 // std::string -> std::filesystem::path; throws std::runtime_error
                 // (-> RuntimeError) if the file can't be opened or indexed.
                 return std::make_unique<gio::fasta_index>(path);
             }),
             py::arg("path"),
             "Open a FASTA file and load (or create) its .fai index.")
        .def("fetch",
             py::overload_cast<const std::string&, std::size_t, std::size_t>(
                 &gio::fasta_index::fetch, py::const_),
             py::arg("name"), py::arg("start"), py::arg("end"),
             R"pbdoc(
                 fetch(name, start, end) -> str

                 Bases of sequence `name` over the 0-based half-open region
                 [start, end). Raises IndexError if `name` is unknown or the
                 region is invalid (start >= end, or beyond htslib's limit).
             )pbdoc")
        .def("fetch",
             py::overload_cast<const std::string&>(&gio::fasta_index::fetch,
                                                   py::const_),
             py::arg("name"),
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
             [](const gio::fasta_index& idx) {
                 std::vector<std::string> out;
                 out.reserve(idx.sequence_count());
                 for (std::size_t i = 0; i < idx.sequence_count(); ++i) {
                     out.push_back(idx.sequence_name(i));
                 }
                 return out;
             },
             "List of all sequence names, in index order.")
        .def("__len__", &gio::fasta_index::sequence_count)
        .def("__contains__", &gio::fasta_index::has_sequence, py::arg("name"));
}