/*
 * Bindings for the FASTA/FASTQ sequence reader: FastaEntry (value type) and
 * FastaReader (single-pass iterator). Mirrors genogrove/io/fasta_reader.hpp.
 *
 * FASTA records are named sequences, not genomic intervals, so this is a
 * standalone reader (no grove integration / coordinate derivation). The reader
 * handles both FASTA (`>` headers) and FASTQ (`@` headers + per-base quality),
 * auto-detected per file. Random-access (fasta_index) is deferred.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // std::optional<std::string> quality -> str | None

#include <memory>
#include <string>

#include <genogrove/io/fasta_reader.hpp>

namespace py = pybind11;
namespace gio = genogrove::io;

inline void bind_fasta_entry(py::module_& m) {
    py::class_<gio::fasta_entry>(m, "FastaEntry", R"pbdoc(
        A single FASTA/FASTQ record: a named nucleotide sequence.

        For FASTQ records, `quality` holds the per-base quality string; for FASTA
        it is None. `name` is the header text up to the first whitespace, and
        `comment` is the rest of the header line.
    )pbdoc")
        .def(py::init<>())
        .def(py::init<std::string, std::string>(),
             py::arg("name"), py::arg("sequence"))
        .def_readwrite("name", &gio::fasta_entry::name,
                       "Sequence name (header text up to the first whitespace).")
        .def_readwrite("comment", &gio::fasta_entry::comment,
                       "Optional description (rest of the header line).")
        .def_readwrite("sequence", &gio::fasta_entry::sequence,
                       "Nucleotide sequence.")
        .def_readwrite("quality", &gio::fasta_entry::quality,
                       "Per-base quality string (FASTQ only; None for FASTA).")
        .def("is_fastq",
             [](const gio::fasta_entry& e) { return e.quality.has_value(); },
             "Whether this record carries quality scores (i.e. came from FASTQ).")
        .def("__len__",
             [](const gio::fasta_entry& e) { return e.sequence.size(); },
             "Length of the sequence.")
        .def("__repr__", [](const gio::fasta_entry& e) {
            return "FastaEntry(name='" + e.name + "', len=" +
                   std::to_string(e.sequence.size()) +
                   (e.quality ? ", fastq=True)" : ")");
        });
}

inline void bind_fasta_reader(py::module_& m) {
    // FastaEntry must already be registered (bind_fasta_entry).
    py::class_<gio::fasta_reader>(m, "FastaReader", R"pbdoc(
        A single-pass iterator over the records of a FASTA or FASTQ file.

        Iterate it directly to get FastaEntry objects::

            for rec in pygenogrove.FastaReader("genome.fa"):
                print(rec.name, len(rec))

        FASTA (`>`) and FASTQ (`@`) are auto-detected, and plain or
        gzip/BGZF-compressed (`.gz`) inputs are accepted. The reader owns an
        htslib handle and is single-pass — it cannot be restarted or iterated twice.

        Parameters
        ----------
        path : str
            Path to the FASTA/FASTQ file. A missing/unreadable file raises.
        skip_empty_sequences : bool, optional
            Skip records whose sequence is empty (default False).
    )pbdoc")
        .def(py::init([](const std::string& path, bool skip_empty_sequences) {
                 gio::fasta_reader_options opts;
                 opts.skip_empty_sequences = skip_empty_sequences;
                 return std::make_unique<gio::fasta_reader>(path, opts);
             }),
             py::arg("path"), py::arg("skip_empty_sequences") = false)
        .def("__iter__",
             [](gio::fasta_reader& r) -> gio::fasta_reader& { return r; })
        .def("__next__",
             [](gio::fasta_reader& r) {
                 gio::fasta_entry entry;
                 if (!r.read_next(entry)) {
                     throw py::stop_iteration();
                 }
                 return entry;
             })
        .def("get_error_message", &gio::fasta_reader::get_error_message,
             "Error message from the most recent read; empty on clean EOF.")
        .def("get_current_line", &gio::fasta_reader::get_current_line,
             "1-based physical line number consumed so far; 0 before the first read.");
}