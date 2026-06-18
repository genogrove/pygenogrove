/*
 * Bindings for the SAM/BAM alignment reader: SamFlags / AlignmentFlags /
 * SamEntry (value types) and BamReader (single-pass iterator). Mirrors
 * genogrove/io/bam_reader.hpp.
 *
 * sam_entry is not serializable (it holds optionals / a tag variant / a CIGAR
 * vector), so there is no typed `BamGrove`. The intended flow is to load
 * alignments into the universal Grove as JSON: build the strand-aware key with
 * SamEntry.to_coordinate() and a payload with SamEntry.to_dict() (or your own
 * dict).
 *
 * v1 scope: the core record + flags + CIGAR string form + reader. The CIGAR
 * element list, paired-end mate info, and auxiliary tags are not yet exposed.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include <genogrove/data_type/genomic_coordinate.hpp>
#include <genogrove/io/bam_reader.hpp>

namespace py = pybind11;
namespace gio = genogrove::io;
namespace gdt = genogrove::data_type;

inline void bind_sam_entry(py::module_& m) {
    // SAM FLAG bit constants (e.g. SamFlags.REVERSE) for use with has_flag().
    py::class_<gio::sam_flags> sf(m, "SamFlags",
                                  "SAM/BAM FLAG bit constants (for has_flag()).");
    sf.attr("PAIRED") = gio::sam_flags::PAIRED;
    sf.attr("PROPER_PAIR") = gio::sam_flags::PROPER_PAIR;
    sf.attr("UNMAPPED") = gio::sam_flags::UNMAPPED;
    sf.attr("MATE_UNMAPPED") = gio::sam_flags::MATE_UNMAPPED;
    sf.attr("REVERSE") = gio::sam_flags::REVERSE;
    sf.attr("MATE_REVERSE") = gio::sam_flags::MATE_REVERSE;
    sf.attr("READ1") = gio::sam_flags::READ1;
    sf.attr("READ2") = gio::sam_flags::READ2;
    sf.attr("SECONDARY") = gio::sam_flags::SECONDARY;
    sf.attr("QCFAIL") = gio::sam_flags::QCFAIL;
    sf.attr("DUPLICATE") = gio::sam_flags::DUPLICATE;
    sf.attr("SUPPLEMENTARY") = gio::sam_flags::SUPPLEMENTARY;

    py::class_<gio::alignment_flags>(m, "AlignmentFlags", R"pbdoc(
        The SAM/BAM bitwise FLAG, with convenience accessors.
    )pbdoc")
        .def(py::init<>())
        .def(py::init<uint16_t>(), py::arg("flags"))
        .def("value", &gio::alignment_flags::value, "Raw 16-bit FLAG value.")
        .def("has_flag", &gio::alignment_flags::has_flag, py::arg("flag"),
             "Whether a specific FLAG bit (see SamFlags) is set.")
        .def("is_paired", &gio::alignment_flags::is_paired)
        .def("is_proper_pair", &gio::alignment_flags::is_proper_pair)
        .def("is_unmapped", &gio::alignment_flags::is_unmapped)
        .def("is_mate_unmapped", &gio::alignment_flags::is_mate_unmapped)
        .def("is_reverse", &gio::alignment_flags::is_reverse)
        .def("is_mate_reverse", &gio::alignment_flags::is_mate_reverse)
        .def("is_read1", &gio::alignment_flags::is_read1)
        .def("is_read2", &gio::alignment_flags::is_read2)
        .def("is_secondary", &gio::alignment_flags::is_secondary)
        .def("is_qc_fail", &gio::alignment_flags::is_qc_fail)
        .def("is_duplicate", &gio::alignment_flags::is_duplicate)
        .def("is_supplementary", &gio::alignment_flags::is_supplementary)
        .def("__repr__", [](const gio::alignment_flags& f) {
            return "AlignmentFlags(" + std::to_string(f.value()) + ")";
        });

    py::class_<gio::sam_entry>(m, "SamEntry", R"pbdoc(
        A single SAM/BAM alignment record (htslib-native coordinates: `start` is
        0-based, `end` is 0-based exclusive).

        Load alignments into the universal Grove as JSON, deriving the strand-aware
        key from the record::

            g = pygenogrove.Grove()
            for aln in pygenogrove.BamReader("reads.bam"):
                g.insert(aln.chrom, aln.to_coordinate(), aln.to_dict())
    )pbdoc")
        .def(py::init<>())
        .def_readwrite("qname", &gio::sam_entry::qname, "Query/read name (QNAME)")
        .def_readwrite("chrom", &gio::sam_entry::chrom, "Reference name (RNAME)")
        .def_readwrite("start", &gio::sam_entry::start,
                       "0-based start (htslib-native, from POS)")
        .def_readwrite("end", &gio::sam_entry::end,
                       "0-based exclusive end (POS + CIGAR ref length)")
        .def_readwrite("mapq", &gio::sam_entry::mapq, "Mapping quality (0-255)")
        .def_readwrite("sequence", &gio::sam_entry::sequence, "Read sequence (SEQ)")
        .def_readwrite("quality", &gio::sam_entry::quality,
                       "ASCII quality scores (QUAL)")
        .def_readwrite("flags", &gio::sam_entry::flags, "The AlignmentFlags (FLAG)")
        .def_property_readonly(
            "cigar",
            [](const gio::sam_entry& e) { return e.cigar_string_repr(); },
            "CIGAR string form, e.g. '100M' ('*' if none).")
        .def("get_strand", &gio::sam_entry::get_strand,
             "Strand char from the FLAG: '+' forward, '-' reverse, '.' unmapped.")
        .def("is_primary", &gio::sam_entry::is_primary,
             "Not secondary and not supplementary.")
        .def("is_mapped", &gio::sam_entry::is_mapped, "The read is mapped.")
        .def("consumes_reference", &gio::sam_entry::consumes_reference,
             "Covers >= 1 reference base (False for unmapped / pure soft-clip).")
        .def("to_coordinate",
             [](const gio::sam_entry& e) {
                 if (!e.consumes_reference()) {
                     throw std::invalid_argument(
                         "SamEntry covers no reference bases (unmapped or "
                         "zero-ref-consuming CIGAR); cannot derive a coordinate");
                 }
                 // htslib half-open [start, end) -> closed [start, end - 1]
                 return gdt::genomic_coordinate(e.get_strand(), e.start, e.end - 1);
             },
             R"pbdoc(
                 Derive the strand-aware GenomicCoordinate key for this alignment
                 (strand from the FLAG; 0-based half-open [start, end) -> closed
                 [start, end-1]). Raises ValueError if the read covers no reference
                 bases (unmapped or a zero-ref-consuming CIGAR).
             )pbdoc")
        .def("to_dict",
             [](const gio::sam_entry& e) {
                 py::dict d;
                 d["qname"] = e.qname;
                 d["mapq"] = e.mapq;
                 d["strand"] = std::string(1, e.get_strand());
                 d["cigar"] = e.cigar_string_repr();
                 d["flags"] = e.flags.value();
                 d["is_primary"] = e.is_primary();
                 d["is_mapped"] = e.is_mapped();
                 return d;
             },
             "A dict of the core alignment fields (qname, mapq, strand, cigar, "
             "flags, is_primary, is_mapped) — convenient as a Grove JSON payload.")
        .def("__repr__", [](const gio::sam_entry& e) {
            return "SamEntry(qname='" + e.qname + "', chrom='" + e.chrom +
                   "', start=" + std::to_string(e.start) +
                   ", end=" + std::to_string(e.end) +
                   ", strand='" + std::string(1, e.get_strand()) + "')";
        });
}

inline void bind_bam_reader(py::module_& m) {
    // SamEntry must already be registered (bind_sam_entry) — BamReader yields it.
    py::class_<gio::bam_reader>(m, "BamReader", R"pbdoc(
        A single-pass iterator over the alignments of a SAM/BAM file.

        Iterate it directly to get SamEntry objects::

            for aln in pygenogrove.BamReader("reads.bam"):
                ...

        SAM/BAM are auto-detected by htslib. The reader owns an htslib handle and
        is single-pass — it cannot be restarted or iterated twice.

        Parameters
        ----------
        path : str
            Path to the SAM/BAM file. A missing/unreadable file raises.
        skip_unmapped : bool, optional
            Skip unmapped reads (default True).
        skip_secondary, skip_supplementary, skip_qc_fail, skip_duplicates : bool, optional
            Skip the corresponding alignment categories (default False).
        min_mapq : int, optional
            Minimum mapping quality; reads below it are skipped (0 = no filter).
    )pbdoc")
        .def(py::init([](const std::string& path, bool skip_unmapped,
                         bool skip_secondary, bool skip_supplementary,
                         bool skip_qc_fail, bool skip_duplicates,
                         uint8_t min_mapq) {
                 gio::bam_reader_options opts;
                 opts.skip_unmapped = skip_unmapped;
                 opts.skip_secondary = skip_secondary;
                 opts.skip_supplementary = skip_supplementary;
                 opts.skip_qc_fail = skip_qc_fail;
                 opts.skip_duplicates = skip_duplicates;
                 opts.min_mapq = min_mapq;
                 return std::make_unique<gio::bam_reader>(path, opts);
             }),
             py::arg("path"), py::arg("skip_unmapped") = true,
             py::arg("skip_secondary") = false,
             py::arg("skip_supplementary") = false,
             py::arg("skip_qc_fail") = false, py::arg("skip_duplicates") = false,
             py::arg("min_mapq") = 0)
        .def("__iter__", [](gio::bam_reader& r) -> gio::bam_reader& { return r; })
        .def("__next__",
             [](gio::bam_reader& r) {
                 gio::sam_entry entry;
                 if (!r.read_next(entry)) {
                     throw py::stop_iteration();
                 }
                 return entry;
             },
             // The read (htslib decode / disk I/O) touches no Python objects;
             // pybind reacquires the GIL before converting the returned entry.
             py::call_guard<py::gil_scoped_release>())
        .def("get_error_message", &gio::bam_reader::get_error_message,
             "Error message from the most recent read; empty on clean EOF.")
        .def("get_current_line", &gio::bam_reader::get_current_line,
             "Records consumed so far (advances on skipped/filtered records too).");
}