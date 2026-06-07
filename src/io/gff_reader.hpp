/*
 * Bindings for the GFF/GTF value types defined in genogrove io/gff_reader.hpp:
 * the gff_entry record, its detected gff_format, and the GTF helper accessors.
 *
 * Exposed so a data-carrying grove (grove<interval, gff_entry>) can store,
 * query, and serialize GFF/GTF records. The gff_reader file iterator itself is
 * bound separately (future work).
 *
 * NOTE: <genogrove/io/gff_reader.hpp> transitively pulls in <htslib/bgzf.h>,
 * but htslib is already a transitive build dependency. We only touch the POD
 * struct/enum here, never the reader/BGZF handle.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <string_view>

#include <genogrove/io/gff_reader.hpp>

namespace py = pybind11;
namespace gio = genogrove::io;

inline void bind_gff_entry(py::module_& m) {
    py::enum_<gio::gff_format>(m, "GffFormat",
        "The format of a GFF record, auto-detected during parsing.")
        .value("GFF3", gio::gff_format::GFF3, "GFF3 (key=value attributes)")
        .value("GTF", gio::gff_format::GTF, "GTF/GTF2 (key \"value\" attributes)")
        .value("UNKNOWN", gio::gff_format::UNKNOWN, "Format not yet determined");

    py::class_<gio::gff_entry>(m, "GffEntry", R"pbdoc(
        A single GFF3/GTF record (one line of a GFF file).

        Coordinates are GFF-native: 1-based, both endpoints inclusive. This is
        distinct from both `Interval` (0-based closed) and `BedEntry` (0-based
        half-open), so convert deliberately when building the grove key.
        Optional columns are None when absent ('.').

        Parameters
        ----------
        seqid : str
            Sequence/chromosome name (column 1).
        start : int
            1-based inclusive start (column 4).
        end : int
            1-based inclusive end (column 5).
        type : str
            Feature type, e.g. "gene", "exon", "CDS" (column 3).
    )pbdoc")
        .def(py::init<>())
        .def(py::init<std::string, size_t, size_t, std::string>(),
             py::arg("seqid"), py::arg("start"), py::arg("end"), py::arg("type"))
        .def_readwrite("seqid", &gio::gff_entry::seqid,
                       "Sequence/chromosome name (column 1)")
        .def_readwrite("source", &gio::gff_entry::source,
                       "Source of the feature (column 2)")
        .def_readwrite("type", &gio::gff_entry::type,
                       "Feature type, e.g. gene/exon/CDS (column 3)")
        .def_readwrite("start", &gio::gff_entry::start,
                       "1-based inclusive start (column 4)")
        .def_readwrite("end", &gio::gff_entry::end,
                       "1-based inclusive end (column 5)")
        .def_readwrite("score", &gio::gff_entry::score,
                       "Optional[float] score (column 6)")
        .def_readwrite("strand", &gio::gff_entry::strand,
                       "Optional[str] strand — a single character ('+', '-', "
                       "'.', '?'); assigning an empty or multi-character string "
                       "raises ValueError, None clears it (column 7)")
        .def_readwrite("phase", &gio::gff_entry::phase,
                       "Optional[int] CDS phase (0, 1, or 2; column 8)")
        .def_readwrite("attributes", &gio::gff_entry::attributes,
                       "dict[str, str] of column-9 key/value attributes "
                       "(assigned/returned by copy)")
        .def_readwrite("format", &gio::gff_entry::format,
                       "GffFormat detected for this record")
        .def("is_gtf", &gio::gff_entry::is_gtf,
             "True if this record was parsed as GTF.")
        .def("is_gff3", &gio::gff_entry::is_gff3,
             "True if this record was parsed as GFF3.")
        .def("get_attribute", &gio::gff_entry::get_attribute, py::arg("key"),
             "Return the value of a column-9 attribute, or None if absent.")
        .def("get_gene_id", &gio::gff_entry::get_gene_id,
             "GTF `gene_id` attribute, or None.")
        .def("get_transcript_id", &gio::gff_entry::get_transcript_id,
             "GTF `transcript_id` attribute, or None.")
        .def("get_exon_number", &gio::gff_entry::get_exon_number,
             "GTF `exon_number` attribute as int, or None.")
        .def("get_gene_name", &gio::gff_entry::get_gene_name,
             "GTF `gene_name` attribute, or None.")
        .def("get_gene_biotype", &gio::gff_entry::get_gene_biotype,
             "GTF `gene_biotype`/`gene_type` attribute, or None.")
        .def("__repr__", [](const gio::gff_entry& e) {
            return "GffEntry(seqid='" + e.seqid + "', start=" +
                   std::to_string(e.start) + ", end=" + std::to_string(e.end) +
                   ", type='" + e.type + "')";
        });
}