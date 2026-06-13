/*
 * Bindings for the file-type detector: the Filetype / CompressionType enums and
 * FiletypeDetector. Mirrors genogrove/io/filetype_detector.hpp.
 *
 * FiletypeDetector.detect_filetype(path) returns a (Filetype, CompressionType)
 * tuple, inferred from the path's extension and magic bytes.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // std::tuple -> Python tuple

#include <string>

#include <genogrove/io/filetype_detector.hpp>

namespace py = pybind11;
namespace gio = genogrove::io;

inline void bind_filetype_detector(py::module_& m) {
    py::enum_<gio::filetype>(m, "Filetype", "Detected file format.")
        .value("BED", gio::filetype::BED)
        .value("BEDGRAPH", gio::filetype::BEDGRAPH)
        .value("GFF", gio::filetype::GFF)
        .value("GTF", gio::filetype::GTF)
        .value("VCF", gio::filetype::VCF)
        .value("SAM", gio::filetype::SAM)
        .value("BAM", gio::filetype::BAM)
        .value("FASTA", gio::filetype::FASTA)
        .value("FASTQ", gio::filetype::FASTQ)
        .value("GG", gio::filetype::GG)
        .value("UNKNOWN", gio::filetype::UNKNOWN);

    py::enum_<gio::compression_type>(m, "CompressionType", "Detected compression.")
        .value("NONE", gio::compression_type::NONE)
        .value("GZIP", gio::compression_type::GZIP)
        .value("BZIP2", gio::compression_type::BZIP2)
        .value("XZ", gio::compression_type::XZ)
        .value("ZSTD", gio::compression_type::ZSTD)
        .value("LZ4", gio::compression_type::LZ4)
        .value("UNKNOWN", gio::compression_type::UNKNOWN);

    py::class_<gio::filetype_detector>(m, "FiletypeDetector", R"pbdoc(
        Detects a file's format and compression from its extension + magic bytes.
    )pbdoc")
        .def(py::init<>())
        .def("detect_filetype",
             [](gio::filetype_detector& d, const std::string& path) {
                 return d.detect_filetype(path);  // std::string -> fs::path
             },
             py::arg("path"),
             R"pbdoc(
                 detect_filetype(path) -> (Filetype, CompressionType)

                 Infer the file format and compression of `path`. Returns a tuple,
                 e.g. (Filetype.BED, CompressionType.GZIP) for "peaks.bed.gz".
                 Unrecognized inputs yield Filetype.UNKNOWN / CompressionType.UNKNOWN.
             )pbdoc");
}