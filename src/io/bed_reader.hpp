/*
 * Bindings for the BED value types defined in genogrove io/bed_reader.hpp:
 * bed_entry and its sub-structures (block_info, thick_info, rgb_color).
 *
 * These are plain data carriers parsed from BED files, exposed so a
 * data-carrying grove (grove<interval, bed_entry>) can store, query, and
 * serialize them. The bed_reader file iterator itself is bound separately
 * (future work).
 *
 * NOTE: <genogrove/io/bed_reader.hpp> transitively pulls in <htslib/bgzf.h>,
 * but htslib is already a transitive build dependency (genogrove links it
 * PUBLIC). We only touch the POD structs here, never the reader/BGZF handle.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include <genogrove/io/bed_reader.hpp>

namespace py = pybind11;
namespace gio = genogrove::io;

inline void bind_bed_entry(py::module_& m) {
    // Register the inner structs BEFORE BedEntry, whose optional fields
    // reference them as bound types.

    py::class_<gio::block_info>(m, "BlockInfo", R"pbdoc(
        BED12 block (exon) structure: parallel arrays of block sizes and
        block starts (starts are relative to the entry's chromStart).
    )pbdoc")
        .def(py::init<>())
        .def(py::init<int, std::vector<size_t>, std::vector<size_t>>(),
             py::arg("count"), py::arg("sizes"), py::arg("starts"))
        .def_readwrite("count", &gio::block_info::count,
                       "Number of blocks (BED blockCount)")
        .def_readwrite("sizes", &gio::block_info::sizes,
                       "Block sizes (list[int]; assigned/returned by copy)")
        .def_readwrite("starts", &gio::block_info::starts,
                       "Block starts relative to chromStart (list[int]; by copy)")
        .def("__repr__", [](const gio::block_info& b) {
            return "BlockInfo(count=" + std::to_string(b.count) +
                   ", sizes=[" + std::to_string(b.sizes.size()) + " items])";
        });

    py::class_<gio::thick_info>(m, "ThickInfo", R"pbdoc(
        BED thick-drawing range (thickStart, thickEnd).
    )pbdoc")
        .def(py::init<>())
        .def(py::init<uint64_t, uint64_t>(), py::arg("start"), py::arg("end"))
        .def_readwrite("start", &gio::thick_info::start, "thickStart")
        .def_readwrite("end", &gio::thick_info::end, "thickEnd")
        .def("__repr__", [](const gio::thick_info& t) {
            return "ThickInfo(start=" + std::to_string(t.start) +
                   ", end=" + std::to_string(t.end) + ")";
        });

    py::class_<gio::rgb_color>(m, "RgbColor", R"pbdoc(
        BED itemRgb display color. Each channel is an int in [0, 255].
    )pbdoc")
        .def(py::init<>())
        .def(py::init<uint8_t, uint8_t, uint8_t>(),
             py::arg("red"), py::arg("green"), py::arg("blue"))
        .def_readwrite("red", &gio::rgb_color::red, "Red channel (0-255)")
        .def_readwrite("green", &gio::rgb_color::green, "Green channel (0-255)")
        .def_readwrite("blue", &gio::rgb_color::blue, "Blue channel (0-255)")
        .def("__repr__", [](const gio::rgb_color& c) {
            return "RgbColor(" + std::to_string(c.red) + ", " +
                   std::to_string(c.green) + ", " + std::to_string(c.blue) + ")";
        });

    py::class_<gio::bed_entry>(m, "BedEntry", R"pbdoc(
        A single BED record.

        Coordinates are BED-native: 0-based, half-open [start, end). This is
        the raw record's own coordinate system and is distinct from the closed
        [start, end] of Interval used as the grove key. Optional fields are
        None when absent.

        Parameters
        ----------
        chrom : str
        start : int
            0-based start (BED chromStart).
        end : int
            0-based exclusive end (BED chromEnd).
    )pbdoc")
        .def(py::init<>())
        .def(py::init<std::string, size_t, size_t>(),
             py::arg("chrom"), py::arg("start"), py::arg("end"))
        .def_readwrite("chrom", &gio::bed_entry::chrom, "Chromosome name")
        .def_readwrite("start", &gio::bed_entry::start,
                       "0-based start position (BED chromStart)")
        .def_readwrite("end", &gio::bed_entry::end,
                       "0-based exclusive end position (BED chromEnd)")
        .def_readwrite("name", &gio::bed_entry::name,
                       "Optional[str] feature name (BED4+)")
        .def_readwrite("score", &gio::bed_entry::score,
                       "Optional[int] score (BED5+)")
        .def_readwrite("strand", &gio::bed_entry::strand,
                       "Optional[str] strand — a single character ('+', '-', "
                       "'.'); assigning an empty or multi-character string "
                       "raises ValueError, None clears it (BED6+)")
        .def_readwrite("thickness", &gio::bed_entry::thickness,
                       "Optional[ThickInfo] thick range (BED7+)")
        .def_readwrite("item_rgb", &gio::bed_entry::item_rgb,
                       "Optional[RgbColor] display color (BED9+)")
        .def_readwrite("blocks", &gio::bed_entry::blocks,
                       "Optional[BlockInfo] block structure (BED12)")
        .def("__repr__", [](const gio::bed_entry& e) {
            return "BedEntry(chrom='" + e.chrom + "', start=" +
                   std::to_string(e.start) + ", end=" +
                   std::to_string(e.end) + ")";
        });
}

inline void bind_bed_reader(py::module_& m) {
    // BedEntry must already be registered (bind_bed_entry) — BedReader yields it.
    py::class_<gio::bed_reader>(m, "BedReader", R"pbdoc(
        A single-pass iterator over the records of a BED file.

        Iterate it directly to get BedEntry objects::

            for entry in pygenogrove.BedReader("peaks.bed"):
                ...

        Plain and BGZF/gzip-compressed (`.gz`) files are both accepted (the
        format is auto-detected). The reader owns an htslib file handle and is
        single-pass — it cannot be restarted or iterated twice.

        Parameters
        ----------
        path : str
            Path to the BED file. A missing/unreadable file raises an exception.
        skip_invalid_lines : bool, optional
            If False (default), a malformed line raises RuntimeError mid-iteration.
            If True, malformed lines are skipped silently.
    )pbdoc")
        .def(py::init([](const std::string& path, bool skip_invalid_lines) {
                 gio::bed_reader_options opts;
                 opts.skip_invalid_lines = skip_invalid_lines;
                 return std::make_unique<gio::bed_reader>(path, opts);
             }),
             py::arg("path"), py::arg("skip_invalid_lines") = false)
        .def("__iter__", [](gio::bed_reader& r) -> gio::bed_reader& { return r; })
        .def("__next__", [](gio::bed_reader& r) {
            gio::bed_entry entry;
            if (!r.read_next(entry)) {
                throw py::stop_iteration();
            }
            return entry;
        })
        .def("get_error_message", &gio::bed_reader::get_error_message,
             "Error message from the most recent read; empty on clean EOF.")
        .def("get_current_line", &gio::bed_reader::get_current_line,
             "1-based physical line number consumed so far (comments/blanks "
             "count); 0 before the first read.");
}