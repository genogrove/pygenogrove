/*
 * Binding for gdt::genomic_coordinate — the stranded genomic interval key type.
 * Mirrors genogrove/include/genogrove/data_type/genomic_coordinate.hpp.
 *
 * Like Interval but strand-aware: overlap requires both coordinate overlap AND
 * strand compatibility ('*' is a wildcard matching any strand). Sorting is
 * coordinate-first (start, then end, then strand with order * < . < + < -).
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <string>

#include <genogrove/data_type/genomic_coordinate.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

inline void bind_genomic_coordinate(py::module_& m) {
    py::class_<gdt::genomic_coordinate>(m, "GenomicCoordinate", R"pbdoc(
        A stranded genomic interval: closed [start, end] coordinates (0-based,
        both inclusive) plus a strand.

        Strand is one of:
        - '+' : forward / plus strand
        - '-' : reverse / minus strand
        - '.' : no strand information (strand-agnostic)
        - '*' : wildcard — matches any strand in overlap queries

        Overlap requires BOTH coordinate overlap AND strand compatibility:
        equal strands overlap, and '*' matches any strand. Sorting is
        coordinate-first (start, then end, then strand).

        Parameters
        ----------
        strand : str
            One of '+', '-', '.', '*' (a single character).
        start : int
            Start position (0-based, inclusive).
        end : int
            End position (0-based, inclusive).
    )pbdoc")
        .def(py::init<>())
        .def(py::init<char, size_t, size_t>(),
             py::arg("strand"),
             py::arg("start"),
             py::arg("end"))
        .def_property_readonly("strand",
                     &gdt::genomic_coordinate::get_strand,
                     "Strand character (read-only — see set_strand)")
        .def_property_readonly("start",
                     &gdt::genomic_coordinate::get_start,
                     "Start position (read-only — see set_range)")
        .def_property_readonly("end",
                     &gdt::genomic_coordinate::get_end,
                     "End position (read-only — see set_range)")
        .def("set_range", &gdt::genomic_coordinate::set_range,
             py::arg("start"), py::arg("end"),
             R"pbdoc(
                 Atomically set both endpoints.

                 Do NOT call this on a coordinate that has already been inserted
                 into a Grove — mutating a stored key silently corrupts B+ tree
                 ordering. Use this only on coordinates not yet inserted (e.g.
                 queries you intend to reuse).
             )pbdoc")
        .def("set_strand", &gdt::genomic_coordinate::set_strand,
             py::arg("strand"),
             R"pbdoc(
                 Set the strand ('+', '-', '.', '*').

                 Same warning as set_range: do NOT call this on a coordinate
                 already inserted into a Grove — strand participates in B+ tree
                 ordering, so mutating a stored key corrupts it silently.
             )pbdoc")
        .def("__str__", &gdt::genomic_coordinate::to_string)
        .def("__repr__", [](const gdt::genomic_coordinate& c) {
            return std::string("GenomicCoordinate('") + c.get_strand() + "', " +
                   std::to_string(c.get_start()) + ", " +
                   std::to_string(c.get_end()) + ")";
        })
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def(py::self == py::self)
        .def_static("overlaps", &gdt::genomic_coordinate::overlaps,
                   py::arg("a"),
                   py::arg("b"),
                   R"pbdoc(
                       Check if two coordinates overlap: coordinates must
                       intersect AND strands must be compatible (equal, or one is
                       the wildcard '*').
                   )pbdoc");
}