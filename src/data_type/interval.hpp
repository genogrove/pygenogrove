/*
 * Binding for gdt::interval — the genomic interval key type.
 * Mirrors genogrove/include/genogrove/data_type/interval.hpp.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <string>

#include <genogrove/data_type/interval.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

inline void bind_interval(py::module_& m) {
    py::class_<gdt::interval>(m, "Interval", R"pbdoc(
        A genomic interval with closed [start, end] coordinates (0-based, both inclusive).

        Parameters
        ----------
        start : int
            Start position of the interval (0-based, inclusive)
        end : int
            End position of the interval (0-based, inclusive)
    )pbdoc")
        .def(py::init<>())
        .def(py::init<size_t, size_t>(),
             py::arg("start"),
             py::arg("end"))
        .def_property_readonly("start",
                     &gdt::interval::get_start,
                     "Start position of the interval (read-only — see set_range)")
        .def_property_readonly("end",
                     &gdt::interval::get_end,
                     "End position of the interval (read-only — see set_range)")
        .def("set_range", &gdt::interval::set_range,
             py::arg("start"), py::arg("end"),
             R"pbdoc(
                 Atomically set both endpoints.

                 Do NOT call this on an interval that has already been inserted into
                 a Grove — mutating a stored key silently corrupts B+ tree ordering.
                 Use this only on intervals not yet inserted (e.g. queries you intend
                 to reuse).
             )pbdoc")
        .def("__str__", &gdt::interval::to_string)
        .def("__repr__", [](const gdt::interval& i) {
            return "Interval(" + std::to_string(i.get_start()) + ", " +
                   std::to_string(i.get_end()) + ")";
        })
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def(py::self == py::self)
        .def_static("overlaps", &gdt::interval::overlaps,
                   py::arg("a"),
                   py::arg("b"),
                   "Check if two intervals overlap (closed-interval semantics)");
}