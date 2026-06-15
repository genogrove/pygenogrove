/*
 * Binding for gdt::numeric — a simple integer point key type.
 * Mirrors genogrove/include/genogrove/data_type/numeric.hpp.
 *
 * Point semantics (not a range): overlap is exact equality, aggregate is max.
 * Useful for non-genomic B+ tree use (ids, timestamps, …). Pairs with
 * NumericGrove / NumericKey / NumericQueryResult.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <string>

#include <genogrove/data_type/numeric.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

inline void bind_numeric(py::module_& m) {
    py::class_<gdt::numeric>(m, "Numeric", R"pbdoc(
        A simple integer point key: wraps a single int value (not a range).

        Overlap is exact equality (numeric(5) overlaps numeric(5) only), so a
        NumericGrove behaves as a B+ tree for point lookups. Ordering is by the
        integer value.

        Parameters
        ----------
        value : int
            The wrapped integer value. Defaults to INT_MIN (the aggregation
            sentinel) when omitted.
    )pbdoc")
        .def(py::init<>())
        .def(py::init<int>(), py::arg("value"))
        .def_property_readonly("value", &gdt::numeric::get_value,
                               "The wrapped integer value (read-only — see set_value)")
        .def("set_value", &gdt::numeric::set_value, py::arg("value"),
             R"pbdoc(
                 Set the integer value.

                 Do NOT call this on a Numeric already inserted into a Grove —
                 the value participates in B+ tree ordering, so mutating a stored
                 key silently corrupts it. Use only on values not yet inserted
                 (e.g. queries you intend to reuse).
             )pbdoc")
        .def("__str__", &gdt::numeric::to_string)
        .def("__repr__", [](const gdt::numeric& n) {
            return "Numeric(" + std::to_string(n.get_value()) + ")";
        })
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def(py::self == py::self)
        .def_static("overlaps", &gdt::numeric::overlaps,
                    py::arg("a"), py::arg("b"),
                    "Check if two Numerics overlap — true iff they are equal.");
}