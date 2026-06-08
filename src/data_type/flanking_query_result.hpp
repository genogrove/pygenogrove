/*
 * Binding for gdt::flanking_query_result<interval, DataT> — the result of
 * Grove.flanking(). Mirrors genogrove data_type/flanking_query_result.hpp.
 *
 * Holds the nearest non-overlapping key on each side of a query (predecessor /
 * successor); either may be None. The keys point into grove storage.
 */
#pragma once

#include <pybind11/pybind11.h>

#include <genogrove/data_type/interval.hpp>
#include <genogrove/data_type/flanking_query_result.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

template <typename DataT>
void bind_flanking_query_result(py::module_& m, const char* name) {
    using fqr_t = gdt::flanking_query_result<gdt::interval, DataT>;

    py::class_<fqr_t>(m, name, R"pbdoc(
        Result of a Grove.flanking() query: the nearest non-overlapping keys on
        either side of the query, in the grove's sort order.

        `predecessor` is the closest key entirely before the query (largest end
        with no overlap); `successor` is the closest key entirely after it
        (smallest start with no overlap). Either is None if no such key exists.
        Distance is computed by the caller from the key values, e.g. for a query
        Q and predecessor P: ``Q.start - P.value.end - 1`` (closed coordinates).
    )pbdoc")
        .def_property_readonly(
            "predecessor",
            [](const fqr_t& r) { return r.get_predecessor(); },
            py::return_value_policy::reference_internal,
            "Nearest non-overlapping key before the query (a Key), or None.")
        .def_property_readonly(
            "successor",
            [](const fqr_t& r) { return r.get_successor(); },
            py::return_value_policy::reference_internal,
            "Nearest non-overlapping key after the query (a Key), or None.");
}