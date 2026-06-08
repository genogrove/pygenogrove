/*
 * Binding for gdt::query_result<KeyT, DataT> — the container returned by
 * Grove.intersect(). Mirrors genogrove data_type/query_result.hpp. Generic over
 * the key type KeyT (instantiated per concrete key type from the grove binding).
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <genogrove/data_type/query_result.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

template <typename KeyT, typename DataT>
void bind_query_result(py::module_& m, const char* name) {
    using qr_t = gdt::query_result<KeyT, DataT>;

    py::class_<qr_t>(m, name, R"pbdoc(
        Result of an intersect() query: the query interval plus the matching keys.
    )pbdoc")
        .def_property_readonly("query", &qr_t::get_query,
                               "The query interval used for this search")
        .def_property_readonly("keys", &qr_t::get_keys,
                               py::return_value_policy::reference_internal,
                               "List of matching keys (pointers into the grove)")
        .def("__len__", [](const qr_t& qr) { return qr.get_keys().size(); })
        .def("__iter__", [](const qr_t& qr) {
            return py::make_iterator(qr.get_keys().begin(), qr.get_keys().end());
        }, py::keep_alive<0, 1>());
}