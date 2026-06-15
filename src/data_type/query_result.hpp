/*
 * Binding for gdt::query_result<KeyT, DataT> — the container returned by
 * Grove.intersect(). Mirrors genogrove data_type/query_result.hpp. Generic over
 * the key type KeyT (instantiated per concrete key type from the grove binding).
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <genogrove/data_type/query_result.hpp>

#include "key_list.hpp"

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
        .def_property_readonly("keys",
                               [](py::object self) {
                                   // Pin each Key to this QueryResult (which keeps
                                   // its Grove alive) so an extracted Key can't
                                   // dangle after the list is dropped — issue #37.
                                   return pinned_key_list(
                                       self.cast<const qr_t&>().get_keys(), self);
                               },
                               "List of matching keys; each Key keeps this result "
                               "(and its Grove) alive.")
        .def("__len__", [](const qr_t& qr) { return qr.get_keys().size(); })
        .def("__iter__", [](const qr_t& qr) {
            return py::make_iterator(qr.get_keys().begin(), qr.get_keys().end());
        }, py::keep_alive<0, 1>());
}