/*
 * Binding for gdt::key<interval, DataT> — a key wrapping an interval plus an
 * optional associated-data payload. Mirrors genogrove data_type/key.hpp,
 * fixed to the `interval` key type (named interval_key so other key types —
 * numeric, genomic_coordinate, kmer — can get their own bindings later; see
 * issue #1).
 *
 * One template covers both the dataless key (DataT = void) and data-carrying
 * keys; the `.data` accessor is only added when DataT is non-void.
 */
#pragma once

#include <pybind11/pybind11.h>

#include <type_traits>

#include <genogrove/data_type/interval.hpp>
#include <genogrove/data_type/key.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

template <typename DataT>
void bind_interval_key(py::module_& m, const char* name) {
    using key_t = gdt::key<gdt::interval, DataT>;

    auto cls = py::class_<key_t>(m, name, R"pbdoc(
        A key wrapping an interval stored in the grove structure.

        Returned by Grove.insert() and yielded by QueryResult iteration. Wraps a
        pointer into the grove's storage, so the key remains valid only as long
        as the originating Grove is alive (the key keeps the Grove alive).
    )pbdoc")
        .def_property_readonly(
            "value",
            [](const key_t& k) { return k.get_value(); },
            "The interval value of this key (returned by value/copy, so mutating "
            "it cannot corrupt the grove's B+ tree ordering)")
        .def("__str__", [](const key_t& k) { return k.to_string(); });

    if constexpr (!std::is_void_v<DataT>) {
        cls.def_property_readonly(
            "data",
            [](key_t& k) -> DataT& { return k.get_data(); },
            py::return_value_policy::reference_internal,
            "The associated data payload — a live, mutable reference into grove "
            "storage. Unlike .value, the data is not part of the B+ tree "
            "ordering, so mutating it in place is safe.");
    }
}