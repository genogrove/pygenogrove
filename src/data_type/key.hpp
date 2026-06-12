/*
 * Binding for gdt::key<KeyT, DataT> — a key wrapping a key value (interval,
 * genomic_coordinate, …) plus an optional associated-data payload. Mirrors
 * genogrove data_type/key.hpp.
 *
 * Generic over the key type KeyT: instantiated per concrete key type (and per
 * DataT) from the grove binding, producing a distinct Python class each time
 * (e.g. Key, BedKey, GenomicCoordinateKey). One template covers both the
 * dataless key (DataT = void) and data-carrying keys; the `.data` accessor is
 * only added when DataT is non-void.
 */
#pragma once

#include <pybind11/pybind11.h>

#include <type_traits>

#include <genogrove/data_type/key.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

template <typename KeyT, typename DataT>
void bind_key(py::module_& m, const char* name) {
    using key_t = gdt::key<KeyT, DataT>;

    auto cls = py::class_<key_t>(m, name, R"pbdoc(
        A key wrapping a key value stored in the grove structure.

        Returned by Grove.insert() and yielded by QueryResult iteration. Wraps a
        pointer into the grove's storage, so the key remains valid only as long
        as the originating Grove is alive (the key keeps the Grove alive).
    )pbdoc")
        .def_property_readonly(
            "value",
            [](const key_t& k) { return k.get_value(); },
            "The key value (returned by value/copy, so mutating it cannot "
            "corrupt the grove's B+ tree ordering)")
        .def("__str__", [](const key_t& k) { return k.to_string(); });

    // Every grove carries a payload (the universal Grove's is JSON; BedKey/GffKey
    // carry a typed record), so .data is always present.
    cls.def_property_readonly(
        "data",
        [](key_t& k) -> DataT& { return k.get_data(); },
        py::return_value_policy::reference_internal,
        "The associated data payload (not part of B+ tree ordering). On the typed "
        "BedKey/GffKey it is a live, mutable reference into grove storage "
        "(mutating it in place is safe). On the universal Grove the payload is "
        "JSON, so .data returns a freshly decoded copy each access — mutating that "
        "copy does not persist; re-insert to change it.");
}