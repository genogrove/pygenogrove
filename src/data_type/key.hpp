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

#include <pybind11/operators.h>
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
        .def("__str__", [](const key_t& k) { return k.to_string(); })
        .def("__repr__", [name](const key_t& k) {
            return std::string(name) + "(" + k.to_string() + ")";
        })
        // Comparisons delegate to the wrapped value (ignoring .data), mirroring
        // key_t's own operator==/</> semantics — matches the B+ tree's notion
        // of identity.
        .def(py::self == py::self)
        .def(py::self < py::self)
        .def(py::self > py::self)
        // Round-trips through Python to reuse the value type's own __hash__
        // rather than re-deriving each type's mixing logic here. A measurable
        // cost on bulk hashing (e.g. set(query_result)) would be the signal to
        // switch to a direct C++ call instead.
        .def("__hash__", [](const key_t& k) { return py::hash(py::cast(k.get_value())); });

    // Every grove carries a payload (the universal Grove's is JSON; BedKey/GffKey
    // carry a typed record), so .data is always present.
    cls.def_property_readonly(
        "data",
        [](key_t& k) -> DataT& { return k.get_data(); },
        // reference_internal only bites for a pybind11-wrapped return value
        // (true for BedKey/GffKey's live reference). json_value's custom caster
        // decodes straight to a native Python object, so it's a no-op there —
        // see the docstring below.
        py::return_value_policy::reference_internal,
        "The associated data payload (not part of B+ tree ordering). On the typed "
        "BedKey/GffKey it is a live, mutable reference into grove storage "
        "(mutating it in place is safe). On the universal Grove the payload is "
        "JSON, so .data returns a freshly decoded copy each access — mutating that "
        "copy does not persist; re-insert to change it.");
}