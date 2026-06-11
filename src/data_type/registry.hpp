/*
 * Binding for gdt::registry<Key, Tag, Payload> — a process-wide singleton that
 * interns values into small stable integer ids. Mirrors genogrove
 * data_type/registry.hpp.
 *
 * Generic over the registry type: instantiated per concrete (Key[, Tag, Payload])
 * from bindings.cpp. The single-arg vs (key, payload) intern() form is switched
 * on registry::key_is_payload with `if constexpr`. The Tag template parameter is
 * a C++ compile-time discriminator for independent same-type pools; it has no
 * Python analogue, so Python sees one default-tag singleton per bound class.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // std::optional<id> -> int | None

#include <fstream>
#include <stdexcept>
#include <string>

#include <genogrove/data_type/registry.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

template <typename Key, typename Tag = void, typename Payload = Key>
void bind_registry(py::module_& m, const char* name) {
    using reg_t = gdt::registry<Key, Tag, Payload>;
    using id_type = typename reg_t::id_type;

    auto cls = py::class_<reg_t>(m, name, R"pbdoc(
        A process-wide singleton that interns values into small, stable integer
        ids (deduplicated). Access the single instance with instance().

        Global state — use reset() (or clear()) to wipe it, e.g. between tests.
    )pbdoc");

    cls.def_static("instance", &reg_t::instance,
                   py::return_value_policy::reference,
                   "Return the process-wide singleton instance.");

    if constexpr (reg_t::key_is_payload) {
        cls.def("intern",
                [](reg_t& r, const Key& value) { return r.intern(value); },
                py::arg("value"),
                R"pbdoc(
                    Intern a value and return its stable integer id. Idempotent:
                    the same value always returns the same id (deduplicated).
                )pbdoc");
    } else {
        cls.def("intern",
                [](reg_t& r, const Key& key, const Payload& payload) {
                    return r.intern(key, payload);
                },
                py::arg("key"), py::arg("payload"),
                R"pbdoc(
                    Intern key -> payload and return key's stable id. First write
                    wins: re-interning an existing key keeps its original payload.
                )pbdoc");
    }

    cls.def("find",
            [](const reg_t& r, const Key& key) { return r.find(key); },
            py::arg("key"),
            "Return the id for key if interned, otherwise None. Does not insert.")
       .def("get",
            [](const reg_t& r, id_type id) { return r.get(id); },
            py::arg("id"),
            "Return the payload for id. Raises IndexError if id is invalid.")
       .def("contains",
            [](const reg_t& r, id_type id) { return r.contains(id); },
            py::arg("id"),
            "Whether id refers to a valid entry.")
       .def("size", &reg_t::size, "Number of interned entries.")
       .def("__len__", &reg_t::size)
       .def("empty", &reg_t::empty, "Whether the registry has no entries.")
       .def("clear", &reg_t::clear,
            "Remove all interned data; ids restart from 0 afterward.")
       .def_static("reset", &reg_t::reset,
                   "Clear the singleton (convenience for e.g. test isolation; "
                   "equivalent to instance().clear()).");

    cls.def("serialize",
            [](const reg_t& r, const std::string& path) {
                std::ofstream os(path, std::ios::binary);
                if (!os) {
                    throw std::runtime_error(
                        "Failed to open file for writing: " + path);
                }
                r.serialize(os);
                if (!os) {
                    throw std::runtime_error(
                        "Failed to write registry to file: " + path);
                }
            },
            py::arg("path"),
            "Serialize the registry's (key, payload) entries to a binary file.")
       .def_static("deserialize",
            [](const std::string& path) -> reg_t& {
                std::ifstream is(path, std::ios::binary);
                if (!is) {
                    throw std::runtime_error(
                        "Failed to open file for reading: " + path);
                }
                return reg_t::deserialize(is);
            },
            py::arg("path"),
            py::return_value_policy::reference,
            R"pbdoc(
                Load entries from a file written with serialize() INTO the
                singleton (replacing its current data) and return it.
            )pbdoc");

    // Sentinel id (= max uint32) returned where "no id" is meaningful.
    cls.attr("null_id") = reg_t::null_id;
}