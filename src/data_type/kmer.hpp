/*
 * Binding for gdt::kmer — a DNA k-mer key type (2-bit encoded, k <= 32).
 * Mirrors genogrove/include/genogrove/data_type/kmer.hpp.
 *
 * Point semantics like numeric: overlap is exact equality (same sequence AND
 * same length), aggregate is max-by-encoding. Only canonical bases A/C/G/T are
 * supported (case-insensitive). Pairs with KmerGrove / KmerKey / KmerQueryResult.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <cstdint>
#include <functional>
#include <string>

#include <genogrove/data_type/kmer.hpp>

namespace py = pybind11;
namespace gdt = genogrove::data_type;

inline void bind_kmer(py::module_& m) {
    py::class_<gdt::kmer> cls(m, "Kmer", R"pbdoc(
        A DNA k-mer: a length-k sequence over {A, C, G, T}, stored as a compact
        2-bit encoding (so k <= 32). Immutable.

        Overlap is exact equality — same bases AND same length — so a KmerGrove
        behaves as a k-mer dictionary for membership lookups. Ordering is by
        length first, then by the 2-bit encoding (lexicographic A < C < G < T).

        Parameters
        ----------
        sequence : str
            A DNA sequence over A/C/G/T (case-insensitive), length 0..32.
    )pbdoc");

    cls.def(py::init<>())
        .def(py::init([](const std::string& sequence) { return gdt::kmer(sequence); }),
             py::arg("sequence"),
             "Build a k-mer from a DNA sequence (A/C/G/T, case-insensitive). "
             "Raises ValueError on an invalid base or length > 32.")
        .def(py::init<uint64_t, uint8_t>(), py::arg("encoding"), py::arg("k"),
             "Build a k-mer from a precomputed 2-bit encoding and length k "
             "(1..32). Raises ValueError if k > 32.")
        .def_property_readonly("encoding", &gdt::kmer::get_encoding,
                               "The 2-bit encoding as a 64-bit integer.")
        .def_property_readonly("k", &gdt::kmer::get_k, "The k-mer length.")
        .def("__len__", &gdt::kmer::get_k)
        .def("__str__", &gdt::kmer::to_string)
        .def("__repr__", [](const gdt::kmer& km) {
            return "Kmer('" + km.to_string() + "')";
        })
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def(py::self == py::self)
        .def("__hash__", [](const gdt::kmer& km) {
            // == compares both encoding AND k, so the hash must mix both.
            std::size_t h = std::hash<uint64_t>{}(km.get_encoding());
            h ^= std::hash<uint8_t>{}(km.get_k()) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        })
        .def_static("overlaps", &gdt::kmer::overlaps,
                    py::arg("a"), py::arg("b"),
                    "Check if two k-mers overlap — true iff identical (same bases "
                    "and same length).")
        .def_static("is_valid", &gdt::kmer::is_valid, py::arg("sequence"),
                    "Whether a sequence contains only A/C/G/T (case-insensitive).");

    cls.attr("max_k") = gdt::kmer::max_k;  // 32
}