/*
 * Binding for ggs::grove<interval, DataT> — the B+ tree container.
 * Mirrors genogrove structure/grove/grove.hpp, fixed to the `interval` key
 * type (named interval_grove so other key types — numeric, genomic_coordinate,
 * kmer — can get their own grove bindings later; see issue #1).
 *
 * A single template covers both the dataless grove (DataT = void, exposed as
 * Grove/Key/QueryResult) and data-carrying groves (e.g. DataT = bed_entry,
 * exposed as BedGrove/BedKey/BedQueryResult). The only differences between the
 * two are the insert / add_external_key signatures (which gain a data argument)
 * and the Key.data accessor; these are switched with `if constexpr`.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

#include <genogrove/data_type/interval.hpp>
#include <genogrove/data_type/key.hpp>
#include <genogrove/structure/grove/grove.hpp>

#include "../data_type/interval_key.hpp"
#include "../data_type/query_result.hpp"

namespace py = pybind11;
namespace ggs = genogrove::structure;
namespace gdt = genogrove::data_type;

template <typename DataT>
void bind_interval_grove(py::module_& m, const char* grove_name,
                         const char* key_name, const char* qr_name) {
    using grove_t = ggs::grove<gdt::interval, DataT>;
    using key_t = gdt::key<gdt::interval, DataT>;

    // The grove's Key and QueryResult instantiations must be registered first,
    // since insert()/intersect() take and return them.
    bind_interval_key<DataT>(m, key_name);
    bind_query_result<DataT>(m, qr_name);

    auto cls = py::class_<grove_t>(m, grove_name, R"pbdoc(
        A B+ tree container for efficient genomic interval storage and querying.

        The grove supports multi-index operations, where each index (e.g., chromosome)
        maintains its own B+ tree structure.

        Parameters
        ----------
        order : int, optional
            Maximum branching factor of the B+ tree (default: 3, minimum: 3).
            Controls the maximum number of keys per node (order - 1).
    )pbdoc")
        .def(py::init<>())
        .def(py::init<int>(), py::arg("order"))
        .def("__str__", [grove_name](const grove_t& g) {
            return std::string(grove_name) + "(size=" +
                   std::to_string(g.indexed_vertex_count()) + ")";
        })
        .def("__repr__", [grove_name](const grove_t& g) {
            return std::string(grove_name) + "(order=" +
                   std::to_string(g.get_order()) +
                   ", size=" + std::to_string(g.indexed_vertex_count()) + ")";
        })
        .def("__len__", &grove_t::indexed_vertex_count)
        .def("size", &grove_t::indexed_vertex_count,
             "Number of indexed intervals across all indices (alias of len)")
        .def("indexed_vertex_count", &grove_t::indexed_vertex_count,
             "Number of indexed intervals (B+ tree leaf keys)")
        .def("get_order", &grove_t::get_order,
             "Get the order (branching factor) of the B+ tree");

    // ---- Insert (data argument only present when DataT is non-void) ----
    if constexpr (std::is_void_v<DataT>) {
        cls.def("insert",
                [](grove_t& g, const std::string& index,
                   const gdt::interval& interval) {
                    key_t key(interval);
                    return g.insert(index, key);
                },
                py::arg("index"), py::arg("interval"),
                py::return_value_policy::reference_internal,
                R"pbdoc(
                    Insert an interval into the grove at the specified index.

                    Parameters
                    ----------
                    index : str
                        The index name (e.g., chromosome name like "chr1")
                    interval : Interval
                        The interval to insert (copied into the grove)

                    Returns
                    -------
                    Key
                        Stable reference to the inserted key. Remains valid as
                        long as the Grove is alive.
                )pbdoc");
    } else {
        cls.def("insert",
                [](grove_t& g, const std::string& index,
                   const gdt::interval& interval, DataT data) {
                    return g.insert_data(index, interval, std::move(data));
                },
                py::arg("index"), py::arg("interval"), py::arg("data"),
                py::return_value_policy::reference_internal,
                R"pbdoc(
                    Insert an interval with associated data at the given index.

                    Parameters
                    ----------
                    index : str
                        The index name (e.g., chromosome name like "chr1")
                    interval : Interval
                        The interval key (copied into the grove). Drives B+ tree
                        ordering — do not mutate it after insertion.
                    data : object
                        The associated data payload (copied into the grove).

                    Returns
                    -------
                    Key
                        Stable reference to the inserted key. Its .data payload
                        is freely mutable; its .value (interval) is not.
                )pbdoc");
    }

    // ---- Queries (identical for both cases) ----
    // keep_alive<0, 1>: the returned QueryResult (and the keys it yields) hold
    // pointers into the grove's storage, so the grove must outlive the result.
    cls.def("intersect",
            py::overload_cast<const gdt::interval&>(&grove_t::intersect),
            py::arg("query"), py::keep_alive<0, 1>(),
            R"pbdoc(
                Find all intervals that overlap with the query across all indices.
            )pbdoc")
       .def("intersect",
            py::overload_cast<const gdt::interval&, std::string_view>(
                &grove_t::intersect),
            py::arg("query"), py::arg("index"), py::keep_alive<0, 1>(),
            R"pbdoc(
                Find all intervals that overlap with the query in a specific index.
            )pbdoc")

        // ---- Graph overlay (directed edges between keys) ----
        .def("add_edge",
             [](grove_t& g, key_t* source, key_t* target) {
                 g.add_edge(source, target);
             },
             py::arg("source"), py::arg("target"),
             R"pbdoc(
                 Add a directed edge from source to target.

                 source and target must be Keys belonging to this Grove (returned
                 by insert(), add_external_key(), or yielded by a QueryResult).
                 Raises ValueError if either is None.
             )pbdoc")
        .def("remove_edge",
             [](grove_t& g, key_t* source, key_t* target) {
                 return g.remove_edge(source, target);
             },
             py::arg("source"), py::arg("target"),
             "Remove the directed edge from source to target. Returns True if an "
             "edge was removed, False if it did not exist.")
        .def("has_edge",
             [](const grove_t& g, const key_t* source, const key_t* target) {
                 return g.has_edge(source, target);
             },
             py::arg("source"), py::arg("target"),
             "Return True if a directed edge from source to target exists.")
        .def("get_neighbors",
             [](grove_t& g, key_t* source) {
                 return g.get_neighbors(source);
             },
             py::arg("source"),
             py::return_value_policy::reference_internal,
             R"pbdoc(
                 Return the list of target Keys directly reachable from source.

                 The returned Keys point into this Grove's storage and remain valid
                 only while the Grove is alive.
             )pbdoc")
        .def("out_degree",
             [](const grove_t& g, const key_t* source) {
                 return g.out_degree(source);
             },
             py::arg("source"),
             "Number of outgoing edges from source.")
        .def("edge_count", &grove_t::edge_count,
             "Total number of directed edges in the graph overlay.")
        .def("vertex_count_with_edges", &grove_t::vertex_count_with_edges,
             "Number of keys that have at least one outgoing edge.");

    // ---- External (graph-only) key — gains a data argument when non-void ----
    if constexpr (std::is_void_v<DataT>) {
        cls.def("add_external_key",
                [](grove_t& g, const gdt::interval& interval) {
                    return g.add_external_key(interval);
                },
                py::arg("interval"),
                py::return_value_policy::reference_internal,
                R"pbdoc(
                    Add a key that lives outside the B+ tree index but can
                    participate in the graph overlay (e.g. an enhancer linked to
                    indexed exons).

                    The interval is copied into the Grove. Returns a stable Key
                    that remains valid as long as the Grove is alive. External
                    keys are not returned by intersect() queries.
                )pbdoc");
    } else {
        cls.def("add_external_key",
                [](grove_t& g, const gdt::interval& interval, DataT data) {
                    return g.add_external_key(interval, std::move(data));
                },
                py::arg("interval"), py::arg("data"),
                py::return_value_policy::reference_internal,
                R"pbdoc(
                    Add a key (interval + data) that lives outside the B+ tree
                    index but can participate in the graph overlay.

                    Both the interval and the data are copied into the Grove.
                    Returns a stable Key that remains valid as long as the Grove
                    is alive. External keys are not returned by intersect()
                    queries. (Note: this takes a data argument, unlike the
                    dataless Grove.add_external_key.)
                )pbdoc");
    }

    // ---- Serialization (zlib-compressed .gg binary) ----
    cls.def("serialize",
            [](const grove_t& g, const std::string& path) {
                std::ofstream os(path, std::ios::binary);
                if (!os) {
                    throw std::runtime_error(
                        "Failed to open file for writing: " + path);
                }
                g.serialize(os);
                if (!os) {
                    throw std::runtime_error(
                        "Failed to write grove to file: " + path);
                }
            },
            py::arg("path"),
            R"pbdoc(
                Serialize the Grove (intervals + associated data + graph overlay)
                to a zlib-compressed binary file at the given path.
            )pbdoc")
       .def_static("deserialize",
            [](const std::string& path) {
                std::ifstream is(path, std::ios::binary);
                if (!is) {
                    throw std::runtime_error(
                        "Failed to open file for reading: " + path);
                }
                return grove_t::deserialize(is);
            },
            py::arg("path"),
            R"pbdoc(
                Load a Grove previously written with serialize(). Returns a new
                Grove with the same intervals, associated data, and graph edges.
            )pbdoc");
}