/*
 * Binding for ggs::grove<KeyT, DataT> — the B+ tree container. Mirrors genogrove
 * structure/grove/grove.hpp. Generic over the key type KeyT (interval,
 * genomic_coordinate, …): instantiated per concrete key type from bindings.cpp,
 * producing a distinct Python class each time (Grove, GenomicCoordinateGrove, …).
 *
 * A single template covers the dataless grove (DataT = void, e.g. Grove) and
 * data-carrying groves (e.g. DataT = bed_entry → BedGrove). The differences are
 * the insert / add_external_key signatures (which gain a data argument), the
 * Key.data accessor, and the interval-only entry-deriving insert overloads;
 * these are switched with `if constexpr`.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <genogrove/data_type/interval.hpp>
#include <genogrove/data_type/key.hpp>
#include <genogrove/structure/grove/grove.hpp>

#include "../data_type/key.hpp"
#include "../data_type/query_result.hpp"
#include "../data_type/flanking_query_result.hpp"
#include "../io/entry_interval.hpp"

namespace py = pybind11;
namespace ggs = genogrove::structure;
namespace gdt = genogrove::data_type;

template <typename KeyT, typename DataT>
void bind_grove(py::module_& m, const char* grove_name,
                const char* key_name, const char* qr_name,
                const char* fr_name) {
    using grove_t = ggs::grove<KeyT, DataT>;
    using key_t = gdt::key<KeyT, DataT>;

    // The grove's Key, QueryResult and FlankingResult instantiations must be
    // registered first, since insert()/intersect()/flanking() use them.
    bind_key<KeyT, DataT>(m, key_name);
    bind_query_result<KeyT, DataT>(m, qr_name);
    bind_flanking_query_result<KeyT, DataT>(m, fr_name);

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
                   const KeyT& interval) {
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
                   const KeyT& interval, DataT data) {
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

        // ---- Sorted / bulk insertion (fast paths; non-void data only) ----
        cls.def("insert_sorted",
                [](grove_t& g, const std::string& index,
                   const KeyT& interval, DataT data) {
                    return g.insert_data(index, interval, std::move(data),
                                         ggs::sorted);
                },
                py::arg("index"), py::arg("interval"), py::arg("data"),
                py::return_value_policy::reference_internal,
                R"pbdoc(
                    Insert one (interval, data) record on the optimized sorted
                    path (rightmost-append, no tree traversal).

                    PRECONDITION: the interval must be greater than every interval
                    already present in this index. Insert in ascending order.
                    Violating this corrupts B+ tree ordering (queries silently
                    return wrong results). Use plain insert() if unsure.
                )pbdoc");

        cls.def("insert_bulk",
                [](grove_t& g, const std::string& index,
                   std::vector<std::pair<KeyT, DataT>> items,
                   bool presorted) {
                    if (presorted) {
                        return g.insert_data(index, items, ggs::sorted, ggs::bulk);
                    }
                    return g.insert_data(index, std::move(items), ggs::bulk);
                },
                py::arg("index"), py::arg("items"), py::arg("presorted") = false,
                py::return_value_policy::reference_internal,
                R"pbdoc(
                    Bulk-insert many (interval, data) records at once. 10-20x
                    faster than repeated insert() for large datasets — an empty
                    index is built bottom-up in O(n); a non-empty index appends.

                    Parameters
                    ----------
                    index : str
                        The index name (e.g., chromosome name).
                    items : list[tuple[Interval, object]]
                        The (interval, data) records to insert.
                    presorted : bool, optional
                        If False (default) the records are sorted by interval
                        first. If True, they are assumed already sorted ascending
                        (skips the sort — fastest).

                    Returns
                    -------
                    list[Key]
                        Stable key handles, in insertion (sorted) order.

                    PRECONDITION (appending to a non-empty index): every inserted
                    interval must be greater than every existing interval in the
                    index. Violating this corrupts B+ tree ordering.
                )pbdoc");

        // ---- Entry-deriving OVERLOADS of insert / insert_bulk: pass a file
        //      entry (or a list of them) and the Interval key is derived from
        //      the entry's native coordinates, so you never hand-convert
        //      (BED 0-based half-open, GFF 1-based inclusive). pybind resolves
        //      these against the explicit (interval, data) forms by signature.
        //      Only for the interval key type (the conversion yields an interval)
        //      and entry data types with a known conversion. ----
        if constexpr (std::is_same_v<KeyT, gdt::interval> &&
                      has_entry_interval<DataT>) {
            cls.def("insert",
                    [](grove_t& g, const std::string& index, DataT entry) {
                        KeyT iv = interval_from_entry(entry);
                        return g.insert_data(index, iv, std::move(entry));
                    },
                    py::arg("index"), py::arg("entry"),
                    py::return_value_policy::reference_internal,
                    R"pbdoc(
                        insert(index, entry) -> Key

                        Overload that takes a single file entry and derives the
                        Interval key from its native coordinates (BED half-open
                        [s, e) -> [s, e-1]; GFF 1-based [s, e] -> [s-1, e-1]).
                        The entry keeps its native coordinates as the payload.
                    )pbdoc");

            cls.def("insert_bulk",
                    [](grove_t& g, const std::string& index,
                       std::vector<DataT> entries, bool presorted) {
                        std::vector<std::pair<KeyT, DataT>> items;
                        items.reserve(entries.size());
                        for (auto& entry : entries) {
                            KeyT iv = interval_from_entry(entry);
                            items.emplace_back(iv, std::move(entry));
                        }
                        if (presorted) {
                            return g.insert_data(index, items, ggs::sorted,
                                                 ggs::bulk);
                        }
                        return g.insert_data(index, std::move(items), ggs::bulk);
                    },
                    py::arg("index"), py::arg("entries"),
                    py::arg("presorted") = false,
                    py::return_value_policy::reference_internal,
                    R"pbdoc(
                        insert_bulk(index, entries, presorted=False) -> list[Key]

                        Overload that takes a list of bare file entries (instead
                        of (Interval, data) tuples) and derives each Interval key
                        from the entry's native coordinates. Same append
                        precondition as the explicit form.
                    )pbdoc");
        }
    }

    // ---- Queries (identical for both cases) ----
    // keep_alive<0, 1>: the returned QueryResult (and the keys it yields) hold
    // pointers into the grove's storage, so the grove must outlive the result.
    cls.def("intersect",
            py::overload_cast<const KeyT&>(&grove_t::intersect),
            py::arg("query"), py::keep_alive<0, 1>(),
            R"pbdoc(
                Find all intervals that overlap with the query across all indices.
            )pbdoc")
       .def("intersect",
            py::overload_cast<const KeyT&, std::string_view>(
                &grove_t::intersect),
            py::arg("query"), py::arg("index"), py::keep_alive<0, 1>(),
            R"pbdoc(
                Find all intervals that overlap with the query in a specific index.
            )pbdoc")

        // ---- Flanking (nearest non-overlapping neighbours) ----
        .def("flanking",
             [](const grove_t& g, const KeyT& query,
                const std::string& index) {
                 return g.flanking(query, index);
             },
             py::arg("query"), py::arg("index"), py::keep_alive<0, 1>(),
             R"pbdoc(
                 Find the nearest non-overlapping keys on either side of the query
                 within an index (the predecessor and successor).

                 Returns a FlankingResult with `.predecessor` / `.successor`, each
                 a Key or None. Keys that overlap the query are excluded; for
                 nested intervals the predecessor is the one with the largest end
                 (smallest gap), not the sort-order maximum. Compute the gap
                 distance from the returned key, e.g.
                 `query.start - result.predecessor.value.end - 1`.
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
                [](grove_t& g, const KeyT& interval) {
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
                [](grove_t& g, const KeyT& interval, DataT data) {
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