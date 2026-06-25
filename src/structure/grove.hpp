/*
 * Binding for ggs::grove<KeyT, DataT, EdgeT> — the B+ tree container. Mirrors
 * genogrove structure/grove/grove.hpp. Generic over the key type KeyT, a
 * (non-void) data payload DataT, and an optional edge-metadata type EdgeT
 * (void by default): instantiated per concrete type tuple from bindings.cpp,
 * producing a distinct Python class each time (Grove =
 * grove<genomic_coordinate, json_value, json_value>, BedGrove =
 * grove<genomic_coordinate, bed_entry>, …).
 *
 * Every grove carries a payload. Three type-dependent variations are switched
 * with `if constexpr`: the insert/add_external_key `data` argument defaults to
 * None for the JSON payload (grove_data_optional); the entry-deriving
 * insert(index, entry) overloads exist only for the genomic_coordinate key with
 * a derivable entry type; and the labelled-edge methods (add_edge with a payload,
 * get_edges, get_neighbors_if, link_with) exist only when EdgeT is non-void.
 */
#pragma once

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <functional>
#include <fstream>
#include <optional>
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
#include "../data_type/key_list.hpp"
#include "../data_type/query_result.hpp"
#include "../data_type/flanking_query_result.hpp"
#include "../io/entry_interval.hpp"

namespace py = pybind11;
namespace ggs = genogrove::structure;
namespace gdt = genogrove::data_type;

// Customization point: data types whose insert()/add_external_key() `data`
// argument should default to an absent payload (so dataless inserts can omit
// it). Specialized to true for the JSON payload type in bindings.cpp; false for
// typed payloads (bed_entry/gff_entry), where a default makes no sense.
template <typename>
inline constexpr bool grove_data_optional = false;

template <typename KeyT, typename DataT, typename EdgeT = void>
void bind_grove(py::module_& m, const char* grove_name,
                const char* key_name, const char* qr_name,
                const char* fr_name) {
    using grove_t = ggs::grove<KeyT, DataT, EdgeT>;
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

    // ---- Insert (every grove carries a data payload) ----
    {
        auto insert_fn = [](grove_t& g, const std::string& index,
                            const KeyT& key, DataT data) {
            return g.insert_data(index, key, std::move(data));
        };
        const char* insert_doc = R"pbdoc(
                    Insert a key with an associated data payload at the given index.

                    Parameters
                    ----------
                    index : str
                        The index name (e.g., chromosome name like "chr1")
                    key : GenomicCoordinate
                        The key (copied into the grove). Drives B+ tree ordering —
                        do not mutate it after insertion.
                    data : object
                        The associated data payload (copied into the grove). On the
                        universal Grove this is any JSON-serializable value
                        (dict / list / scalar / None) and defaults to None.

                    Returns
                    -------
                    Key
                        Stable reference to the inserted key.
                )pbdoc";
        if constexpr (grove_data_optional<DataT>) {
            cls.def("insert", insert_fn,
                    py::arg("index"), py::arg("key"), py::arg("data") = DataT{},
                    py::return_value_policy::reference_internal, insert_doc);
        } else {
            cls.def("insert", insert_fn,
                    py::arg("index"), py::arg("key"), py::arg("data"),
                    py::return_value_policy::reference_internal, insert_doc);
        }

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
                [](py::object self, const std::string& index,
                   std::vector<std::pair<KeyT, DataT>> items,
                   bool presorted) {
                    auto& g = self.cast<grove_t&>();
                    // The bulk tree build is a long C++ loop touching no Python
                    // objects (items were already converted to C++); release the
                    // GIL for it, then build the result list with it reacquired.
                    std::vector<key_t*> keys;
                    {
                        py::gil_scoped_release rel;
                        keys = presorted
                                   ? g.insert_data(index, items, ggs::sorted, ggs::bulk)
                                   : g.insert_data(index, std::move(items), ggs::bulk);
                    }
                    // Pin each returned Key to the Grove so extracted keys can't
                    // dangle after the list is dropped — issue #37.
                    return pinned_key_list(keys, self);
                },
                py::arg("index"), py::arg("items"), py::arg("presorted") = false,
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
        //      Only for the genomic_coordinate key type (the conversion yields a
        //      stranded coordinate) and entry data types with a known conversion. ----
        if constexpr (std::is_same_v<KeyT, gdt::genomic_coordinate> &&
                      has_entry_coordinate<DataT>) {
            cls.def("insert",
                    [](grove_t& g, const std::string& index, DataT entry) {
                        KeyT k = genomic_coordinate_from_entry(entry);
                        return g.insert_data(index, k, std::move(entry));
                    },
                    py::arg("index"), py::arg("entry"),
                    py::return_value_policy::reference_internal,
                    R"pbdoc(
                        insert(index, entry) -> Key

                        Overload that takes a single file entry and derives the
                        GenomicCoordinate key from its native coordinates + strand
                        (BED half-open [s, e) -> [s, e-1]; GFF 1-based [s, e] ->
                        [s-1, e-1]; strand from the BED6/GFF strand column, or '.'
                        if absent). The entry keeps its native coordinates as the
                        payload.
                    )pbdoc");

            cls.def("insert_bulk",
                    [](py::object self, const std::string& index,
                       std::vector<DataT> entries, bool presorted) {
                        auto& g = self.cast<grove_t&>();
                        std::vector<std::pair<KeyT, DataT>> items;
                        items.reserve(entries.size());
                        for (auto& entry : entries) {
                            KeyT k = genomic_coordinate_from_entry(entry);
                            items.emplace_back(k, std::move(entry));
                        }
                        // Release the GIL around the pure-C++ bulk build, then
                        // build the result list with it reacquired.
                        std::vector<key_t*> keys;
                        {
                            py::gil_scoped_release rel;
                            keys = presorted
                                       ? g.insert_data(index, items, ggs::sorted,
                                                       ggs::bulk)
                                       : g.insert_data(index, std::move(items),
                                                       ggs::bulk);
                        }
                        // Pin each returned Key to the Grove — issue #37.
                        return pinned_key_list(keys, self);
                    },
                    py::arg("index"), py::arg("entries"),
                    py::arg("presorted") = false,
                    R"pbdoc(
                        insert_bulk(index, entries, presorted=False) -> list[Key]

                        Overload that takes a list of bare file entries (instead
                        of (GenomicCoordinate, data) tuples) and derives each
                        GenomicCoordinate key from the entry's native coordinates
                        + strand. Same append precondition as the explicit form.
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
        .def("flanking",
             [](const grove_t& g, const KeyT& query, const std::string& index,
                std::function<bool(const KeyT&, const KeyT&)> is_compatible) {
                 // The predicate calls back into Python, so the GIL must be held
                 // for the whole query — do NOT release it here.
                 return g.flanking(query, index, std::move(is_compatible));
             },
             py::arg("query"), py::arg("index"), py::arg("is_compatible"),
             py::keep_alive<0, 1>(),
             R"pbdoc(
                 flanking(query, index, is_compatible) -> FlankingResult

                 Predicate-filtered flanking: like flanking(query, index), but only
                 candidate keys for which `is_compatible(candidate, query)` returns
                 True are considered as neighbours. `candidate` and `query` are key
                 values (e.g. GenomicCoordinate); the predicate is applied at every
                 leaf candidate before the overlap/distance checks.

                 The canonical use is strand-aware neighbours on a
                 GenomicCoordinateGrove — the nearest non-overlapping key on the
                 same strand:

                     g.flanking(q, "chr1",
                                lambda cand, q: cand.strand == q.strand)

                 (Internal-node pruning ignores the predicate, so subtrees holding
                 only incompatible keys are still traversed and filtered at the
                 leaves — correct, just not pruned. Exceptions raised by the
                 predicate propagate out.)
             )pbdoc")

        // ---- Graph overlay (directed edges between keys) ----
        .def("add_edge",
             [](grove_t& g, key_t* source, key_t* target) {
                 g.add_edge(source, target);
             },
             py::arg("source").none(false), py::arg("target").none(false),
             R"pbdoc(
                 Add a directed edge from source to target.

                 source and target must be Keys belonging to this Grove (returned
                 by insert(), add_external_key(), or yielded by a QueryResult).
                 Raises TypeError if either is None.
             )pbdoc")
        .def("remove_edge",
             [](grove_t& g, key_t* source, key_t* target) {
                 return g.remove_edge(source, target);
             },
             py::arg("source").none(false), py::arg("target").none(false),
             "Remove the directed edge from source to target. Returns True if an "
             "edge was removed, False if it did not exist.")
        .def("has_edge",
             [](const grove_t& g, const key_t* source, const key_t* target) {
                 return g.has_edge(source, target);
             },
             py::arg("source").none(false), py::arg("target").none(false),
             "Return True if a directed edge from source to target exists.")
        .def("get_neighbors",
             [](py::object self, key_t* source) {
                 // Pin each Key to the Grove so an extracted neighbor can't
                 // dangle after the list is dropped — issue #37.
                 return pinned_key_list(
                     self.cast<grove_t&>().get_neighbors(source), self);
             },
             py::arg("source").none(false),
             R"pbdoc(
                 Return the list of target Keys directly reachable from source.

                 The returned Keys point into this Grove's storage and remain valid
                 only while the Grove is alive.
             )pbdoc")
        .def("out_degree",
             [](const grove_t& g, const key_t* source) {
                 return g.out_degree(source);
             },
             py::arg("source").none(false),
             "Number of outgoing edges from source.")
        .def("edge_count", &grove_t::edge_count,
             "Total number of directed edges in the graph overlay.")
        .def("vertex_count_with_edges", &grove_t::vertex_count_with_edges,
             "Number of keys that have at least one outgoing edge.")

        // ---- Graph edge removal / bulk linking ----
        .def("remove_edges_from",
             [](grove_t& g, key_t* source) {
                 return g.remove_edges_from(source);
             },
             py::arg("source").none(false),
             "Remove all outgoing edges from source. Returns the number removed.")
        .def("remove_edges_to",
             [](grove_t& g, key_t* target) {
                 return g.remove_edges_to(target);
             },
             py::arg("target").none(false),
             "Remove all incoming edges to target (O(E) scan over the graph). "
             "Returns the number removed.")
        .def("remove_all_edges",
             [](grove_t& g, key_t* key) { return g.remove_all_edges(key); },
             py::arg("key").none(false),
             "Remove every edge touching key, incoming and outgoing. Returns the "
             "total number removed.")
        .def("clear_graph", &grove_t::clear_graph,
             "Remove all edges from the graph overlay. The keys themselves are "
             "left intact.")
        .def("graph_empty", &grove_t::graph_empty,
             "Return True if the graph overlay holds no edges.")
        .def("link_if",
             [](grove_t& g, const std::vector<key_t*>& keys,
                std::function<bool(key_t*, key_t*)> predicate) {
                 // The predicate calls back into Python — keep the GIL held.
                 g.link_if(keys, std::move(predicate));
             },
             py::arg("keys"), py::arg("predicate"),
             R"pbdoc(
                 link_if(keys, predicate) -> None

                 Add an (unlabelled) directed edge between each adjacent pair
                 (keys[i], keys[i+1]) for which predicate(keys[i], keys[i+1])
                 returns True. `keys` is typically the list returned by a bulk
                 insert; the predicate receives two Keys. Use link_with() to attach
                 edge metadata.
             )pbdoc")

        // ---- Vertex / storage counts ----
        .def("vertex_count", &grove_t::vertex_count,
             "Total number of keys in the grove: indexed (B+ tree) plus external "
             "(graph-only) keys, including isolated ones with no edges.")
        .def("external_vertex_count", &grove_t::external_vertex_count,
             "Number of external (graph-only) keys — those added with "
             "add_external_key(), not indexed in any B+ tree.")
        .def("key_storage_size", &grove_t::key_storage_size,
             R"pbdoc(
                 Total slots in the indexed-key storage: live leaf data keys +
                 internal B+ tree separator keys + dead slots left behind by
                 remove_key(). Grows across insert/remove cycles until compact()
                 reclaims the dead slots — so it is >= indexed_vertex_count().
             )pbdoc");

    // ---- Predicate-filtered edge removal (every grove; #33) ----
    // genogrove's remove_edges_if takes a generic predicate over `const edge&`
    // ({ target, metadata }); we adapt it to a Python callable. The predicate
    // re-enters Python, so the GIL stays held (no call_guard). The Python-facing
    // signature differs by whether this grove's edges carry metadata.
    if constexpr (std::is_void_v<EdgeT>) {
        cls.def("remove_edges_if",
                [](grove_t& g, std::function<bool(key_t*)> predicate) {
                    return g.remove_edges_if([&predicate](const auto& e) -> bool {
                        return predicate(e.target);
                    });
                },
                py::arg("predicate"),
                R"pbdoc(
                    remove_edges_if(predicate) -> int

                    Remove every edge whose target satisfies
                    predicate(target: Key) -> bool, returning the number removed.
                    (This grove's edges carry no metadata, so the predicate gets
                    only the target Key; the universal Grove also passes the edge
                    metadata.)
                )pbdoc");
    } else {
        cls.def("remove_edges_if",
                [](grove_t& g,
                   std::function<bool(key_t*, const EdgeT&)> predicate) {
                    return g.remove_edges_if([&predicate](const auto& e) -> bool {
                        return predicate(e.target, e.metadata);
                    });
                },
                py::arg("predicate"),
                R"pbdoc(
                    remove_edges_if(predicate) -> int

                    Remove every edge for which
                    predicate(target: Key, metadata: object) -> bool returns True,
                    returning the number removed. The predicate receives the target
                    Key and the decoded edge metadata (the value passed to
                    add_edge).
                )pbdoc");
    }

    // ---- Key removal + storage compaction ----
    cls.def("remove_key",
            [](grove_t& g, const std::string& index, key_t* key) {
                return g.remove_key(index, key);
            },
            py::arg("index"), py::arg("key").none(true),
            R"pbdoc(
                Remove a key from the index's B+ tree, rebalancing as needed.

                Returns True if the key was found and removed, False otherwise
                (including a None key or an unknown index). All graph edges
                touching the key (incoming and outgoing) are also removed. The key
                remains in storage as a dead slot (not freed) — other Keys keep
                their pointers; only compact() reclaims the slot (and invalidates
                pointers — see its warning).
            )pbdoc")
       .def("compact", &grove_t::compact,
            R"pbdoc(
                Reclaim the dead storage slots left by remove_key() (storage
                shrinks to exactly the live key count).

                WARNING: this INVALIDATES every Key previously returned for this
                grove's indexed keys (by insert(), insert_bulk(), or yielded from
                intersect()/flanking()) — they become dangling and must NOT be
                used afterward (doing so is undefined behaviour). After compact(),
                re-discover keys via a fresh intersect()/flanking() query. Keys
                from add_external_key() are NOT affected.
            )pbdoc");

    // ---- External (graph-only) key (coordinate + data payload) ----
    {
        auto ext_fn = [](grove_t& g, const KeyT& key, DataT data) {
            return g.add_external_key(key, std::move(data));
        };
        const char* ext_doc = R"pbdoc(
                    Add a key (coordinate + data) that lives outside the B+ tree
                    index but can participate in the graph overlay.

                    Both the key and the data are copied into the Grove. Returns a
                    stable Key that remains valid as long as the Grove is alive.
                    External keys are not returned by intersect() queries. On the
                    universal Grove the data defaults to None.
                )pbdoc";
        if constexpr (grove_data_optional<DataT>) {
            cls.def("add_external_key", ext_fn,
                    py::arg("key"), py::arg("data") = DataT{},
                    py::return_value_policy::reference_internal, ext_doc);
        } else {
            cls.def("add_external_key", ext_fn,
                    py::arg("key"), py::arg("data"),
                    py::return_value_policy::reference_internal, ext_doc);
        }
    }

    // ---- Labelled edges (only on groves whose edge type is non-void; on the
    //      universal Grove the metadata is any JSON-serializable value) ----
    if constexpr (!std::is_void_v<EdgeT>) {
        cls.def("add_edge",
                [](grove_t& g, key_t* source, key_t* target, EdgeT data) {
                    g.add_edge(source, target, std::move(data));
                },
                py::arg("source").none(false), py::arg("target").none(false),
                py::arg("data"),
                R"pbdoc(
                    add_edge(source, target, data) -> None

                    Add a directed edge from source to target carrying a metadata
                    payload (any JSON-serializable value on the universal Grove).
                    This is an overload of add_edge(source, target); the two-argument
                    form attaches a None payload. Raises TypeError if either key is
                    None.
                )pbdoc");
        cls.def("get_edges",
                [](const grove_t& g, const key_t* source) {
                    return g.get_edges(source);
                },
                py::arg("source").none(false),
                R"pbdoc(
                    get_edges(source) -> list

                    The metadata payloads of all outgoing edges from source, in edge
                    order (parallel to get_neighbors(source)). Edges added without a
                    payload yield None.
                )pbdoc");
        cls.def("get_edge_list",
                [](py::object self, const key_t* source) {
                    const auto& g = self.cast<const grove_t&>();
                    py::list out;
                    for (const auto& e : g.get_edge_list(source)) {
                        // Pin each target Key to the Grove — issue #37.
                        out.append(py::make_tuple(
                            py::cast(e.target,
                                     py::return_value_policy::reference_internal,
                                     self),
                            py::cast(e.metadata)));
                    }
                    return out;
                },
                py::arg("source").none(false),
                R"pbdoc(
                    get_edge_list(source) -> list[tuple[Key, object]]

                    The outgoing edges from source as (target Key, metadata)
                    pairs — i.e. the zip of get_neighbors(source) and
                    get_edges(source). Edges added without a payload yield None
                    metadata. Each returned Key keeps this Grove alive.
                )pbdoc");
        cls.def("get_neighbors_if",
                [](py::object self, key_t* source,
                   std::function<bool(const EdgeT&)> predicate) {
                    // The predicate calls back into Python — keep the GIL held.
                    // Pin each Key to the Grove so extracted keys can't dangle
                    // after the list is dropped — issue #37.
                    return pinned_key_list(
                        self.cast<const grove_t&>().get_neighbors_if(
                            source, std::move(predicate)),
                        self);
                },
                py::arg("source").none(false), py::arg("predicate"),
                R"pbdoc(
                    get_neighbors_if(source, predicate) -> list[Key]

                    Target Keys of the outgoing edges from source whose edge
                    metadata satisfies predicate(metadata). The predicate receives
                    the decoded payload (e.g. the dict you stored). The returned Keys
                    point into this Grove's storage and are valid only while it is
                    alive.
                )pbdoc");
        cls.def("link_with",
                [](grove_t& g, const std::vector<key_t*>& keys,
                   std::function<std::optional<EdgeT>(key_t*, key_t*)> predicate) {
                    // The predicate calls back into Python — keep the GIL held.
                    g.link_if(keys, std::move(predicate));
                },
                py::arg("keys"), py::arg("predicate"),
                R"pbdoc(
                    link_with(keys, predicate) -> None

                    Like link_if(), but the predicate returns the edge metadata to
                    attach, or None to skip: predicate(keys[i], keys[i+1]) ->
                    Optional[object], applied to each adjacent pair. Use this to
                    label edges built over a bulk-inserted run of keys.
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
            // File write + zlib touches no Python objects (JSON payloads are
            // stored as strings); GIL released for the duration.
            py::call_guard<py::gil_scoped_release>(),
            R"pbdoc(
                Serialize the Grove (intervals + associated data + graph overlay)
                to a zlib-compressed binary file at the given path.
            )pbdoc")
       .def("to_sif",
            [](const grove_t& g, const std::string& path) {
                std::ofstream os(path);
                if (!os) {
                    throw std::runtime_error(
                        "Failed to open file for writing: " + path);
                }
                g.grove_to_sif(os);
                if (!os) {
                    throw std::runtime_error(
                        "Failed to write SIF to file: " + path);
                }
            },
            py::arg("path"),
            // Pure C++ tree walk + text write (key.to_string()); no Python.
            py::call_guard<py::gil_scoped_release>(),
            R"pbdoc(
                Write the grove to a SIF (Simple Interaction Format) text file for
                graph visualization (e.g. Cytoscape). Emits the B+ tree structure
                (`nodelink` / `leaflink` lines) and the graph-overlay edges
                (`keylink` lines) as tab-separated interactions; an empty grove
                writes an empty file.

                Note: line and index iteration order are not stable across runs.
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
            // File read + zlib + tree rebuild is pure C++ (payloads stay encoded
            // strings, decoded lazily on key.data); the returned Grove is wrapped
            // after the GIL is reacquired.
            py::call_guard<py::gil_scoped_release>(),
            R"pbdoc(
                Load a Grove previously written with serialize(). Returns a new
                Grove with the same intervals, associated data, and graph edges.
            )pbdoc");
}