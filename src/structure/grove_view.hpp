/*
 * Binding for ggs::grove_view<KeyT, DataT, EdgeT> — a read-only, *partial*
 * reader over a serialized (format 0.2) grove. Mirrors genogrove
 * structure/grove/grove_view.hpp. Where Grove.deserialize() eagerly loads every
 * block, a GroveView pages in only the blocks a query walks (cached, no
 * eviction) — query a large on-disk .gg without loading it whole.
 *
 * One template instantiated per concrete type tuple from bindings.cpp, producing
 * a distinct Python class each time (GroveView = grove_view<genomic_coordinate,
 * json_value, json_value>, BedGroveView = grove_view<genomic_coordinate,
 * bed_entry>, …). It reuses the Key / QueryResult classes already registered by
 * the matching bind_grove<KeyT, DataT, EdgeT>, so it registers nothing new.
 *
 * The surface is query-only: open / intersect / get_neighbors plus the
 * blocks_loaded / block_count partial-load counters. There is no insert or
 * serialize — a view never mutates the grove.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>
#include <string_view>
#include <type_traits>

#include <genogrove/data_type/key.hpp>
#include <genogrove/structure/grove/grove_view.hpp>

#include "../data_type/key_list.hpp"
#include "../data_type/query_result.hpp"

namespace py = pybind11;
namespace ggs = genogrove::structure;
namespace gdt = genogrove::data_type;

template <typename KeyT, typename DataT, typename EdgeT = void>
void bind_grove_view(py::module_& m, const char* view_name) {
    using view_t = ggs::grove_view<KeyT, DataT, EdgeT>;
    using key_t = gdt::key<KeyT, DataT>;

    py::class_<view_t>(m, view_name, R"pbdoc(
        A read-only, partial reader over a serialized (format 0.2) .gg grove.

        Unlike Grove.deserialize(), which loads the whole grove into memory, a
        GroveView reads only the block directory up front, then pages in
        individual blocks on demand as a query descends the tree — and caches
        them for the view's lifetime (no eviction). Use it to query a large
        on-disk index without loading it whole.

        Query-only: it has no insert() or serialize(). Not thread-safe (a query
        mutates the block cache). The Keys it returns point into the view's own
        storage and are valid only while the GroveView is alive.

        Create one with GroveView.open(path); a file written by Grove.serialize()
        is read directly (data_offset=0).
    )pbdoc")
        // Non-copyable AND non-movable (owns the file handle + a z_stream), so
        // there is no py::init: open() is the only entry point. The prvalue it
        // returns direct-initializes the heap object via guaranteed copy elision
        // (no move ctor needed), and pybind's default unique_ptr holder adopts
        // it.
        .def_static(
            "open",
            [](const std::string& path, std::streamoff data_offset) {
                return std::unique_ptr<view_t>(
                    new view_t(view_t::open(path, data_offset)));
            },
            py::arg("path"), py::arg("data_offset") = 0,
            R"pbdoc(
                open(path, data_offset=0) -> GroveView

                Open a serialized grove for partial reading. `path` is a file
                written by Grove.serialize() (a bare grove stream; data_offset=0).
                Pass a non-zero data_offset only for a .gg embedded after a
                leading header (e.g. a genogrove CLI index). Raises RuntimeError
                if the file cannot be opened, the magic is wrong, the source is
                not seekable, or the directory is malformed.
            )pbdoc")

        // keep_alive<0, 1>: the returned QueryResult (and the Keys it yields)
        // point into the view's block cache, so the view must outlive the
        // result. No call_guard — a query pages in blocks (mutating the cache)
        // and the view is not thread-safe, so the GIL stays held.
        .def(
            "intersect",
            [](view_t& v, const KeyT& query) { return v.intersect(query); },
            py::arg("query"), py::keep_alive<0, 1>(),
            R"pbdoc(
                Find all intervals overlapping the query across all indices,
                loading only the blocks the search touches.
            )pbdoc")
        .def(
            "intersect",
            [](view_t& v, const KeyT& query, std::string_view index) {
                return v.intersect(query, index);
            },
            py::arg("query"), py::arg("index"), py::keep_alive<0, 1>(),
            R"pbdoc(
                Find all intervals overlapping the query within a single index,
                loading only the blocks on the descent path and overlapping
                leaves.
            )pbdoc")

        .def(
            "get_neighbors",
            [](py::object self, key_t* source) {
                // Pin each Key to the view so an extracted neighbor can't dangle
                // after the list is dropped — issue #37. Loads each target's
                // block on demand (the cross-chromosome hop).
                return pinned_key_list(
                    self.cast<view_t&>().get_neighbors(source), self);
            },
            py::arg("source").none(false),
            R"pbdoc(
                Return the target Keys directly reachable from source via graph
                edges, paging in each target's block on demand. `source` must be a
                Key this GroveView produced (via intersect() or a prior
                get_neighbors()). The returned Keys are valid only while the view
                is alive. Raises ValueError if source is None.
            )pbdoc")

        .def("blocks_loaded", &view_t::blocks_loaded,
             "Number of blocks paged in so far — proof a query loaded only part "
             "of the file (compare with block_count()).")
        .def("block_count", &view_t::block_count,
             "Total number of blocks in the serialized grove (0 for an empty "
             "grove).");
}