/*
 * pinned_key_list — build a Python list of Key objects from a range of key
 * pointers, tying EACH returned Key's lifetime to `parent` (reference_internal).
 *
 * Why this exists (issue #37): genogrove returns vectors of raw key pointers
 * (e.g. grove::get_neighbors, the bulk insert_data path, query_result::get_keys)
 * that point into a Grove's pointer-stable std::deque. Casting such a vector
 * with py::return_value_policy::reference_internal pins only the resulting list
 * to its parent — the individual Key elements get no keep-alive. A Key pulled
 * out of that list and outliving the list (and every other handle to the owning
 * Grove) then dangles: a use-after-free reachable from pure Python.
 *
 * Building the list element-by-element and casting each Key with the parent as
 * its keep-alive nurse fixes that: an extracted Key keeps `parent` alive, and
 * `parent` (the Grove, or a QueryResult that itself keeps its Grove alive) keeps
 * the underlying storage alive. The list stays indexable / len()-able, so the
 * public API is unchanged.
 */
#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

// `keys` is any range yielding key pointers (e.g. std::vector<key<...>*>).
// `parent` is the owning object whose lifetime each Key must extend — the Grove
// for grove methods, the QueryResult for QueryResult.keys.
template <typename Range>
py::list pinned_key_list(const Range& keys, py::handle parent) {
    py::list out;
    for (auto* k : keys) {
        out.append(py::cast(k, py::return_value_policy::reference_internal, parent));
    }
    return out;
}