/*
 * read_next_guarded — shared __next__ implementation for the streaming file
 * readers (BamReader/BedReader/GffReader/FastaReader/VcfReader).
 *
 * Each reader's __next__ releases the GIL around read_next() (pure I/O/parse,
 * touches no Python objects), but read_next() mutates reader-owned,
 * unsynchronized state (record/EOF counters, reused decode buffers). Two
 * Python threads driving the same reader object concurrently would race on
 * that state — real memory corruption, reachable purely from Python (issue
 * #84). Guard against it with a flag checked-then-set while the GIL is still
 * held, before releasing it for the read: the GIL itself serializes that
 * check-and-set, so no separate lock is needed. A concurrent call sees the
 * flag already set and raises RuntimeError instead of racing.
 */
#pragma once

#include <pybind11/pybind11.h>

#include <stdexcept>
#include <string>

namespace py = pybind11;

template <typename Reader, typename Entry>
Entry read_next_guarded(py::object self, const char* class_name) {
    static constexpr const char* flag_attr = "_pygg_reading";
    if (py::hasattr(self, flag_attr) && self.attr(flag_attr).cast<bool>()) {
        throw std::runtime_error(
            std::string(class_name) +
            ".__next__ called concurrently from another thread — drive one "
            "reader per thread.");
    }
    self.attr(flag_attr) = true;
    auto& r = self.cast<Reader&>();
    Entry entry;
    bool has_next;
    try {
        py::gil_scoped_release rel;
        has_next = r.read_next(entry);
    } catch (...) {
        self.attr(flag_attr) = false;
        throw;
    }
    self.attr(flag_attr) = false;
    if (!has_next) {
        throw py::stop_iteration();
    }
    return entry;
}