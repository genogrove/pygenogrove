/*
 * guarded_reader / read_next_guarded / guarded_getter — thread-safety guard
 * shared by the streaming file readers. read_next() and the getters below it
 * all touch reader-owned, unsynchronized state; a concurrent call from another
 * thread would race on it. The Python class binds guarded_reader<R>, which adds
 * an atomic "reading" flag the Python side cannot see or tamper with; the
 * atomic exchange makes the check-and-set race-free without relying on the GIL.
 */
#pragma once

#include <pybind11/pybind11.h>

#include <atomic>
#include <stdexcept>
#include <string>

namespace py = pybind11;

template <typename Reader>
struct guarded_reader : Reader {
    using Reader::Reader;
    std::atomic<bool> reading{false};
};

namespace pygg_reader_guard_detail {
[[noreturn]] inline void throw_busy(const char* class_name, const char* what) {
    throw std::runtime_error(
        std::string(class_name) + "." + what +
        " called concurrently with __next__ from another thread — drive "
        "one reader per thread.");
}
}  // namespace pygg_reader_guard_detail

template <typename Reader, typename Entry>
Entry read_next_guarded(guarded_reader<Reader>& r, const char* class_name) {
    if (r.reading.exchange(true)) {
        pygg_reader_guard_detail::throw_busy(class_name, "__next__");
    }
    Entry entry;
    bool has_next;
    try {
        py::gil_scoped_release rel;
        has_next = r.read_next(entry);
    } catch (...) {
        r.reading = false;
        throw;
    }
    r.reading = false;
    if (!has_next) {
        throw py::stop_iteration();
    }
    return entry;
}

// Guards a plain getter (get_current_line/get_error_message/__repr__) against
// running concurrently with an in-flight __next__ on another thread.
template <typename Reader, typename Base, typename Ret>
Ret guarded_getter(const guarded_reader<Reader>& r, const char* class_name,
                   const char* method_name, Ret (Base::*method)() const) {
    if (r.reading.load()) {
        pygg_reader_guard_detail::throw_busy(class_name, method_name);
    }
    return (r.*method)();
}