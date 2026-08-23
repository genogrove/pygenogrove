/*
 * read_next_guarded / guarded_getter — thread-safety guard shared by the
 * streaming file readers. read_next() and the getters below it all touch
 * reader-owned, unsynchronized state; a concurrent call from another thread
 * would race on it (issue #84). A flag checked-then-set while the GIL is
 * still held closes that — the GIL serializes the check-and-set, no lock
 * needed. Assumes GIL-enabled CPython (not free-threaded/PEP 703 builds,
 * which this project doesn't target).
 */
#pragma once

#include <pybind11/pybind11.h>

#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace pygg_reader_guard_detail {
inline constexpr const char* flag_attr = "_pygg_reading";

inline void check_not_reading(py::object& self, const char* class_name,
                               const char* what) {
    if (py::hasattr(self, flag_attr) && self.attr(flag_attr).cast<bool>()) {
        throw std::runtime_error(
            std::string(class_name) + "." + what +
            " called concurrently with __next__ from another thread — drive "
            "one reader per thread.");
    }
}
}  // namespace pygg_reader_guard_detail

template <typename Reader, typename Entry>
Entry read_next_guarded(py::object self, const char* class_name) {
    using namespace pygg_reader_guard_detail;
    check_not_reading(self, class_name, "__next__");
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

// Guards a plain getter (get_current_line/get_error_message) against running
// concurrently with an in-flight __next__ on another thread.
template <typename Reader, typename Ret>
Ret guarded_getter(py::object self, const char* class_name,
                    const char* method_name, Ret (Reader::*method)() const) {
    pygg_reader_guard_detail::check_not_reading(self, class_name, method_name);
    return (self.cast<Reader&>().*method)();
}