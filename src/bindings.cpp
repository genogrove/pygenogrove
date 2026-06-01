#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>

#include <genogrove/data_type/interval.hpp>
#include <genogrove/data_type/key.hpp>
#include <genogrove/data_type/query_result.hpp>
#include <genogrove/structure/grove/grove.hpp>

namespace py = pybind11;
namespace ggs = genogrove::structure;
namespace gdt = genogrove::data_type;

PYBIND11_MODULE(pygenogrove, m) {
    m.doc() = R"pbdoc(
        pygenogrove: Python bindings for the genogrove C++ library

        A specialized B+ tree data structure optimized for genomic interval storage and querying.
    )pbdoc";

    py::class_<gdt::interval>(m, "Interval", R"pbdoc(
        A genomic interval with closed [start, end] coordinates (0-based, both inclusive).

        Parameters
        ----------
        start : int
            Start position of the interval (0-based, inclusive)
        end : int
            End position of the interval (0-based, inclusive)
    )pbdoc")
        .def(py::init<>())
        .def(py::init<size_t, size_t>(),
             py::arg("start"),
             py::arg("end"))
        .def_property_readonly("start",
                     &gdt::interval::get_start,
                     "Start position of the interval (read-only — see set_range)")
        .def_property_readonly("end",
                     &gdt::interval::get_end,
                     "End position of the interval (read-only — see set_range)")
        .def("set_range", &gdt::interval::set_range,
             py::arg("start"), py::arg("end"),
             R"pbdoc(
                 Atomically set both endpoints.

                 Do NOT call this on an interval that has already been inserted into
                 a Grove — mutating a stored key silently corrupts B+ tree ordering.
                 Use this only on intervals not yet inserted (e.g. queries you intend
                 to reuse).
             )pbdoc")
        .def("__str__", &gdt::interval::to_string)
        .def("__repr__", [](const gdt::interval& i) {
            return "Interval(" + std::to_string(i.get_start()) + ", " +
                   std::to_string(i.get_end()) + ")";
        })
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def(py::self == py::self)
        .def_static("overlaps", &gdt::interval::overlaps,
                   py::arg("a"),
                   py::arg("b"),
                   "Check if two intervals overlap (closed-interval semantics)");

    py::class_<gdt::key<gdt::interval>>(m, "Key", R"pbdoc(
        A key object wrapping an interval in the grove structure.

        Returned by Grove.insert() and yielded by QueryResult iteration. Wraps a
        pointer into the grove's storage, so the Key remains valid only as long
        as the originating Grove is alive.
    )pbdoc")
        .def_property_readonly("value",
                              [](const gdt::key<gdt::interval>& k) {
                                  return k.get_value();
                              },
                              "The interval value of this key (returned by value)")
        .def("__str__", [](const gdt::key<gdt::interval>& k) {
            return k.to_string();
        });

    py::class_<gdt::query_result<gdt::interval>>(m, "QueryResult", R"pbdoc(
        Result of a query operation containing matching intervals.
    )pbdoc")
        .def_property_readonly("query",
                              &gdt::query_result<gdt::interval>::get_query,
                              "The query interval used for this search")
        .def_property_readonly("keys",
                              &gdt::query_result<gdt::interval>::get_keys,
                              py::return_value_policy::reference_internal,
                              "List of matching keys (pointers into the grove)")
        .def("__len__", [](const gdt::query_result<gdt::interval>& qr) {
            return qr.get_keys().size();
        })
        .def("__iter__", [](const gdt::query_result<gdt::interval>& qr) {
            return py::make_iterator(qr.get_keys().begin(), qr.get_keys().end());
        }, py::keep_alive<0, 1>());

    py::class_<ggs::grove<gdt::interval>>(m, "Grove", R"pbdoc(
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
        .def("__str__", [](const ggs::grove<gdt::interval>& g) {
            return "Grove(size=" + std::to_string(g.indexed_vertex_count()) + ")";
        })
        .def("__repr__", [](const ggs::grove<gdt::interval>& g) {
            return "Grove(order=" + std::to_string(g.get_order()) +
                   ", size=" + std::to_string(g.indexed_vertex_count()) + ")";
        })
        .def("__len__", &ggs::grove<gdt::interval>::indexed_vertex_count)
        .def("size", &ggs::grove<gdt::interval>::indexed_vertex_count,
             "Number of indexed intervals across all indices (alias of len)")
        .def("indexed_vertex_count",
             &ggs::grove<gdt::interval>::indexed_vertex_count,
             "Number of indexed intervals (B+ tree leaf keys)")
        .def("get_order", &ggs::grove<gdt::interval>::get_order,
             "Get the order (branching factor) of the B+ tree")
        .def("insert", [](ggs::grove<gdt::interval>& g,
                         const std::string& index,
                         const gdt::interval& interval) {
            gdt::key<gdt::interval> key(interval);
            return g.insert(index, key);
        }, py::arg("index"), py::arg("interval"),
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
                   Stable reference to the inserted key. Remains valid as long
                   as the Grove is alive.
           )pbdoc")
        .def("intersect",
             py::overload_cast<const gdt::interval&>(
                 &ggs::grove<gdt::interval>::intersect),
             py::arg("query"),
             R"pbdoc(
                 Find all intervals that overlap with the query across all indices.
             )pbdoc")
        .def("intersect",
             py::overload_cast<const gdt::interval&, std::string_view>(
                 &ggs::grove<gdt::interval>::intersect),
             py::arg("query"), py::arg("index"),
             R"pbdoc(
                 Find all intervals that overlap with the query in a specific index.
             )pbdoc");

    m.attr("__version__") = "0.1.0";
}
