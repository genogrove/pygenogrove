#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

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

    // Bind interval class
    py::class_<gdt::interval>(m, "Interval", R"pbdoc(
        A genomic interval with start and end coordinates.

        Parameters
        ----------
        start : int
            Start position of the interval (0-based)
        end : int
            End position of the interval (exclusive)
    )pbdoc")
        .def(py::init<>())
        .def(py::init<size_t, size_t>(),
             py::arg("start"),
             py::arg("end"))
        .def_property("start",
                     &gdt::interval::get_start,
                     &gdt::interval::set_start,
                     "Start position of the interval")
        .def_property("end",
                     &gdt::interval::get_end,
                     &gdt::interval::set_end,
                     "End position of the interval")
        .def("__str__", [](const gdt::interval& i) {
            return std::to_string(i.get_start()) + "-" + std::to_string(i.get_end());
        })
        .def("__repr__", [](const gdt::interval& i) {
            return "Interval(" + std::to_string(i.get_start()) + ", " +
                   std::to_string(i.get_end()) + ")";
        })
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def(py::self == py::self)
        .def_static("overlap", &gdt::interval::overlap,
                   py::arg("a"),
                   py::arg("b"),
                   "Check if two intervals overlap");

    // Bind key class for intervals
    py::class_<gdt::key<gdt::interval>>(m, "Key", R"pbdoc(
        A key object wrapping an interval in the grove structure.

        This is typically returned by insert operations and can be used
        to build graph relationships.
    )pbdoc")
        .def_property_readonly("value",
                              [](const gdt::key<gdt::interval>& k) {
                                  return k.get_value();
                              },
                              "The interval value of this key")
        .def("__str__", [](const gdt::key<gdt::interval>& k) {
            return k.to_string();
        });

    // Bind query_result class
    py::class_<gdt::query_result<gdt::interval>>(m, "QueryResult", R"pbdoc(
        Result of a query operation containing matching intervals.
    )pbdoc")
        .def_property_readonly("query",
                              &gdt::query_result<gdt::interval>::get_query,
                              "The query interval used for this search")
        .def_property_readonly("keys",
                              [](const gdt::query_result<gdt::interval>& qr) {
                                  return qr.get_keys();
                              },
                              "List of matching keys")
        .def("__len__", [](const gdt::query_result<gdt::interval>& qr) {
            return qr.get_keys().size();
        })
        .def("__iter__", [](const gdt::query_result<gdt::interval>& qr) {
            return py::make_iterator(qr.get_keys().begin(), qr.get_keys().end());
        }, py::keep_alive<0, 1>());

    // Bind grove class for intervals
    py::class_<ggs::grove<gdt::interval>>(m, "Grove", R"pbdoc(
        A B+ tree container for efficient genomic interval storage and querying.

        The grove supports multi-index operations, where each index (e.g., chromosome)
        maintains its own B+ tree structure.

        Parameters
        ----------
        order : int, optional
            Maximum branching factor of the B+ tree (default: 3)
            Controls the maximum number of keys per node (order - 1)
    )pbdoc")
        .def(py::init<>())
        .def(py::init<int>(), py::arg("order"))
        .def("__str__", [](const ggs::grove<gdt::interval>& g) {
            return "Grove(size=" + std::to_string(g.size()) + ")";
        })
        .def("__repr__", [](const ggs::grove<gdt::interval>& g) {
            return "Grove(order=" + std::to_string(g.get_order()) +
                   ", size=" + std::to_string(g.size()) + ")";
        })
        .def("size", &ggs::grove<gdt::interval>::size,
             "Get total number of intervals across all indices")
        .def("get_order", &ggs::grove<gdt::interval>::get_order,
             "Get the order (branching factor) of the B+ tree")
        .def("insert", [](ggs::grove<gdt::interval>& g,
                         const std::string& index,
                         const gdt::interval& interval) {
            gdt::key<gdt::interval> key(interval);
            return g.insert(index, key);
        }, py::arg("index"), py::arg("interval"),
           py::return_value_policy::reference,
           R"pbdoc(
               Insert an interval into the grove at the specified index.

               Parameters
               ----------
               index : str
                   The index name (e.g., chromosome name like "chr1")
               interval : Interval
                   The interval to insert

               Returns
               -------
               Key
                   Pointer to the inserted key in the tree
           )pbdoc")
        .def("insert_sorted", [](ggs::grove<gdt::interval>& g,
                                const std::string& index,
                                const gdt::interval& interval) {
            gdt::key<gdt::interval> key(interval);
            return g.insert_sorted(index, key);
        }, py::arg("index"), py::arg("interval"),
           py::return_value_policy::reference,
           R"pbdoc(
               Insert a pre-sorted interval for optimal performance.

               Assumes the interval is greater than all existing intervals in the index.
               Significantly faster than regular insert() for sorted data.

               Parameters
               ----------
               index : str
                   The index name (e.g., chromosome name like "chr1")
               interval : Interval
                   The interval to insert (must be > all existing intervals)

               Returns
               -------
               Key
                   Pointer to the inserted key in the tree
           )pbdoc")
        .def("intersect",
             py::overload_cast<gdt::interval&>(&ggs::grove<gdt::interval>::intersect),
             py::arg("query"),
             R"pbdoc(
                 Find all intervals that overlap with the query across all indices.

                 Parameters
                 ----------
                 query : Interval
                     The query interval to search for

                 Returns
                 -------
                 QueryResult
                     Result object containing all overlapping intervals
             )pbdoc")
        .def("intersect",
             py::overload_cast<gdt::interval&, const std::string&>(
                 &ggs::grove<gdt::interval>::intersect),
             py::arg("query"), py::arg("index"),
             R"pbdoc(
                 Find all intervals that overlap with the query in a specific index.

                 Parameters
                 ----------
                 query : Interval
                     The query interval to search for
                 index : str
                     The index name (e.g., chromosome name) to search within

                 Returns
                 -------
                 QueryResult
                     Result object containing all overlapping intervals from the index
             )pbdoc");

    m.attr("__version__") = "0.1.0";
}