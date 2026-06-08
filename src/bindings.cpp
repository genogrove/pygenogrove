/*
 * pygenogrove — Python bindings for the genogrove C++ library.
 *
 * This translation unit is just the module entry point: it wires together the
 * per-module binding helpers, which are organized to mirror the genogrove
 * source tree (data_type/, io/, structure/).
 */
#include <pybind11/pybind11.h>

#include "data_type/interval.hpp"
#include "io/bed_reader.hpp"
#include "io/gff_reader.hpp"
#include "structure/interval_grove.hpp"

namespace py = pybind11;
namespace gio = genogrove::io;

PYBIND11_MODULE(pygenogrove, m) {
    m.doc() = R"pbdoc(
        pygenogrove: Python bindings for the genogrove C++ library

        A specialized B+ tree data structure optimized for genomic interval storage and querying.
    )pbdoc";

    // Key type shared by every grove instantiation.
    bind_interval(m);

    // Dataless interval grove: grove<interval> exposed as Grove / Key /
    // QueryResult / FlankingResult.
    bind_interval_grove<void>(m, "Grove", "Key", "QueryResult", "FlankingResult");

    // BED value types, then the data-carrying grove<interval, bed_entry>
    // exposed as BedGrove / BedKey / BedQueryResult / BedFlankingResult. BedEntry
    // must be registered before the grove references it. BedReader yields BedEntry.
    bind_bed_entry(m);
    bind_interval_grove<gio::bed_entry>(m, "BedGrove", "BedKey", "BedQueryResult",
                                        "BedFlankingResult");
    bind_bed_reader(m);

    // GFF/GTF value types, then the data-carrying grove<interval, gff_entry>
    // exposed as GffGrove / GffKey / GffQueryResult / GffFlankingResult.
    bind_gff_entry(m);
    bind_interval_grove<gio::gff_entry>(m, "GffGrove", "GffKey", "GffQueryResult",
                                        "GffFlankingResult");
    bind_gff_reader(m);

    m.attr("__version__") = "0.1.0";
}