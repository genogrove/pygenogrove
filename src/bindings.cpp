/*
 * pygenogrove — Python bindings for the genogrove C++ library.
 *
 * This translation unit is just the module entry point: it wires together the
 * per-module binding helpers, which are organized to mirror the genogrove
 * source tree (data_type/, io/, structure/).
 */
#include <pybind11/pybind11.h>

#include <genogrove/config/version.hpp>

#include "data_type/interval.hpp"
#include "io/bed_reader.hpp"
#include "io/gff_reader.hpp"
#include "structure/grove.hpp"

namespace py = pybind11;
namespace gio = genogrove::io;

// Stringify-and-join genogrove's integer version macros into "MAJOR.MINOR.PATCH".
#define PYGENOGROVE_STR2(x) #x
#define PYGENOGROVE_STR(x) PYGENOGROVE_STR2(x)
#define PYGENOGROVE_GENOGROVE_VERSION                                          \
    PYGENOGROVE_STR(genogrove_VERSION_MAJOR)                                   \
    "." PYGENOGROVE_STR(genogrove_VERSION_MINOR) "." PYGENOGROVE_STR(          \
        genogrove_VERSION_PATCH)

PYBIND11_MODULE(pygenogrove, m) {
    m.doc() = R"pbdoc(
        pygenogrove: Python bindings for the genogrove C++ library

        A specialized B+ tree data structure optimized for genomic interval storage and querying.
    )pbdoc";

    // Interval key value type, shared by every interval-keyed grove.
    bind_interval(m);

    // Dataless interval grove: grove<interval> exposed as Grove / Key /
    // QueryResult / FlankingResult.
    bind_grove<gdt::interval, void>(m, "Grove", "Key", "QueryResult",
                                    "FlankingResult");

    // BED value types, then the data-carrying grove<interval, bed_entry>
    // exposed as BedGrove / BedKey / BedQueryResult / BedFlankingResult. BedEntry
    // must be registered before the grove references it. BedReader yields BedEntry.
    bind_bed_entry(m);
    bind_grove<gdt::interval, gio::bed_entry>(m, "BedGrove", "BedKey",
                                              "BedQueryResult", "BedFlankingResult");
    bind_bed_reader(m);

    // GFF/GTF value types, then the data-carrying grove<interval, gff_entry>
    // exposed as GffGrove / GffKey / GffQueryResult / GffFlankingResult.
    bind_gff_entry(m);
    bind_grove<gdt::interval, gio::gff_entry>(m, "GffGrove", "GffKey",
                                              "GffQueryResult", "GffFlankingResult");
    bind_gff_reader(m);

    // __version__ is single-sourced from pyproject.toml via CMake; __genogrove_version__
    // reports the genogrove the wheel was built against (independent SemVer — the two
    // version lines move on their own cadence).
    m.attr("__version__") = PYGENOGROVE_VERSION;
    m.attr("__genogrove_version__") = PYGENOGROVE_GENOGROVE_VERSION;
}