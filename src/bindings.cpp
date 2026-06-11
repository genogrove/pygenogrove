/*
 * pygenogrove — Python bindings for the genogrove C++ library.
 *
 * This translation unit is just the module entry point: it wires together the
 * per-module binding helpers, which are organized to mirror the genogrove
 * source tree (data_type/, io/, structure/).
 */
#include <pybind11/pybind11.h>

#include <genogrove/config/version.hpp>

#include "data_type/genomic_coordinate.hpp"
#include "data_type/interval.hpp"
#include "data_type/registry.hpp"
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

    // GenomicCoordinate key value type (stranded interval), then the dataless
    // grove<genomic_coordinate> exposed as GenomicCoordinateGrove /
    // GenomicCoordinateKey / GenomicCoordinateQueryResult /
    // GenomicCoordinateFlankingResult. Overlap is strand-aware ('*' = wildcard).
    bind_genomic_coordinate(m);
    bind_grove<gdt::genomic_coordinate, void>(
        m, "GenomicCoordinateGrove", "GenomicCoordinateKey",
        "GenomicCoordinateQueryResult", "GenomicCoordinateFlankingResult");

    // BED value types, then the data-carrying grove<interval, bed_entry>
    // exposed as BedGrove / BedKey / BedQueryResult / BedFlankingResult. BedEntry
    // must be registered before the grove references it. BedReader yields BedEntry.
    bind_bed_entry(m);
    bind_grove<gdt::interval, gio::bed_entry>(m, "BedGrove", "BedKey",
                                              "BedQueryResult", "BedFlankingResult");
    // Strand-aware data-carrying grove: grove<genomic_coordinate, bed_entry>.
    // No entry-deriving insert (that's gated to interval keys); construct the
    // GenomicCoordinate explicitly and insert it with a BedEntry payload.
    bind_grove<gdt::genomic_coordinate, gio::bed_entry>(
        m, "GenomicCoordinateBedGrove", "GenomicCoordinateBedKey",
        "GenomicCoordinateBedQueryResult", "GenomicCoordinateBedFlankingResult");
    bind_bed_reader(m);

    // GFF/GTF value types, then the data-carrying grove<interval, gff_entry>
    // exposed as GffGrove / GffKey / GffQueryResult / GffFlankingResult.
    bind_gff_entry(m);
    bind_grove<gdt::interval, gio::gff_entry>(m, "GffGrove", "GffKey",
                                              "GffQueryResult", "GffFlankingResult");
    // Strand-aware data-carrying grove: grove<genomic_coordinate, gff_entry>.
    bind_grove<gdt::genomic_coordinate, gio::gff_entry>(
        m, "GenomicCoordinateGffGrove", "GenomicCoordinateGffKey",
        "GenomicCoordinateGffQueryResult", "GenomicCoordinateGffFlankingResult");
    bind_gff_reader(m);

    // String interning registry: registry<std::string> exposed as StringRegistry
    // (a process-wide singleton, key == payload). Tagged / key->payload registry
    // variants are additional bind_registry<...> instantiations, not yet exposed.
    bind_registry<std::string>(m, "StringRegistry");

    // __version__ is single-sourced from pyproject.toml via CMake; __genogrove_version__
    // reports the genogrove the wheel was built against (independent SemVer — the two
    // version lines move on their own cadence).
    m.attr("__version__") = PYGENOGROVE_VERSION;
    m.attr("__genogrove_version__") = PYGENOGROVE_GENOGROVE_VERSION;
}