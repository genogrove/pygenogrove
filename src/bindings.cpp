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
#include "data_type/json_value.hpp"
#include "data_type/registry.hpp"
#include "io/bam_reader.hpp"
#include "io/bed_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/gff_reader.hpp"
#include "structure/grove.hpp"

namespace py = pybind11;
namespace gio = genogrove::io;

// The universal Grove's JSON payload defaults to None (an absent payload), so
// dataless inserts can be written `g.insert(index, coord)`.
template <>
inline constexpr bool grove_data_optional<pygg::json_value> = true;

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

    // GenomicCoordinate — the standard, strand-aware key value type. Overlap is
    // strand-aware ('*' = wildcard, '.' = a concrete unstranded value).
    bind_genomic_coordinate(m);

    // The universal Grove: grove<genomic_coordinate, json_value>. Stores an
    // arbitrary JSON-serializable Python object (dict / list / scalar / None) as
    // the payload — insert(index, coord, data); key.data round-trips it back. It
    // serializes to a .gg whose payload is JSON text, so a C++
    // grove<genomic_coordinate, std::string> can still read the file.
    bind_grove<gdt::genomic_coordinate, pygg::json_value>(
        m, "Grove", "Key", "QueryResult", "FlankingResult");

    // Typed data-carrying groves over genomic_coordinate, kept for full C++
    // interop (typed binary .gg) and the BED/GFF helper accessors. BedEntry /
    // GffEntry must be registered before the groves reference them; the readers
    // yield them and the entry-deriving insert(index, entry) derives a stranded
    // coordinate from the BED6/GFF strand column.
    bind_bed_entry(m);
    bind_grove<gdt::genomic_coordinate, gio::bed_entry>(
        m, "BedGrove", "BedKey", "BedQueryResult", "BedFlankingResult");
    bind_bed_reader(m);

    bind_gff_entry(m);
    bind_grove<gdt::genomic_coordinate, gio::gff_entry>(
        m, "GffGrove", "GffKey", "GffQueryResult", "GffFlankingResult");
    bind_gff_reader(m);

    // SAM/BAM alignment reader: SamFlags / AlignmentFlags / SamEntry value types
    // and the BamReader iterator. sam_entry isn't serializable, so there's no
    // typed BamGrove — load alignments into the universal Grove via
    // SamEntry.to_coordinate() + .to_dict() (see the SamEntry docstring).
    bind_sam_entry(m);
    bind_bam_reader(m);

    // FASTA/FASTQ sequence reader: FastaEntry value type + FastaReader iterator.
    // Standalone (named sequences, not intervals — no grove integration).
    bind_fasta_entry(m);
    bind_fasta_reader(m);

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