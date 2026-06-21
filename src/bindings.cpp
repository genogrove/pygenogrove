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
#include "data_type/kmer.hpp"
#include "data_type/numeric.hpp"
#include "data_type/registry.hpp"
#include "io/bam_reader.hpp"
#include "io/bed_reader.hpp"
#include "io/fasta_index.hpp"
#include "io/fasta_reader.hpp"
#include "io/filetype_detector.hpp"
#include "io/gff_reader.hpp"
#include "io/vcf_reader.hpp"
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

    // The universal Grove: grove<genomic_coordinate, json_value, json_value>.
    // Stores an arbitrary JSON-serializable Python object (dict / list / scalar /
    // None) as the payload — insert(index, coord, data); key.data round-trips it
    // back. Graph edges also carry a JSON-serializable payload (add_edge(s, t,
    // data) / get_edges / get_neighbors_if / link_with). It serializes to a .gg
    // whose payload and edge metadata are JSON text, so a C++
    // grove<genomic_coordinate, std::string, std::string> can still read the file
    // (an edgeless .gg is also readable by grove<genomic_coordinate, std::string>).
    bind_grove<gdt::genomic_coordinate, pygg::json_value, pygg::json_value>(
        m, "Grove", "Key", "QueryResult", "FlankingResult");

    // Alternative point key types (overlap = exact equality, no range semantics).
    // Each gets the same universal surface as Grove — optional JSON payload
    // (insert without data -> None), labelled JSON edges, serialization. Numeric
    // wraps an int (ids / timestamps); Kmer is a 2-bit-encoded DNA k-mer (k <= 32,
    // a membership dictionary).
    bind_numeric(m);
    bind_grove<gdt::numeric, pygg::json_value, pygg::json_value>(
        m, "NumericGrove", "NumericKey", "NumericQueryResult",
        "NumericFlankingResult");

    bind_kmer(m);
    bind_grove<gdt::kmer, pygg::json_value, pygg::json_value>(
        m, "KmerGrove", "KmerKey", "KmerQueryResult", "KmerFlankingResult");

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

    // VCF/BCF variant reader: SampleGenotype / VcfEntry value types + VcfReader
    // iterator. Like SAM, vcf_entry isn't serializable (variant-valued INFO /
    // nested samples), so there's no typed VcfGrove — load into the universal
    // Grove via VcfEntry.to_coordinate() + .to_dict().
    bind_vcf_reader(m);

    // FASTA/FASTQ sequence reader: FastaEntry value type + FastaReader iterator.
    // Standalone (named sequences, not intervals — no grove integration).
    bind_fasta_entry(m);
    bind_fasta_reader(m);

    // Random-access FASTA: FastaIndex(path) over htslib faidx — fetch a region's
    // bases by (name, start, end). Pairs with GenomicCoordinate (fetch a
    // feature's sequence); requires a writable directory to build a missing .fai.
    bind_fasta_index(m);

    // File-type detector: Filetype / CompressionType enums + FiletypeDetector.
    bind_filetype_detector(m);

    // Universal interning registry: registry<std::string, void, json_value>
    // exposed as Registry — a process-wide singleton mapping a string identity
    // (gene_id / chrom / transcript_id) to any JSON-serializable payload. One
    // bound class covers every payload shape, mirroring how the universal Grove
    // uses json_value instead of a per-type template zoo. The two-arg
    // intern(key, payload) form comes from the generic template; the single-arg
    // intern(value) sugar below stores the string as its own payload so
    // get(id) returns it back (plain string interning).
    using registry_t = gdt::registry<std::string, void, pygg::json_value>;
    auto reg = bind_registry<std::string, void, pygg::json_value>(m, "Registry");
    reg.def(
        "intern",
        [](registry_t& r, const std::string& value) {
            return r.intern(value, py::cast<pygg::json_value>(py::str(value)));
        },
        py::arg("value"),
        R"pbdoc(
            Intern a string as both key and payload and return its stable id.
            Idempotent (deduplicated); get(id) returns the string back.
            Convenience for plain string interning — use intern(key, payload)
            to attach a distinct JSON payload.
        )pbdoc");

    // __version__ is single-sourced from pyproject.toml via CMake; __genogrove_version__
    // reports the genogrove the wheel was built against (independent SemVer — the two
    // version lines move on their own cadence).
    m.attr("__version__") = PYGENOGROVE_VERSION;
    m.attr("__genogrove_version__") = PYGENOGROVE_GENOGROVE_VERSION;
}