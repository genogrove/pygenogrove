/*
 * interval_from_entry — derive the grove's canonical Interval key (0-based,
 * closed [start, end]) from a file-format entry's native coordinates.
 *
 * Each file format uses a different coordinate convention; this is the one
 * place that conversion lives, so the entry-deriving overloads of
 * Grove.insert()/insert_bulk() are foolproof and users never hand-convert:
 *   - BED  is 0-based half-open [start, end)   -> Interval(start, end - 1)
 *   - GFF  is 1-based inclusive [start, end]   -> Interval(start - 1, end - 1)
 *
 * Adding a conversion overload here automatically enables the entry-deriving
 * insert overloads for that data type (gated by the has_entry_interval concept).
 */
#pragma once

#include <concepts>
#include <optional>
#include <stdexcept>
#include <string>

#include <genogrove/data_type/genomic_coordinate.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/bed_reader.hpp>
#include <genogrove/io/gff_reader.hpp>

namespace gdt = genogrove::data_type;
namespace gio = genogrove::io;

// BED: 0-based half-open [start, end) -> 0-based closed [start, end - 1].
// A real BED feature has end > start; an empty/inverted range has no closed
// Interval representation, so reject it explicitly (rather than letting `end-1`
// underflow into a huge interval, or producing an inverted one).
inline gdt::interval interval_from_entry(const gio::bed_entry& e) {
    if (e.end <= e.start) {
        throw std::invalid_argument(
            "BedEntry has an empty/invalid range [" + std::to_string(e.start) +
            ", " + std::to_string(e.end) +
            "); cannot derive a closed Interval key (BED end must be > start)");
    }
    return gdt::interval(e.start, e.end - 1);
}

// GFF/GTF: 1-based inclusive [start, end] -> 0-based closed [start - 1, end - 1].
// Valid GFF coordinates are 1-based with start >= 1 and end >= start; reject
// anything else explicitly.
inline gdt::interval interval_from_entry(const gio::gff_entry& e) {
    if (e.start == 0 || e.end < e.start) {
        throw std::invalid_argument(
            "GffEntry has an invalid 1-based range [" + std::to_string(e.start) +
            ", " + std::to_string(e.end) +
            "]; cannot derive a closed Interval key (need start >= 1, end >= start)");
    }
    return gdt::interval(e.start - 1, e.end - 1);
}

// True for data types that have an interval_from_entry overload (so the grove
// can derive a key from them). Gates the entry-deriving insert overloads.
template <typename T>
concept has_entry_interval = requires(const T& e) {
    { interval_from_entry(e) } -> std::same_as<gdt::interval>;
};

// ---------------------------------------------------------------------------
// genomic_coordinate (stranded) key derivation — the grove's standard key.
// Reuses interval_from_entry for the (validated) coordinate conversion and
// layers the strand on top: a BED6/GFF strand column maps straight through;
// a missing strand or GFF '?' becomes '.' (unstranded). genomic_coordinate's
// own ctor only accepts '+','-','.','*', so anything else normalizes to '.'.
inline char strand_char_from_entry(std::optional<char> strand) {
    if (strand && (*strand == '+' || *strand == '-' || *strand == '.' ||
                   *strand == '*')) {
        return *strand;
    }
    return '.';
}

inline gdt::genomic_coordinate genomic_coordinate_from_entry(const gio::bed_entry& e) {
    gdt::interval iv = interval_from_entry(e);
    return gdt::genomic_coordinate(strand_char_from_entry(e.strand),
                                   iv.get_start(), iv.get_end());
}

inline gdt::genomic_coordinate genomic_coordinate_from_entry(const gio::gff_entry& e) {
    gdt::interval iv = interval_from_entry(e);
    return gdt::genomic_coordinate(strand_char_from_entry(e.strand),
                                   iv.get_start(), iv.get_end());
}

// True for data types that can derive a genomic_coordinate key. Gates the
// entry-deriving insert overloads on the (now standard) genomic_coordinate grove.
template <typename T>
concept has_entry_coordinate = requires(const T& e) {
    { genomic_coordinate_from_entry(e) } -> std::same_as<gdt::genomic_coordinate>;
};