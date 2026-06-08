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

#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/bed_reader.hpp>
#include <genogrove/io/gff_reader.hpp>

namespace gdt = genogrove::data_type;
namespace gio = genogrove::io;

// BED: 0-based half-open [start, end) -> 0-based closed [start, end - 1].
inline gdt::interval interval_from_entry(const gio::bed_entry& e) {
    return gdt::interval(e.start, e.end == 0 ? 0 : e.end - 1);
}

// GFF/GTF: 1-based inclusive [start, end] -> 0-based closed [start - 1, end - 1].
inline gdt::interval interval_from_entry(const gio::gff_entry& e) {
    return gdt::interval(e.start == 0 ? 0 : e.start - 1,
                         e.end == 0 ? 0 : e.end - 1);
}

// True for data types that have an interval_from_entry overload (so the grove
// can derive a key from them). Gates the entry-deriving insert overloads.
template <typename T>
concept has_entry_interval = requires(const T& e) {
    { interval_from_entry(e) } -> std::same_as<gdt::interval>;
};