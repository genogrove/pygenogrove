"""
Tests for the entry-deriving overloads of insert / insert_bulk on data-carrying
groves: passing a bare file entry (or a list of them) derives the grove's
canonical 0-based closed Interval key from the entry's native coordinates, so
callers never hand-convert.

  - BED is 0-based half-open [start, end)  -> Interval(start, end - 1)
  - GFF is 1-based inclusive [start, end]  -> Interval(start - 1, end - 1)

These overloads live on BedGrove/GffGrove alongside the explicit
insert(index, interval, data) / insert_bulk(index, [(interval, data)]) forms;
pybind resolves them by signature.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_bed_insert_entry_derives_closed_interval():
    """insert(index, BedEntry) — BED [0, 100) half-open -> Interval(0, 99) closed."""
    pg = _pg()
    g = pg.BedGrove(8)
    key = g.insert("chr1", pg.BedEntry("chr1", 0, 100))
    assert key.value == pg.Interval(0, 99)
    assert key.value.start == 0 and key.value.end == 99
    # payload keeps its native half-open coordinates
    assert key.data.start == 0 and key.data.end == 100


def test_gff_insert_entry_derives_closed_interval():
    """insert(index, GffEntry) — GFF [1, 100] 1-based -> Interval(0, 99) closed."""
    pg = _pg()
    g = pg.GffGrove(8)
    key = g.insert("chr1", pg.GffEntry("chr1", 1, 100, "gene"))
    assert key.value == pg.Interval(0, 99)
    # payload keeps its native 1-based coordinates
    assert key.data.start == 1 and key.data.end == 100


def test_bed_and_gff_same_region_yield_same_key():
    """The worked example: a BED and a GFF record for the SAME 100 bases derive
    the SAME Interval key, despite different native coordinate conventions."""
    pg = _pg()
    bed_key = pg.BedGrove(8).insert("chr1", pg.BedEntry("chr1", 0, 100))   # [0,100)
    gff_key = pg.GffGrove(8).insert("chr1", pg.GffEntry("chr1", 1, 100, "x"))  # [1,100] 1-based
    assert bed_key.value == gff_key.value == pg.Interval(0, 99)


def test_explicit_and_entry_insert_coexist():
    """The 3-arg explicit insert and the 2-arg entry insert both work (overloads)."""
    pg = _pg()
    g = pg.BedGrove(8)
    k_explicit = g.insert("chr1", pg.Interval(10, 20), pg.BedEntry("chr1", 10, 21))
    k_entry = g.insert("chr1", pg.BedEntry("chr1", 100, 151))
    assert k_explicit.value == pg.Interval(10, 20)
    assert k_entry.value == pg.Interval(100, 150)
    assert g.size() == 2


def test_entry_insert_is_queryable():
    """A key inserted from a BED entry is found by an overlapping query, and the
    half-open end is respected (the BED end position itself does not overlap)."""
    pg = _pg()
    g = pg.BedGrove(8)
    e = pg.BedEntry("chr1", 1000, 2000)   # half-open -> Interval(1000, 1999)
    e.name = "geneA"
    g.insert("chr1", e)

    assert len(list(g.intersect(pg.Interval(1500, 1500), "chr1"))) == 1
    assert len(list(g.intersect(pg.Interval(1999, 1999), "chr1"))) == 1
    assert len(list(g.intersect(pg.Interval(2000, 2000), "chr1"))) == 0
    assert list(g.intersect(pg.Interval(1500, 1500), "chr1"))[0].data.name == "geneA"


def test_insert_bulk_entries_derives_keys():
    """insert_bulk(index, [entries]) derives each key from its coordinates."""
    pg = _pg()
    g = pg.BedGrove(8)
    entries = []
    for i in range(30):
        e = pg.BedEntry("chr1", i * 100, i * 100 + 50)   # [s, s+50) -> [s, s+49]
        e.name = f"f{i}"
        entries.append(e)
    keys = g.insert_bulk("chr1", entries, presorted=True)

    assert len(keys) == 30
    assert g.size() == 30
    assert [k.value.start for k in keys] == [i * 100 for i in range(30)]
    assert [k.value.end for k in keys] == [i * 100 + 49 for i in range(30)]
    assert [k.data.name for k in keys] == [f"f{i}" for i in range(30)]


def test_insert_bulk_entries_sorts_unsorted():
    """insert_bulk(entries, presorted=False) sorts by derived interval, keeping pairing."""
    pg = _pg()
    g = pg.BedGrove(8)
    entries = []
    for i in range(20):
        e = pg.BedEntry("chr1", i * 100, i * 100 + 50)
        e.name = f"f{i}"
        entries.append(e)
    keys = g.insert_bulk("chr1", list(reversed(entries)))   # descending

    assert [k.value.start for k in keys] == [i * 100 for i in range(20)]
    assert [k.data.name for k in keys] == [f"f{i}" for i in range(20)]


def test_insert_bulk_entries_gff():
    """insert_bulk on GffGrove derives 0-based closed keys from 1-based coords."""
    pg = _pg()
    g = pg.GffGrove(8)
    entries = [pg.GffEntry("chr1", i * 100 + 1, i * 100 + 100, "gene")
               for i in range(10)]
    keys = g.insert_bulk("chr1", entries, presorted=True)
    assert len(keys) == 10
    assert [k.value.start for k in keys] == [i * 100 for i in range(10)]
    assert [k.value.end for k in keys] == [i * 100 + 99 for i in range(10)]


def test_explicit_pair_bulk_still_works():
    """The explicit (Interval, data) bulk form still resolves alongside the entry form."""
    pg = _pg()
    g = pg.BedGrove(8)
    items = [(pg.Interval(i * 100, i * 100 + 50), pg.BedEntry("chr1", i * 100, i * 100 + 51))
             for i in range(10)]
    keys = g.insert_bulk("chr1", items, presorted=True)
    assert len(keys) == 10
    assert [k.value.start for k in keys] == [i * 100 for i in range(10)]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])