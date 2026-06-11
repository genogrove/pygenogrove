"""
Tests for the data-carrying genomic-coordinate groves:
grove<genomic_coordinate, bed_entry>  -> GenomicCoordinateBedGrove
grove<genomic_coordinate, gff_entry>  -> GenomicCoordinateGffGrove

These pair strand-aware coordinates with a BED/GFF payload. They reuse the
generic bind_grove template, so the surface matches the other data-carrying
groves EXCEPT the entry-deriving insert(index, entry) overloads, which are gated
to interval keys (you build the GenomicCoordinate explicitly).
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _gc(pg, strand, start, end):
    return pg.GenomicCoordinate(strand, start, end)


def _bed(pg, name):
    e = pg.BedEntry("chr1", 0, 1)
    e.name = name
    return e


# ----------------------------------------------------- GenomicCoordinateBedGrove

def test_insert_carries_coordinate_and_data():
    pg = _pg()
    g = pg.GenomicCoordinateBedGrove(8)
    k = g.insert("chr1", _gc(pg, "+", 100, 200), _bed(pg, "feat1"))
    assert g.size() == 1
    assert k.value == _gc(pg, "+", 100, 200)
    assert k.value.strand == "+"
    assert k.data.name == "feat1"

    # .data is a live, mutable reference into grove storage.
    k.data.name = "renamed"
    assert k.data.name == "renamed"


def test_strand_aware_intersect_returns_matching_data():
    pg = _pg()
    g = pg.GenomicCoordinateBedGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200), _bed(pg, "plus"))
    g.insert("chr1", _gc(pg, "-", 100, 200), _bed(pg, "minus"))

    res = g.intersect(_gc(pg, "+", 150, 160), "chr1")
    assert [k.data.name for k in res] == ["plus"]


def test_wildcard_query_returns_all_strands_with_data():
    pg = _pg()
    g = pg.GenomicCoordinateBedGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200), _bed(pg, "plus"))
    g.insert("chr1", _gc(pg, "-", 100, 200), _bed(pg, "minus"))

    res = g.intersect(_gc(pg, "*", 150, 160), "chr1")
    assert sorted(k.data.name for k in res) == ["minus", "plus"]


def test_bulk_insert_pairs():
    pg = _pg()
    g = pg.GenomicCoordinateBedGrove(8)
    items = [(_gc(pg, "+", i * 10, i * 10 + 5), _bed(pg, f"f{i}")) for i in range(5)]
    keys = g.insert_bulk("chr1", items)
    assert len(keys) == 5
    assert g.size() == 5
    res = g.intersect(_gc(pg, "+", 0, 1000), "chr1")
    assert sorted(k.data.name for k in res) == [f"f{i}" for i in range(5)]


def test_strand_filtered_flanking_with_data():
    pg = _pg()
    g = pg.GenomicCoordinateBedGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200), _bed(pg, "up_plus"))
    g.insert("chr1", _gc(pg, "-", 300, 400), _bed(pg, "near_minus"))
    q = _gc(pg, "+", 500, 510)

    r = g.flanking(q, "chr1", lambda c, query: c.strand == query.strand)
    assert r.predecessor.data.name == "up_plus"
    assert r.successor is None


def test_serialization_roundtrip_preserves_strand_and_data(tmp_path):
    pg = _pg()
    g = pg.GenomicCoordinateBedGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200), _bed(pg, "plus"))
    g.insert("chr1", _gc(pg, "-", 300, 400), _bed(pg, "minus"))

    path = str(tmp_path / "gc_bed.gg")
    g.serialize(path)
    loaded = pg.GenomicCoordinateBedGrove.deserialize(path)

    assert loaded.size() == 2
    res = loaded.intersect(_gc(pg, "+", 150, 160), "chr1")
    assert [(k.value.strand, k.data.name) for k in res] == [("+", "plus")]


def test_no_entry_deriving_insert_overload():
    """The insert(index, entry) form is gated to interval keys, so it's absent."""
    pg = _pg()
    g = pg.GenomicCoordinateBedGrove(8)
    with pytest.raises(TypeError):
        g.insert("chr1", _bed(pg, "x"))  # missing the explicit coordinate


# ----------------------------------------------------- GenomicCoordinateGffGrove

def test_gff_variant_insert_and_strand_query():
    pg = _pg()
    g = pg.GenomicCoordinateGffGrove(8)
    plus = pg.GffEntry("chr1", 1000, 2000, "gene")
    minus = pg.GffEntry("chr1", 1000, 2000, "gene")
    kp = g.insert("chr1", _gc(pg, "+", 1000, 2000), plus)
    g.insert("chr1", _gc(pg, "-", 1000, 2000), minus)

    assert kp.value.strand == "+"
    assert kp.data.type == "gene"

    res = g.intersect(_gc(pg, "+", 1500, 1600), "chr1")
    assert len(res) == 1
    assert [k.value.strand for k in res] == ["+"]