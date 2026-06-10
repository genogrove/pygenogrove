"""
Tests for the dataless genomic-coordinate grove (grove<genomic_coordinate>),
exposed as GenomicCoordinateGrove / GenomicCoordinateKey / ...QueryResult /
...FlankingResult.

Mirrors genogrove/tests/structure/genomic_coordinate_grove_test.cpp over the
bound surface: strand-aware intersect (matching / wildcard / unstranded),
flanking, and a serialization round-trip that preserves strand.
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _gc(pg, strand, start, end):
    return pg.GenomicCoordinate(strand, start, end)


def test_creation_and_insert():
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    assert g.get_order() == 8
    assert g.size() == 0

    key = g.insert("chr1", _gc(pg, "+", 100, 200))
    assert g.size() == 1
    assert key.value == _gc(pg, "+", 100, 200)
    assert key.value.strand == "+"


def test_strand_specific_query_matches_only_same_strand():
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200))
    g.insert("chr1", _gc(pg, "-", 100, 200))
    g.insert("chr1", _gc(pg, ".", 100, 200))

    results = g.intersect(_gc(pg, "+", 150, 160), "chr1")
    strands = [k.value.strand for k in results]
    assert strands == ["+"]


def test_wildcard_query_finds_all_strands():
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200))
    g.insert("chr1", _gc(pg, "-", 100, 200))
    g.insert("chr1", _gc(pg, ".", 100, 200))

    results = g.intersect(_gc(pg, "*", 150, 160), "chr1")
    assert len(results) == 3
    assert sorted(k.value.strand for k in results) == ["+", "-", "."]


def test_unstranded_query_finds_only_unstranded():
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200))
    g.insert("chr1", _gc(pg, ".", 100, 200))

    results = g.intersect(_gc(pg, ".", 150, 160), "chr1")
    strands = [k.value.strand for k in results]
    assert strands == ["."]


def test_no_overlap_different_strands():
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200))

    results = g.intersect(_gc(pg, "-", 150, 160), "chr1")
    assert len(results) == 0


def test_no_spatial_overlap():
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200))

    results = g.intersect(_gc(pg, "+", 300, 400), "chr1")
    assert len(results) == 0


def test_flanking_bracketed_by_two_keys():
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200))
    g.insert("chr1", _gc(pg, "+", 500, 600))

    r = g.flanking(_gc(pg, "+", 300, 400), "chr1")
    assert r.predecessor.value == _gc(pg, "+", 100, 200)
    assert r.successor.value == _gc(pg, "+", 500, 600)


def test_flanking_empty_grove_returns_none():
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    r = g.flanking(_gc(pg, "+", 100, 200), "chr1")
    assert r.predecessor is None
    assert r.successor is None


def test_serialization_roundtrip_preserves_strand(tmp_path):
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200))
    g.insert("chr1", _gc(pg, "-", 300, 400))
    g.insert("chr2", _gc(pg, ".", 150, 250))

    path = str(tmp_path / "coords.gg")
    g.serialize(path)

    loaded = pg.GenomicCoordinateGrove.deserialize(path)
    assert loaded.size() == 3
    assert loaded.get_order() == 8

    # The '+' coordinate is still found only by a compatible-strand query.
    plus = loaded.intersect(_gc(pg, "+", 150, 160), "chr1")
    assert [k.value.strand for k in plus] == ["+"]
    minus = loaded.intersect(_gc(pg, "-", 350, 360), "chr1")
    assert [(k.value.strand, k.value.start, k.value.end) for k in minus] == [("-", 300, 400)]


def test_keys_keep_grove_alive():
    """A key materialized from a query keeps the grove alive (keep_alive)."""
    pg = _pg()
    g = pg.GenomicCoordinateGrove(8)
    g.insert("chr1", _gc(pg, "+", 100, 200))
    results = g.intersect(_gc(pg, "+", 150, 160), "chr1")
    keys = list(results)
    del g, results
    gc.collect()
    assert keys[0].value == _gc(pg, "+", 100, 200)