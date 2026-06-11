"""
Tests for the data-carrying grove BedGrove (grove<interval, bed_entry>).

Mirrors the dataless grove tests (tests/structure/test_dataless_grove.py),
extended to cover the associated BedEntry data that BedGrove keys carry. The
BedKey value/data/lifetime behavior lives in tests/data_type/test_bed_key.py.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_creation():
    """BedGrove(order) reports its order and starts empty."""
    pg = _pg()
    g = pg.BedGrove(100)
    assert g.get_order() == 100
    assert g.size() == 0


def test_default_order():
    """A default BedGrove uses order 3."""
    pg = _pg()
    assert pg.BedGrove().get_order() == 3


def test_str_repr_use_class_name():
    """__str__/__repr__ identify the grove by its Python class name, not 'Grove'."""
    pg = _pg()
    g = pg.BedGrove(7)
    g.insert("chr1", pg.GenomicCoordinate(".", 100, 200), pg.BedEntry("chr1", 100, 200))
    assert str(g).startswith("BedGrove(")
    rep = repr(g)
    assert rep.startswith("BedGrove(")
    assert "7" in rep and "1" in rep  # order and size


def test_insert_carries_value_and_data():
    """insert(index, interval, data) returns a key exposing both .value and .data."""
    pg = _pg()
    g = pg.BedGrove(100)
    key = g.insert("chr1", pg.GenomicCoordinate(".", 100, 200), pg.BedEntry("chr1", 100, 201))

    assert g.size() == 1
    assert key.value.start == 100
    assert key.value.end == 200
    assert key.data.chrom == "chr1"
    assert key.data.start == 100
    assert key.data.end == 201


def test_intersect_preserves_data():
    """Keys returned by intersect() carry the inserted BedEntry payload."""
    pg = _pg()
    g = pg.BedGrove(100)
    g.insert("chr1", pg.GenomicCoordinate(".", 100, 200), pg.BedEntry("chr1", 100, 201))
    g.insert("chr1", pg.GenomicCoordinate(".", 300, 400), pg.BedEntry("chr1", 300, 401))

    hits = list(g.intersect(pg.GenomicCoordinate(".", 150, 350), "chr1"))
    assert len(hits) == 2
    starts = sorted(k.data.start for k in hits)
    assert starts == [100, 300]


def test_multi_index_intersect_counts():
    """intersect with/without an index mirrors the dataless grove semantics."""
    pg = _pg()
    g = pg.BedGrove(100)
    g.insert("chr1", pg.GenomicCoordinate(".", 100, 200), pg.BedEntry("chr1", 100, 200))
    g.insert("chr2", pg.GenomicCoordinate(".", 150, 250), pg.BedEntry("chr2", 150, 250))
    g.insert("chr3", pg.GenomicCoordinate(".", 300, 400), pg.BedEntry("chr3", 300, 400))

    assert len(g.intersect(pg.GenomicCoordinate(".", 175, 225))) == 2          # chr1 + chr2
    assert len(g.intersect(pg.GenomicCoordinate(".", 175, 225), "chr1")) == 1


def test_graph_overlay_with_external_data_key():
    """Graph overlay works on BedGrove, incl. add_external_key(interval, data)."""
    pg = _pg()
    g = pg.BedGrove(5)
    exon = g.insert("chr1", pg.GenomicCoordinate(".", 1000, 1200), pg.BedEntry("chr1", 1000, 1200))
    enhancer = g.add_external_key(pg.GenomicCoordinate(".", 5000, 5500),
                                  pg.BedEntry("chr1", 5000, 5500))
    g.add_edge(exon, enhancer)

    # External key is not indexed.
    assert g.size() == 1
    assert g.has_edge(exon, enhancer)
    assert g.edge_count() == 1

    neighbors = g.get_neighbors(exon)
    assert len(neighbors) == 1
    assert neighbors[0].value.start == 5000
    assert neighbors[0].data.start == 5000

    # External keys are excluded from spatial queries.
    assert len(g.intersect(pg.GenomicCoordinate(".", 5000, 5500))) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])