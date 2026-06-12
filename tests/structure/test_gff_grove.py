"""
Tests for the data-carrying grove GffGrove (grove<interval, gff_entry>).

Mirrors the dataless grove tests and the BedGrove tests, extended for the
GffEntry payload (GTF/GFF3 records with a column-9 attributes dict).
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _gene(pg, seqid="chr1", start=100, end=200, **attrs):
    e = pg.GffEntry(seqid, start, end, "gene")
    if attrs:
        e.attributes = attrs
    return e


def test_creation():
    """GffGrove(order) reports its order and starts empty."""
    pg = _pg()
    g = pg.GffGrove(100)
    assert g.get_order() == 100
    assert g.size() == 0


def test_default_order():
    """A default GffGrove uses order 3."""
    pg = _pg()
    assert pg.GffGrove().get_order() == 3


def test_str_repr_use_class_name():
    """__str__/__repr__ identify the grove by its Python class name."""
    pg = _pg()
    g = pg.GffGrove(7)
    g.insert("chr1", pg.GenomicCoordinate(".", 100, 200), _gene(pg))
    assert str(g).startswith("GffGrove(")
    rep = repr(g)
    assert rep.startswith("GffGrove(")
    assert "7" in rep and "1" in rep  # order and size


def test_insert_carries_value_and_data():
    """insert(index, interval, data) returns a key exposing both .value and .data."""
    pg = _pg()
    g = pg.GffGrove(100)
    key = g.insert("chr1", pg.GenomicCoordinate(".", 100, 200),
                   _gene(pg, gene_id="ENSG1"))

    assert g.size() == 1
    assert key.value.start == 100
    assert key.data.type == "gene"
    assert key.data.get_gene_id() == "ENSG1"


def test_intersect_preserves_data():
    """Keys returned by intersect() carry the inserted GffEntry payload."""
    pg = _pg()
    g = pg.GffGrove(100)
    g.insert("chr1", pg.GenomicCoordinate(".", 100, 200), _gene(pg, "chr1", 100, 200, ID="g1"))
    g.insert("chr1", pg.GenomicCoordinate(".", 300, 400),
             pg.GffEntry("chr1", 300, 400, "exon"))

    hits = list(g.intersect(pg.GenomicCoordinate(".", 150, 350), "chr1"))
    assert len(hits) == 2
    types = sorted(k.data.type for k in hits)
    assert types == ["exon", "gene"]


def test_multi_index_intersect_counts():
    """intersect with/without an index mirrors the dataless grove semantics."""
    pg = _pg()
    g = pg.GffGrove(100)
    g.insert("chr1", pg.GenomicCoordinate(".", 100, 200), _gene(pg, "chr1", 100, 200))
    g.insert("chr2", pg.GenomicCoordinate(".", 150, 250), _gene(pg, "chr2", 150, 250))
    g.insert("chr3", pg.GenomicCoordinate(".", 300, 400), _gene(pg, "chr3", 300, 400))

    assert len(g.intersect(pg.GenomicCoordinate(".", 175, 225))) == 2          # chr1 + chr2
    assert len(g.intersect(pg.GenomicCoordinate(".", 175, 225), "chr1")) == 1


def test_graph_overlay_with_external_data_key():
    """Graph overlay works on GffGrove, incl. add_external_key(interval, data)."""
    pg = _pg()
    g = pg.GffGrove(5)
    exon = g.insert("chr1", pg.GenomicCoordinate(".", 1000, 1200),
                    pg.GffEntry("chr1", 1000, 1200, "exon"))
    enhancer = g.add_external_key(pg.GenomicCoordinate(".", 5000, 5500),
                                  pg.GffEntry("chr1", 5000, 5500, "enhancer"))
    g.add_edge(exon, enhancer)

    assert g.size() == 1            # external key is not indexed
    assert g.has_edge(exon, enhancer)
    assert g.edge_count() == 1

    neighbors = g.get_neighbors(exon)
    assert len(neighbors) == 1
    assert neighbors[0].data.type == "enhancer"

    assert len(g.intersect(pg.GenomicCoordinate(".", 5000, 5500))) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
