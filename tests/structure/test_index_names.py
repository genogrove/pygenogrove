"""
get_index_names() on the eager groves lists the index (chromosome) names, in
unspecified order — the same meaning as GroveView.get_index_names().
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _inserts(pg):
    """(grove class, insert callable taking (grove, index, i)) for every eager grove."""
    gc = lambda i: pg.GenomicCoordinate(".", 10 * i, 10 * i + 5)
    return {
        "Grove": (pg.Grove, lambda g, idx, i: g.insert(idx, gc(i))),
        "BedGrove": (pg.BedGrove, lambda g, idx, i: g.insert(idx, gc(i), pg.BedEntry(idx, 10 * i, 10 * i + 6))),
        "GffGrove": (pg.GffGrove, lambda g, idx, i: g.insert(idx, gc(i), pg.GffEntry(idx, 10 * i + 1, 10 * i + 6, "gene"))),
        "NumericGrove": (pg.NumericGrove, lambda g, idx, i: g.insert(idx, pg.Numeric(i))),
        "KmerGrove": (pg.KmerGrove, lambda g, idx, i: g.insert(idx, pg.Kmer("ACGT"))),
    }


@pytest.mark.parametrize("name", ["Grove", "BedGrove", "GffGrove", "NumericGrove", "KmerGrove"])
def test_get_index_names(name):
    """Every eager grove reports the indices it holds; an empty grove has none."""
    pg = _pg()
    cls, insert = _inserts(pg)[name]
    g = cls()
    assert g.get_index_names() == []

    for i, idx in enumerate(["chr1", "chr2", "chrX"]):
        insert(g, idx, i)
    insert(g, "chr1", 9)                      # a second key in an existing index adds no name
    assert sorted(g.get_index_names()) == ["chr1", "chr2", "chrX"]


def test_get_index_names_survives_round_trip(tmp_path):
    """A deserialized grove lists the same indices as the one that was written."""
    pg = _pg()
    g = pg.Grove()
    g.insert("chr1", pg.GenomicCoordinate("+", 1, 2))
    g.insert("chr2", pg.GenomicCoordinate("-", 3, 4))
    path = str(tmp_path / "g.gg")
    g.serialize(path)
    assert sorted(pg.Grove.deserialize(path).get_index_names()) == ["chr1", "chr2"]