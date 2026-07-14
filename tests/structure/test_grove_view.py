"""
Tests for GroveView — the partial (random-access) reader over a serialized .gg.

Mirrors genogrove/tests/structure/grove_view_test.cpp over the surface the Python
bindings expose: a view opened on a file written by Grove.serialize() returns the
same results as an eager Grove.deserialize(), but pages in only the blocks a query
touches. One view family is bound (GroveView / BedGroveView / …) from a single
template, so the behaviour is exercised on the universal GroveView plus a
BedGroveView smoke test rather than 1:1 across every instantiation.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _coord(pg, start, end):
    return pg.GenomicCoordinate(".", start, end)


def test_view_matches_eager_deserialize(tmp_path):
    pg = _pg()

    g = pg.Grove(4)
    g.insert("chr1", _coord(pg, 100, 200), {"n": "a"})
    g.insert("chr1", _coord(pg, 300, 400), {"n": "b"})
    g.insert("chr2", _coord(pg, 150, 250), {"n": "c"})

    path = str(tmp_path / "view.gg")
    g.serialize(path)

    eager = pg.Grove.deserialize(path)
    view = pg.GroveView.open(path)

    def payloads(qr):
        return sorted(k.data["n"] for k in qr)

    q = _coord(pg, 150, 350)
    assert payloads(view.intersect(q, "chr1")) == payloads(eager.intersect(q, "chr1"))
    assert payloads(view.intersect(q, "chr1")) == ["a", "b"]
    # all-index overload
    assert payloads(view.intersect(_coord(pg, 0, 10_000))) == ["a", "b", "c"]


def test_view_loads_only_part_of_the_file(tmp_path):
    pg = _pg()

    # Order 3 + many intervals -> a multi-block tree, so a point query can prove
    # it loaded only a fraction of the blocks.
    g = pg.Grove(3)
    for i in range(200):
        s = i * 100
        g.insert("chr1", _coord(pg, s, s + 50))

    path = str(tmp_path / "big.gg")
    g.serialize(path)

    view = pg.GroveView.open(path)
    assert view.block_count() > 20, "need a multi-block tree for this to be meaningful"

    hits = view.intersect(_coord(pg, 1000, 1005), "chr1")  # one interval
    assert len(hits) == 1
    assert 0 < view.blocks_loaded() < view.block_count()

    # repr/str reflect the class name and block counts
    assert "GroveView" in repr(view)
    assert str(view.block_count()) in repr(view)
    assert "GroveView" in str(view)


def test_view_traverses_graph_edges(tmp_path):
    pg = _pg()

    g = pg.Grove(3)
    a = g.insert("chr1", _coord(pg, 100, 200))
    b = g.insert("chr1", _coord(pg, 300, 400))
    ext = g.add_external_key(_coord(pg, 5000, 5500))
    g.add_edge(a, b)
    g.add_edge(b, ext)

    path = str(tmp_path / "edges.gg")
    g.serialize(path)

    view = pg.GroveView.open(path)
    src = list(view.intersect(_coord(pg, 100, 200), "chr1"))[0]
    n1 = view.get_neighbors(src)
    assert [k.value.start for k in n1] == [300]

    # follow the chain onto the external key (a separate on-disk block)
    n2 = view.get_neighbors(n1[0])
    assert [k.value.start for k in n2] == [5000]


def test_view_reads_edge_metadata(tmp_path):
    """get_edges / get_neighbors_if surface edge payloads through a view
    without a full deserialize (genogrove #480)."""
    pg = _pg()

    g = pg.Grove(3)
    a = g.insert("chr1", _coord(pg, 100, 200))
    b = g.insert("chr1", _coord(pg, 300, 400))
    c = g.insert("chr1", _coord(pg, 500, 600))
    g.add_edge(a, b, {"w": 1})
    g.add_edge(a, c, {"w": 10})

    path = str(tmp_path / "meta.gg")
    g.serialize(path)

    view = pg.GroveView.open(path)
    src = list(view.intersect(_coord(pg, 100, 200), "chr1"))[0]

    # payloads parallel to get_neighbors, same as the mutable Grove
    paired = {n.value.start: e for n, e in zip(view.get_neighbors(src), view.get_edges(src))}
    assert paired == {300: {"w": 1}, 500: {"w": 10}}

    strong = view.get_neighbors_if(src, lambda meta: meta["w"] >= 10)
    assert [k.value.start for k in strong] == [500]

    # empty for a source with no edges; None for None source
    leaf = list(view.intersect(_coord(pg, 500, 600), "chr1"))[0]
    assert view.get_edges(leaf) == []
    with pytest.raises((TypeError, ValueError)):
        view.get_neighbors_if(None, lambda m: True)
    with pytest.raises((TypeError, ValueError)):
        view.get_edges(None)


def test_view_unknown_index_and_empty_grove(tmp_path):
    pg = _pg()

    populated = str(tmp_path / "pop.gg")
    g = pg.Grove(3)
    g.insert("chr1", _coord(pg, 10, 20))
    g.serialize(populated)
    view = pg.GroveView.open(populated)
    assert len(view.intersect(_coord(pg, 10, 20), "nope")) == 0

    empty = str(tmp_path / "empty.gg")
    pg.Grove(3).serialize(empty)
    ev = pg.GroveView.open(empty)
    assert ev.block_count() == 0
    assert len(ev.intersect(_coord(pg, 1, 2), "chr1")) == 0
    assert len(ev.intersect(_coord(pg, 1, 2))) == 0


def test_view_get_neighbors_rejects_none(tmp_path):
    pg = _pg()
    path = str(tmp_path / "one.gg")
    g = pg.Grove(3)
    g.insert("chr1", _coord(pg, 10, 20))
    g.serialize(path)
    view = pg.GroveView.open(path)
    with pytest.raises((TypeError, ValueError)):
        view.get_neighbors(None)


def test_open_missing_and_corrupt_raise(tmp_path):
    pg = _pg()
    with pytest.raises((RuntimeError, IOError, OSError)):
        pg.GroveView.open(str(tmp_path / "does_not_exist.gg"))

    bad = tmp_path / "garbage.gg"
    bad.write_bytes(b"not a format 0.2 grove stream")
    with pytest.raises((RuntimeError, ValueError, IOError, OSError)):
        pg.GroveView.open(str(bad))


def test_bed_grove_view_smoke(tmp_path):
    """The typed instantiation binds and round-trips like the universal one."""
    pg = _pg()

    g = pg.BedGrove(4)
    a = pg.BedEntry("chr7", 100, 200)
    a.name = "featA"
    b = pg.BedEntry("chr7", 500, 600)
    b.name = "featB"
    g.insert("chr7", _coord(pg, 100, 200), a)
    g.insert("chr7", _coord(pg, 500, 600), b)

    path = str(tmp_path / "bed.gg")
    g.serialize(path)

    view = pg.BedGroveView.open(path)
    hits = list(view.intersect(_coord(pg, 150, 550), "chr7"))
    assert sorted(k.data.name for k in hits) == ["featA", "featB"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])