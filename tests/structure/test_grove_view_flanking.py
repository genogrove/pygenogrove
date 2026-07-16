"""
Tests for GroveView.flanking() — nearest non-overlapping neighbours over a
serialized .gg, paging in only the descent-path blocks (genogrove #483).

The flanking logic itself is covered against the eager Grove in test_flanking.py;
here we only exercise the new view surface: parity with Grove.flanking(), the
predicate overload, keep-alive of the returned Keys, and that a query loads only
part of the file.
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _coord(pg, start, end):
    return pg.GenomicCoordinate(".", start, end)


def test_view_flanking_matches_eager(tmp_path):
    pg = _pg()
    g = pg.Grove(4)
    g.insert("chr1", _coord(pg, 100, 200))
    g.insert("chr1", _coord(pg, 500, 600))
    path = str(tmp_path / "f.gg")
    g.serialize(path)

    eager = pg.Grove.deserialize(path)
    view = pg.GroveView.open(path)

    q = _coord(pg, 300, 400)
    er = eager.flanking(q, "chr1")
    vr = view.flanking(q, "chr1")
    assert vr.predecessor.value == er.predecessor.value == _coord(pg, 100, 200)
    assert vr.successor.value == er.successor.value == _coord(pg, 500, 600)


def test_view_flanking_edges_and_missing_index(tmp_path):
    pg = _pg()
    g = pg.Grove(4)
    g.insert("chr1", _coord(pg, 100, 200))
    path = str(tmp_path / "e.gg")
    g.serialize(path)
    view = pg.GroveView.open(path)

    # query after the only key -> predecessor, no successor
    r = view.flanking(_coord(pg, 500, 600), "chr1")
    assert r.predecessor.value == _coord(pg, 100, 200)
    assert r.successor is None

    # unknown index -> both None
    r = view.flanking(_coord(pg, 500, 600), "chr2")
    assert r.predecessor is None and r.successor is None


def test_view_flanking_predicate(tmp_path):
    pg = _pg()
    g = pg.Grove(4)
    g.insert("chr1", _coord(pg, 100, 200))
    g.insert("chr1", _coord(pg, 300, 400))
    g.insert("chr1", _coord(pg, 500, 600))
    path = str(tmp_path / "p.gg")
    g.serialize(path)
    view = pg.GroveView.open(path)

    q = _coord(pg, 420, 430)
    # keep only keys starting below 250 -> 100-200 is the only predecessor match
    r = view.flanking(q, "chr1", lambda cand, query: cand.start < 250)
    assert r.predecessor.value == _coord(pg, 100, 200)
    assert r.successor is None

    # predicate exceptions propagate
    def boom(cand, query):
        raise ValueError("boom")

    with pytest.raises(ValueError):
        view.flanking(q, "chr1", boom)


def test_view_flanking_keys_keep_view_alive(tmp_path):
    """Chain key -> FlankingResult -> view (keep_alive + reference_internal). A
    regression is a use-after-free crash, not a failed assert."""
    pg = _pg()
    g = pg.Grove(4)
    g.insert("chr1", _coord(pg, 100, 200))
    g.insert("chr1", _coord(pg, 500, 600))
    path = str(tmp_path / "k.gg")
    g.serialize(path)

    view = pg.GroveView.open(path)
    r = view.flanking(_coord(pg, 300, 400), "chr1")
    pred, succ = r.predecessor, r.successor

    del view, r
    gc.collect()

    assert pred.value == _coord(pg, 100, 200)
    assert succ.value == _coord(pg, 500, 600)


def test_view_flanking_loads_only_part_of_the_file(tmp_path):
    pg = _pg()
    g = pg.Grove(3)
    for i in range(200):
        s = i * 100
        g.insert("chr1", _coord(pg, s, s + 50))
    path = str(tmp_path / "big.gg")
    g.serialize(path)

    view = pg.GroveView.open(path)
    assert view.block_count() > 20
    r = view.flanking(_coord(pg, 9970, 9980), "chr1")  # gap near one leaf
    assert r.predecessor is not None
    assert 0 < view.blocks_loaded() < view.block_count()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])