"""
Tests for Grove.flanking() — the nearest non-overlapping predecessor/successor
of a query (a FlankingResult with `.predecessor` / `.successor`, each Key or None).

Ports the interval cases from genogrove grove_flanking_test.cpp. (The strand /
numeric cases there need the genomic_coordinate / numeric key types, not yet
bound.)
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_empty_grove_returns_none():
    pg = _pg()
    g = pg.Grove(8)
    r = g.flanking(pg.Interval(100, 200), "chr1")
    assert r.predecessor is None
    assert r.successor is None


def test_missing_index_returns_none():
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(50, 60))
    r = g.flanking(pg.Interval(100, 200), "chr2")
    assert r.predecessor is None
    assert r.successor is None


def test_query_before_all_keys_only_successor():
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(500, 600))
    r = g.flanking(pg.Interval(100, 200), "chr1")
    assert r.predecessor is None
    assert r.successor.value == pg.Interval(500, 600)


def test_query_after_all_keys_only_predecessor():
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(100, 200))
    r = g.flanking(pg.Interval(500, 600), "chr1")
    assert r.predecessor.value == pg.Interval(100, 200)
    assert r.successor is None


def test_query_bracketed_by_two_keys():
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(100, 200))
    g.insert("chr1", pg.Interval(500, 600))
    r = g.flanking(pg.Interval(300, 400), "chr1")
    assert r.predecessor.value == pg.Interval(100, 200)
    assert r.successor.value == pg.Interval(500, 600)


def test_overlapping_keys_are_skipped():
    """A key overlapping the query is never a flanking neighbour."""
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(50, 60))      # far predecessor
    g.insert("chr1", pg.Interval(100, 350))    # overlaps query -> skipped
    g.insert("chr1", pg.Interval(500, 600))    # far successor
    r = g.flanking(pg.Interval(300, 400), "chr1")
    assert r.predecessor.value == pg.Interval(50, 60)
    assert r.successor.value == pg.Interval(500, 600)


def test_abutting_keys_have_zero_gap():
    """Closed-coordinate abutting keys (gap 0) are valid flanking neighbours."""
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(50, 99))      # ends at 99, abuts query start 100
    g.insert("chr1", pg.Interval(201, 300))    # starts at 201, abuts query end 200
    q = pg.Interval(100, 200)
    r = g.flanking(q, "chr1")
    assert r.predecessor.value == pg.Interval(50, 99)
    assert r.successor.value == pg.Interval(201, 300)
    # gap distances (closed coords) are exactly zero
    assert q.start - r.predecessor.value.end - 1 == 0
    assert r.successor.value.start - q.end - 1 == 0


def test_chooses_closest_predecessor_among_many():
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(0, 10))
    g.insert("chr1", pg.Interval(20, 30))
    g.insert("chr1", pg.Interval(40, 50))      # closest to query
    r = g.flanking(pg.Interval(100, 200), "chr1")
    assert r.predecessor.value == pg.Interval(40, 50)
    assert r.successor is None


def test_nested_intervals_pick_closest_end():
    """With nested upstream intervals, the predecessor is the one with the
    largest end (smallest gap), not the sort-order maximum."""
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(50, 100))     # outer — larger end, closer
    g.insert("chr1", pg.Interval(80, 90))      # inner — larger sort key, farther
    r = g.flanking(pg.Interval(200, 300), "chr1")
    assert r.predecessor.value == pg.Interval(50, 100)


def test_multi_leaf_tree_picks_global_nearest():
    """Order 4 with 16 keys forces splits; pruning must still find the global
    nearest across leaves."""
    pg = _pg()
    g = pg.Grove(4)
    for i in range(16):
        g.insert("chr1", pg.Interval(i * 100, i * 100 + 20))
    r = g.flanking(pg.Interval(440, 460), "chr1")   # gap between [400,420] and [500,520]
    assert r.predecessor is not None
    assert r.successor is not None
    assert r.predecessor.value.end == 420
    assert r.successor.value.start == 500


def test_flanking_keys_keep_grove_alive():
    """Keys returned by flanking() keep the grove alive via the chain
    key -> FlankingResult -> grove (keep_alive on flanking() + reference_internal
    on .predecessor/.successor). A regression would be a use-after-free (crash),
    not a failed assertion. Mirrors test_query_result.py::test_keys_keep_grove_alive."""
    pg = _pg()
    g = pg.Grove(8)
    g.insert("chr1", pg.Interval(100, 200))
    g.insert("chr1", pg.Interval(500, 600))

    r = g.flanking(pg.Interval(300, 400), "chr1")
    pred = r.predecessor
    succ = r.successor

    del g, r            # drop the grove and the result; only the keys remain
    gc.collect()

    assert pred.value == pg.Interval(100, 200)
    assert succ.value == pg.Interval(500, 600)


def test_flanking_carries_data_on_bed_grove():
    """flanking() works on data-carrying groves; the returned keys expose .data."""
    pg = _pg()
    g = pg.BedGrove(8)
    up = pg.BedEntry("chr1", 50, 60)
    up.name = "upstream"
    down = pg.BedEntry("chr1", 500, 600)
    down.name = "downstream"
    g.insert("chr1", up)        # entry-deriving insert
    g.insert("chr1", down)

    r = g.flanking(pg.Interval(300, 400), "chr1")
    assert r.predecessor.data.name == "upstream"
    assert r.successor.data.name == "downstream"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])