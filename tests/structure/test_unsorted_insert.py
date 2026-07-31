"""
Out-of-order ``insert()`` keeps every key spatially queryable.

Regression coverage for pygenogrove #68 / genogrove #517: insert routing used
the node's bounding box (a subtree *minimum* for interval-like keys) instead of
the subtree maximum, so a key that did not sort last landed in a leaf too far
right. ``size()`` and the graph overlay still saw it; ``intersect()`` did not.
Fixed in genogrove v0.25.6 (genogrove #518).

Both the in-memory path and the ``deserialize()`` path are exercised — the
latter is the workflow #68 was filed on, and the only trigger for the subtree
maxima rebuild v0.25.6 added to deserialization.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _gc(pg, start, end):
    return pg.GenomicCoordinate(".", start, end)


# 500 non-overlapping coordinates: [i*100, i*100+50], so an exact query on one
# of them can only ever hit that one key.
COORDS = [(i * 100, i * 100 + 50) for i in range(500)]

# Stride rotation: seven ascending passes, each picking every 7th coordinate.
# Unique (every coordinate inserted exactly once), deterministic, and every
# pass after the first inserts *below* the tree's current maximum — which is
# the misroute case. Do not "simplify" this back into sorted order: appending
# in sorted order is the one insert pattern the bug never affected.
INSERT_ORDER = [c for offset in range(7) for c in COORDS[offset::7]]


def _assert_all_findable(pg, g):
    for start, end in COORDS:
        hits = g.intersect(_gc(pg, start, end), "chr1")
        found = [k for k in hits if (k.value.start, k.value.end) == (start, end)]
        assert len(found) == 1, f"lost key at {start}"


def test_unsorted_insert_stays_queryable():
    pg = _pg()
    g = pg.Grove(3)
    for start, end in INSERT_ORDER:
        g.insert("chr1", _gc(pg, start, end))

    assert g.size() == len(COORDS)
    _assert_all_findable(pg, g)


def test_insert_into_deserialized_grove_stays_queryable(tmp_path):
    """The #68 reproducer: augment a grove loaded from disk, then query it."""
    pg = _pg()
    g = pg.Grove(3)
    for start, end in INSERT_ORDER:
        g.insert("chr1", _gc(pg, start, end))

    path = str(tmp_path / "unsorted.gg")
    g.serialize(path)
    loaded = pg.Grove.deserialize(path)

    # mid-tree inserts into the loaded grove: 25 coordinates interleaved with
    # the existing ones (the +75 offset keeps them clear of every stored key)
    added = [(i * 100 + 75, i * 100 + 80) for i in range(0, 500, 20)]
    for start, end in added:
        loaded.insert("chr1", _gc(pg, start, end))

    assert loaded.size() == len(COORDS) + len(added)
    for start, end in added:
        hits = loaded.intersect(_gc(pg, start, end), "chr1")
        assert any((k.value.start, k.value.end) == (start, end) for k in hits), (
            f"lost key at {start}"
        )
    # the keys that came back from disk are still there too
    _assert_all_findable(pg, loaded)