"""
Tests for the sorted / bulk insertion fast paths on the data-carrying groves
(BedGrove / GffGrove): insert_sorted and insert_bulk.

Mirrors the genogrove key_type_grove_test harness (assert_sorted_insert,
assert_bulk_insert, assert_bulk_insert_returns_handles, assert_edge_creation).
These fast paths require associated data, so they exist on BedGrove/GffGrove but
not the dataless Grove.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _items(pg, n, start=0, step=100):
    """n ascending (Interval, BedEntry) records, named f0..f{n-1}."""
    out = []
    for i in range(n):
        s = start + i * step
        e = pg.BedEntry("chr1", s, s + 50)
        e.name = f"f{i}"
        out.append((pg.Interval(s, s + 50), e))
    return out


def test_insert_sorted_single(tmp_path=None):
    """insert_sorted appends ascending records; all are queryable. Mirrors assert_sorted_insert."""
    pg = _pg()
    g = pg.BedGrove(8)
    keys = [g.insert_sorted("chr1", iv, d) for iv, d in _items(pg, 30)]

    assert len(keys) == 30
    assert g.size() == 30
    assert keys[5].data.name == "f5"
    hits = list(g.intersect(pg.Interval(0, 30 * 100), "chr1"))
    assert len(hits) == 30


def test_insert_bulk_presorted():
    """insert_bulk(presorted=True) inserts a sorted batch in one call. Mirrors assert_bulk_insert."""
    pg = _pg()
    g = pg.BedGrove(8)
    keys = g.insert_bulk("chr1", _items(pg, 50), presorted=True)

    assert len(keys) == 50
    assert g.size() == 50
    assert len(list(g.intersect(pg.Interval(0, 50 * 100), "chr1"))) == 50


def test_insert_bulk_sorts_unsorted_input():
    """insert_bulk (presorted=False, default) sorts the batch by interval first."""
    pg = _pg()
    g = pg.BedGrove(8)
    items = _items(pg, 40)
    keys = g.insert_bulk("chr1", list(reversed(items)))   # descending input

    assert len(keys) == 40
    assert g.size() == 40
    # returned handles are in sorted (insertion) order regardless of input order
    starts = [k.value.start for k in keys]
    assert starts == sorted(starts)
    assert len(list(g.intersect(pg.Interval(0, 40 * 100), "chr1"))) == 40


def test_bulk_returns_handles_with_correct_data():
    """The returned handles map 1:1 to the inserted records. Mirrors assert_bulk_insert_returns_handles."""
    pg = _pg()
    g = pg.BedGrove(8)
    keys = g.insert_bulk("chr1", _items(pg, 12), presorted=True)
    assert [k.data.name for k in keys] == [f"f{i}" for i in range(12)]
    assert [k.value.start for k in keys] == [i * 100 for i in range(12)]


def test_bulk_keys_usable_for_edges():
    """Handles from insert_bulk are real keys — build an edge chain. Mirrors assert_edge_creation."""
    pg = _pg()
    g = pg.BedGrove(8)
    keys = g.insert_bulk("chr1", _items(pg, 10), presorted=True)
    for i in range(len(keys) - 1):
        g.add_edge(keys[i], keys[i + 1])
    assert g.edge_count() == 9
    assert all(g.has_edge(keys[i], keys[i + 1]) for i in range(9))


def test_bulk_equivalent_to_individual_insert():
    """Bulk insert yields the same query results as one-by-one insert()."""
    pg = _pg()
    items = _items(pg, 40)

    g_one = pg.BedGrove(8)
    for iv, d in items:
        g_one.insert("chr1", iv, d)

    g_bulk = pg.BedGrove(8)
    g_bulk.insert_bulk("chr1", items, presorted=True)

    assert g_one.size() == g_bulk.size() == 40
    q = pg.Interval(550, 2550)
    assert len(list(g_one.intersect(q, "chr1"))) == len(list(g_bulk.intersect(q, "chr1")))


def test_bulk_empty_batch():
    """An empty batch inserts nothing and returns no handles."""
    pg = _pg()
    g = pg.BedGrove(8)
    assert g.insert_bulk("chr1", []) == []
    assert g.size() == 0


def test_bulk_append_to_nonempty_index():
    """A second bulk batch with strictly greater intervals appends correctly."""
    pg = _pg()
    g = pg.BedGrove(8)
    g.insert_bulk("chr1", _items(pg, 10), presorted=True)            # 0..950
    g.insert_bulk("chr1", _items(pg, 10, start=5000), presorted=True)  # 5000..
    assert g.size() == 20
    assert len(list(g.intersect(pg.Interval(5000, 5050), "chr1"))) == 1


def test_bulk_on_gff_grove():
    """The fast paths exist on GffGrove too (same template)."""
    pg = _pg()
    g = pg.GffGrove(8)
    items = [(pg.Interval(i * 100, i * 100 + 50),
              pg.GffEntry("chr1", i * 100 + 1, i * 100 + 51, "gene"))
             for i in range(20)]
    keys = g.insert_bulk("chr1", items, presorted=True)
    assert len(keys) == 20
    assert g.size() == 20
    assert keys[0].data.type == "gene"


def test_dataless_grove_has_no_bulk():
    """Bulk/sorted insert require associated data — the dataless Grove lacks them."""
    pg = _pg()
    g = pg.Grove(8)
    assert not hasattr(g, "insert_bulk")
    assert not hasattr(g, "insert_sorted")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])