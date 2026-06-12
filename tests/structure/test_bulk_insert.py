"""
Tests for the sorted / bulk insertion fast paths on the data-carrying groves
(BedGrove / GffGrove): insert_sorted and insert_bulk (explicit Interval keys).

Mirrors the genogrove key_type_grove_test harness (assert_sorted_insert,
assert_bulk_insert, assert_bulk_insert_returns_handles, assert_edge_creation).
These fast paths require associated data, so they exist on BedGrove/GffGrove but
not the dataless Grove.

Pairing is verified by a coordinate-free identity in the payload (``name``), not
by coordinate equality — the Interval key and the BedEntry payload use different
coordinate conventions (closed vs half-open), so a key carries "its own" data
iff the payload's name matches the record it was inserted with.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _items(pg, n, start=0, step=100):
    """n ascending (Interval, BedEntry) records, each tagged name=f"f{i}".

    The BedEntry's half-open [s, s+51) covers the same region as the closed
    Interval [s, s+50] — honest coordinates — but pairing is checked via name.
    """
    out = []
    for i in range(n):
        s = start + i * step
        e = pg.BedEntry("chr1", s, s + 51)
        e.name = f"f{i}"
        out.append((pg.GenomicCoordinate(".", s, s + 50), e))
    return out


def test_insert_sorted_single():
    """insert_sorted appends ascending records; data stays paired. Mirrors assert_sorted_insert."""
    pg = _pg()
    g = pg.BedGrove(8)
    keys = [g.insert_sorted("chr1", iv, d) for iv, d in _items(pg, 30)]

    assert len(keys) == 30
    assert g.size() == 30
    assert [k.data.name for k in keys] == [f"f{i}" for i in range(30)]
    assert len(list(g.intersect(pg.GenomicCoordinate(".", 0, 30 * 100), "chr1"))) == 30


def test_insert_bulk_presorted():
    """insert_bulk(presorted=True) inserts a sorted batch; handles map 1:1.

    n=50 at order 8 spans many leaves, exercising multi-leaf bottom-up build.
    """
    pg = _pg()
    g = pg.BedGrove(8)
    keys = g.insert_bulk("chr1", _items(pg, 50), presorted=True)

    assert len(keys) == 50
    assert g.size() == 50
    assert [k.value.start for k in keys] == [i * 100 for i in range(50)]
    assert [k.data.name for k in keys] == [f"f{i}" for i in range(50)]
    assert len(list(g.intersect(pg.GenomicCoordinate(".", 0, 50 * 100), "chr1"))) == 50


def test_insert_bulk_sorts_unsorted_input():
    """insert_bulk (presorted=False) sorts the batch AND keeps each datum with
    its interval. Feeds descending input so the internal sort must reorder pairs;
    the name check catches a pairing scramble (intervals-sorted alone would not)."""
    pg = _pg()
    g = pg.BedGrove(8)
    items = _items(pg, 40)
    keys = g.insert_bulk("chr1", list(reversed(items)))   # descending input

    assert len(keys) == 40
    assert g.size() == 40
    assert [k.value.start for k in keys] == [i * 100 for i in range(40)]
    # crucially: each interval still carries ITS OWN data after the sort
    assert [k.data.name for k in keys] == [f"f{i}" for i in range(40)]
    assert len(list(g.intersect(pg.GenomicCoordinate(".", 0, 40 * 100), "chr1"))) == 40


def test_presorted_flag_is_honored():
    """presorted changes behavior: handles come back in input order when
    presorted=True (trusts the caller) vs sorted order when False.

    The presorted=True grove is fed descending data on purpose (an intentional
    precondition violation); its tree is left ill-ordered, so it is NOT queried —
    only the well-defined returned-handle order is asserted.
    """
    pg = _pg()
    descending = list(reversed(_items(pg, 20)))
    input_starts = [iv.start for iv, _ in descending]            # 1900, 1800, ...

    g_sorted = pg.BedGrove(8)
    k_sorted = g_sorted.insert_bulk("chr1", descending)          # presorted=False -> sorts
    assert [k.value.start for k in k_sorted] == sorted(input_starts)

    g_trust = pg.BedGrove(8)
    k_trust = g_trust.insert_bulk("chr1", descending, presorted=True)  # trusts input order
    assert [k.value.start for k in k_trust] == input_starts
    assert [k.value.start for k in k_sorted] != [k.value.start for k in k_trust]


def test_bulk_returns_handles_with_correct_data():
    """Handles map 1:1 to inserted records. Mirrors assert_bulk_insert_returns_handles."""
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
    """Bulk insert yields the same query CONTENT as one-by-one insert()."""
    pg = _pg()
    items = _items(pg, 40)

    g_one = pg.BedGrove(8)
    for iv, d in items:
        g_one.insert("chr1", iv, d)

    g_bulk = pg.BedGrove(8)
    g_bulk.insert_bulk("chr1", items, presorted=True)

    assert g_one.size() == g_bulk.size() == 40
    q = pg.GenomicCoordinate(".", 550, 2550)
    names_one = sorted(k.data.name for k in g_one.intersect(q, "chr1"))
    names_bulk = sorted(k.data.name for k in g_bulk.intersect(q, "chr1"))
    assert names_one == names_bulk
    assert len(names_bulk) > 0


def test_bulk_empty_batch():
    """An empty batch inserts nothing and returns no handles."""
    pg = _pg()
    g = pg.BedGrove(8)
    assert g.insert_bulk("chr1", []) == []
    assert g.size() == 0


def test_bulk_append_to_nonempty_index():
    """A second bulk batch with strictly greater intervals appends correctly,
    exercising the rightmost-append path (not bottom-up), and the appended
    handles still carry their own data."""
    pg = _pg()
    g = pg.BedGrove(8)
    g.insert_bulk("chr1", _items(pg, 10), presorted=True)              # 0..950 (bottom-up)
    batch2 = g.insert_bulk("chr1", _items(pg, 10, start=5000), presorted=True)  # append path

    assert g.size() == 20
    assert [k.value.start for k in batch2] == [5000 + i * 100 for i in range(10)]
    assert [k.data.name for k in batch2] == [f"f{i}" for i in range(10)]
    assert len(list(g.intersect(pg.GenomicCoordinate(".", 5000, 5050), "chr1"))) == 1
    # the first batch survives the append
    assert len(list(g.intersect(pg.GenomicCoordinate(".", 0, 950), "chr1"))) == 10


def test_bulk_on_gff_grove():
    """The fast paths exist on GffGrove too; data stays paired (same template)."""
    pg = _pg()
    g = pg.GffGrove(8)
    items = []
    for i in range(20):
        e = pg.GffEntry("chr1", i * 100 + 1, i * 100 + 51, "gene")
        e.attributes = {"idx": str(i)}
        items.append((pg.GenomicCoordinate(".", i * 100, i * 100 + 50), e))
    keys = g.insert_bulk("chr1", items, presorted=True)

    assert len(keys) == 20
    assert g.size() == 20
    assert [k.value.start for k in keys] == [i * 100 for i in range(20)]
    assert [k.data.get_attribute("idx") for k in keys] == [str(i) for i in range(20)]


def test_grove_supports_bulk_with_json_payload():
    """The universal Grove carries data, so it has sorted/bulk insert too."""
    pg = _pg()
    g = pg.Grove(8)
    assert hasattr(g, "insert_bulk")
    assert hasattr(g, "insert_sorted")

    items = [(pg.GenomicCoordinate(".", i * 10, i * 10 + 5), {"i": i}) for i in range(5)]
    keys = g.insert_bulk("chr1", items)
    assert len(keys) == 5
    res = g.intersect(pg.GenomicCoordinate(".", 0, 1000), "chr1")
    assert sorted(k.data["i"] for k in res) == [0, 1, 2, 3, 4]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])