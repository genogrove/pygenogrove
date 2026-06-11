"""
Tests for the dataless grove (grove<interval>, i.e. grove<interval, void, void>).

Mirrors genogrove/tests/structure/dataless_grove_test.cpp (``DatalessGroveTest``):
insert/query, non-overlapping queries, node splits, and basic grove introspection
for the no-associated-data interval grove the Python bindings expose.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_creation():
    """Grove(order) reports its order and starts empty."""
    pg = _pg()
    grove = pg.Grove(100)
    assert grove.get_order() == 100
    assert grove.size() == 0


def test_default_order():
    """A default Grove uses order 3."""
    pg = _pg()
    grove = pg.Grove()
    assert grove.get_order() == 3


def test_insert_and_query():
    """Inserting an interval makes it discoverable and returns its Key."""
    pg = _pg()
    grove = pg.Grove(100)
    interval = pg.GenomicCoordinate(".", 100, 200)

    key = grove.insert("chr1", interval)

    assert grove.size() == 1
    assert key is not None
    assert key.value.start == 100
    assert key.value.end == 200


def test_multiple_inserts():
    """Inserting across indices accumulates the total size."""
    pg = _pg()
    grove = pg.Grove(100)

    grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    grove.insert("chr1", pg.GenomicCoordinate(".", 300, 400))
    grove.insert("chr2", pg.GenomicCoordinate(".", 100, 200))

    assert grove.size() == 3


def test_intersect_specific_index():
    """Querying a specific index returns only that index's overlaps."""
    pg = _pg()
    grove = pg.Grove(100)

    grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    grove.insert("chr1", pg.GenomicCoordinate(".", 300, 400))
    grove.insert("chr2", pg.GenomicCoordinate(".", 150, 250))

    query = pg.GenomicCoordinate(".", 150, 350)
    results = grove.intersect(query, "chr1")

    assert len(results) == 2


def test_intersect_all_indices():
    """Querying with no index searches every index."""
    pg = _pg()
    grove = pg.Grove(100)

    grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    grove.insert("chr2", pg.GenomicCoordinate(".", 150, 250))
    grove.insert("chr3", pg.GenomicCoordinate(".", 300, 400))

    query = pg.GenomicCoordinate(".", 175, 225)
    results = grove.intersect(query)

    # Should find overlaps in chr1 and chr2
    assert len(results) == 2


def test_empty_query():
    """Querying an empty grove returns no hits."""
    pg = _pg()
    grove = pg.Grove(100)
    query = pg.GenomicCoordinate(".", 100, 200)
    results = grove.intersect(query, "chr1")

    assert len(results) == 0


def test_non_overlapping_query_returns_empty():
    """A query disjoint from all stored intervals returns nothing."""
    pg = _pg()
    grove = pg.Grove(100)
    grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    grove.insert("chr1", pg.GenomicCoordinate(".", 300, 400))

    query = pg.GenomicCoordinate(".", 500, 600)
    results = grove.intersect(query, "chr1")

    assert len(results) == 0


def test_len_str_repr():
    """__len__ mirrors size(); __str__/__repr__ expose size and order."""
    pg = _pg()
    grove = pg.Grove(7)
    grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    grove.insert("chr1", pg.GenomicCoordinate(".", 300, 400))

    assert len(grove) == 2
    assert len(grove) == grove.size()
    assert "2" in str(grove)
    rep = repr(grove)
    assert "7" in rep and "2" in rep  # order and size


def test_node_splits_small_order():
    """Order 3 forces many node splits; size() must count only leaf keys
    (separator keys are not indexed vertices) and queries must find all."""
    pg = _pg()
    grove = pg.Grove(3)
    n = 50
    for i in range(n):
        grove.insert("chr1", pg.GenomicCoordinate(".", i * 100, i * 100 + 50))

    assert grove.size() == n
    assert grove.indexed_vertex_count() == n

    # A query spanning everything returns every inserted key.
    hits = grove.intersect(pg.GenomicCoordinate(".", 0, n * 100), "chr1")
    assert len(hits) == n


def test_intersect_adjacent_boundary():
    """intersect uses closed-interval overlap: a query touching a stored
    interval's boundary still matches."""
    pg = _pg()
    grove = pg.Grove(5)
    grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    # Query [200, 300] touches [100, 200] at 200.
    assert len(grove.intersect(pg.GenomicCoordinate(".", 200, 300), "chr1")) == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])