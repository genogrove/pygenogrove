"""
Basic tests for pygenogrove functionality.

These tests verify the core interval and grove operations.
"""

import pytest


def test_imports():
    """Test that the module can be imported.

    This is a hard failure (not a skip) on purpose: if the compiled extension
    is not installed/importable, the rest of the suite would otherwise pass by
    silently skipping every test, hiding a broken build or packaging.
    """
    import pygenogrove
    assert hasattr(pygenogrove, 'Interval')
    assert hasattr(pygenogrove, 'Grove')
    assert hasattr(pygenogrove, 'Key')
    assert hasattr(pygenogrove, 'QueryResult')


class TestInterval:
    """Tests for the Interval class."""

    def test_interval_creation(self):
        """Test creating an interval."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        interval = pg.Interval(100, 200)
        assert interval.start == 100
        assert interval.end == 200

    def test_interval_properties_readonly(self):
        """start/end are read-only. Mutating an inserted interval would corrupt B+ tree ordering."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        interval = pg.Interval(100, 200)
        with pytest.raises(AttributeError):
            interval.start = 150
        with pytest.raises(AttributeError):
            interval.end = 250

    def test_interval_set_range(self):
        """set_range is the only path to mutate an interval (use before insertion)."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        interval = pg.Interval(100, 200)
        interval.set_range(150, 250)
        assert interval.start == 150
        assert interval.end == 250

    def test_interval_overlap(self):
        """Test interval overlap detection."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        interval1 = pg.Interval(100, 200)
        interval2 = pg.Interval(150, 250)
        interval3 = pg.Interval(300, 400)

        # Test overlapping intervals
        assert pg.Interval.overlaps(interval1, interval2)

        # Test non-overlapping intervals
        assert not pg.Interval.overlaps(interval1, interval3)

    def test_interval_string_representation(self):
        """Test interval string conversion."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        interval = pg.Interval(100, 200)
        str_repr = str(interval)
        assert "100" in str_repr
        assert "200" in str_repr

    def test_interval_repr(self):
        """repr should round-trip the constructor form."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        assert repr(pg.Interval(100, 200)) == "Interval(100, 200)"

    def test_interval_default_constructor(self):
        """Default-constructed interval is usable."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        # Should not raise; endpoints are accessible integers.
        iv = pg.Interval()
        assert isinstance(iv.start, int)
        assert isinstance(iv.end, int)

    def test_interval_rejects_inverted(self):
        """end < start is rejected at construction (C++ throws invalid_argument)."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        with pytest.raises(ValueError):
            pg.Interval(200, 100)
        with pytest.raises(ValueError):
            pg.Interval(10, 5)

    def test_interval_accepts_equal_start_end(self):
        """A zero-length interval [p, p] is valid."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        iv = pg.Interval(100, 100)
        assert iv.start == 100
        assert iv.end == 100

    def test_interval_overlap_adjacent_boundary(self):
        """Closed intervals that touch at a boundary overlap: [10,20] vs [20,30]."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        a = pg.Interval(10, 20)
        b = pg.Interval(20, 30)
        assert pg.Interval.overlaps(a, b)
        assert pg.Interval.overlaps(b, a)

    def test_interval_overlap_disjoint_by_one(self):
        """[10,20] and [21,30] do not overlap (closed-interval semantics)."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        assert not pg.Interval.overlaps(pg.Interval(10, 20), pg.Interval(21, 30))

    def test_interval_overlap_contained_and_identical(self):
        """Containment and identity both count as overlap."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        outer, inner = pg.Interval(10, 50), pg.Interval(20, 30)
        assert pg.Interval.overlaps(outer, inner)
        assert pg.Interval.overlaps(inner, outer)
        same = pg.Interval(10, 30)
        assert pg.Interval.overlaps(same, pg.Interval(10, 30))

    def test_interval_comparison_operators(self):
        """<, >, == sort by start, then by end."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        assert pg.Interval(10, 20) < pg.Interval(20, 30)   # by start
        assert pg.Interval(10, 20) < pg.Interval(10, 25)   # equal start, by end
        assert pg.Interval(20, 30) > pg.Interval(10, 20)
        assert pg.Interval(20, 30) == pg.Interval(20, 30)
        assert not (pg.Interval(20, 30) == pg.Interval(20, 40))

    def test_interval_large_coordinates(self):
        """Large (genome-scale) coordinates round-trip and overlap correctly."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        big = pg.Interval(3_000_000_000, 3_000_000_500)
        assert big.start == 3_000_000_000
        assert big.end == 3_000_000_500
        assert pg.Interval.overlaps(big, pg.Interval(3_000_000_250, 3_000_001_000))


class TestGrove:
    """Tests for the Grove class."""

    def test_grove_creation(self):
        """Test creating a grove."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)
        assert grove.get_order() == 100
        assert grove.size() == 0

    def test_grove_default_order(self):
        """Test grove with default order."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove()
        assert grove.get_order() == 3

    def test_grove_insert(self):
        """Test inserting intervals into a grove."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)
        interval = pg.Interval(100, 200)

        key = grove.insert("chr1", interval)

        assert grove.size() == 1
        assert key is not None
        assert key.value.start == 100
        assert key.value.end == 200

    def test_grove_multiple_inserts(self):
        """Test inserting multiple intervals."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)

        grove.insert("chr1", pg.Interval(100, 200))
        grove.insert("chr1", pg.Interval(300, 400))
        grove.insert("chr2", pg.Interval(100, 200))

        assert grove.size() == 3

    def test_grove_intersect_specific_index(self):
        """Test querying intervals in a specific index."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)

        grove.insert("chr1", pg.Interval(100, 200))
        grove.insert("chr1", pg.Interval(300, 400))
        grove.insert("chr2", pg.Interval(150, 250))

        # Query chr1
        query = pg.Interval(150, 350)
        results = grove.intersect(query, "chr1")

        assert len(results) == 2

    def test_grove_intersect_all_indices(self):
        """Test querying intervals across all indices."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)

        grove.insert("chr1", pg.Interval(100, 200))
        grove.insert("chr2", pg.Interval(150, 250))
        grove.insert("chr3", pg.Interval(300, 400))

        # Query all chromosomes
        query = pg.Interval(175, 225)
        results = grove.intersect(query)

        # Should find overlaps in chr1 and chr2
        assert len(results) == 2

    def test_grove_empty_query(self):
        """Test querying an empty grove."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)
        query = pg.Interval(100, 200)
        results = grove.intersect(query, "chr1")

        assert len(results) == 0

    def test_grove_no_overlap_query(self):
        """Test querying with no overlapping intervals."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)
        grove.insert("chr1", pg.Interval(100, 200))
        grove.insert("chr1", pg.Interval(300, 400))

        query = pg.Interval(500, 600)
        results = grove.intersect(query, "chr1")

        assert len(results) == 0

    def test_grove_len_str_repr(self):
        """__len__ mirrors size(); __str__/__repr__ expose size and order."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(7)
        grove.insert("chr1", pg.Interval(100, 200))
        grove.insert("chr1", pg.Interval(300, 400))

        assert len(grove) == 2
        assert len(grove) == grove.size()
        assert "2" in str(grove)
        rep = repr(grove)
        assert "7" in rep and "2" in rep  # order and size

    def test_grove_node_splits_small_order(self):
        """Order 3 forces many node splits; size() must count only leaf keys
        (separator keys are not indexed vertices) and queries must find all."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(3)
        n = 50
        for i in range(n):
            grove.insert("chr1", pg.Interval(i * 100, i * 100 + 50))

        assert grove.size() == n
        assert grove.indexed_vertex_count() == n

        # A query spanning everything returns every inserted key.
        hits = grove.intersect(pg.Interval(0, n * 100), "chr1")
        assert len(hits) == n

    def test_grove_intersect_adjacent_boundary(self):
        """intersect uses closed-interval overlap: a query touching a stored
        interval's boundary still matches."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(5)
        grove.insert("chr1", pg.Interval(100, 200))
        # Query [200, 300] touches [100, 200] at 200.
        assert len(grove.intersect(pg.Interval(200, 300), "chr1")) == 1


class TestKey:
    """Tests for the Key wrapper."""

    def test_key_str(self):
        """Key has a string representation."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(3)
        key = grove.insert("chr1", pg.Interval(100, 200))
        assert isinstance(str(key), str)
        assert str(key) != ""

    def test_key_value_is_a_copy(self):
        """key.value returns a copy; mutating it must not affect stored ordering."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(3)
        key = grove.insert("chr1", pg.Interval(100, 200))

        snapshot = key.value
        snapshot.set_range(0, 5)          # mutate the returned copy
        assert key.value.start == 100     # stored key is unchanged
        assert key.value.end == 200


class TestQueryResult:
    """Tests for QueryResult class."""

    def test_query_result_iteration(self):
        """Test iterating over query results."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)
        grove.insert("chr1", pg.Interval(100, 200))
        grove.insert("chr1", pg.Interval(150, 250))

        query = pg.Interval(175, 225)
        results = grove.intersect(query, "chr1")

        # Test iteration
        count = 0
        for key in results:
            assert key is not None
            assert hasattr(key, 'value')
            count += 1

        assert count == len(results)

    def test_query_result_properties(self):
        """Test query result properties."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(100)
        grove.insert("chr1", pg.Interval(100, 200))

        query = pg.Interval(150, 175)
        results = grove.intersect(query, "chr1")

        # Check query property
        assert results.query.start == 150
        assert results.query.end == 175

        # Check keys property
        assert len(results.keys) == 1

    def test_query_result_empty(self):
        """An empty result has length 0 and yields nothing on iteration."""
        pytest.importorskip("pygenogrove")
        import pygenogrove as pg

        grove = pg.Grove(3)
        grove.insert("chr1", pg.Interval(100, 200))

        results = grove.intersect(pg.Interval(500, 600), "chr1")
        assert len(results) == 0
        assert list(results) == []
        assert len(results.keys) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
