"""
Basic tests for pygenogrove functionality.

These tests verify the core interval and grove operations.
"""

import pytest


def test_imports():
    """Test that the module can be imported."""
    try:
        import pygenogrove
        assert hasattr(pygenogrove, 'Interval')
        assert hasattr(pygenogrove, 'Grove')
        assert hasattr(pygenogrove, 'Key')
        assert hasattr(pygenogrove, 'QueryResult')
    except ImportError as e:
        pytest.skip(f"Module not built yet: {e}")


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
