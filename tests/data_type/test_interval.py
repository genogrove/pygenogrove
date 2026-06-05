"""
Tests for the Interval data type.

Mirrors genogrove/tests/data_type/interval_test.cpp (``intervalTest``),
restricted to the surface currently exposed by the Python bindings.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_default_constructor():
    """Default-constructed interval is usable; endpoints are integers."""
    pg = _pg()
    iv = pg.Interval()
    assert isinstance(iv.start, int)
    assert isinstance(iv.end, int)


def test_parameterized_constructor():
    """Interval(start, end) stores closed [start, end] endpoints."""
    pg = _pg()
    interval = pg.Interval(100, 200)
    assert interval.start == 100
    assert interval.end == 200


def test_constructor_rejects_inverted_interval():
    """end < start is rejected at construction (C++ throws invalid_argument)."""
    pg = _pg()
    with pytest.raises(ValueError):
        pg.Interval(200, 100)
    with pytest.raises(ValueError):
        pg.Interval(10, 5)


def test_constructor_accepts_equal_start_end():
    """A zero-length interval [p, p] is valid."""
    pg = _pg()
    iv = pg.Interval(100, 100)
    assert iv.start == 100
    assert iv.end == 100


def test_equality_operator():
    """== compares both endpoints."""
    pg = _pg()
    assert pg.Interval(20, 30) == pg.Interval(20, 30)
    assert not (pg.Interval(20, 30) == pg.Interval(20, 40))


def test_comparison_operators():
    """<, >, == sort by start, then by end."""
    pg = _pg()
    assert pg.Interval(10, 20) < pg.Interval(20, 30)   # by start
    assert pg.Interval(10, 20) < pg.Interval(10, 25)   # equal start, by end
    assert pg.Interval(20, 30) > pg.Interval(10, 20)
    assert pg.Interval(20, 30) == pg.Interval(20, 30)
    assert not (pg.Interval(20, 30) == pg.Interval(20, 40))


def test_overlap_overlapping():
    """Test interval overlap detection."""
    pg = _pg()
    interval1 = pg.Interval(100, 200)
    interval2 = pg.Interval(150, 250)
    interval3 = pg.Interval(300, 400)

    assert pg.Interval.overlaps(interval1, interval2)
    assert not pg.Interval.overlaps(interval1, interval3)


def test_overlap_adjacent():
    """Closed intervals that touch at a boundary overlap: [10,20] vs [20,30]."""
    pg = _pg()
    a = pg.Interval(10, 20)
    b = pg.Interval(20, 30)
    assert pg.Interval.overlaps(a, b)
    assert pg.Interval.overlaps(b, a)


def test_overlap_disjoint():
    """[10,20] and [21,30] do not overlap (closed-interval semantics)."""
    pg = _pg()
    assert not pg.Interval.overlaps(pg.Interval(10, 20), pg.Interval(21, 30))


def test_overlap_contained_and_identical():
    """Containment and identity both count as overlap."""
    pg = _pg()
    outer, inner = pg.Interval(10, 50), pg.Interval(20, 30)
    assert pg.Interval.overlaps(outer, inner)
    assert pg.Interval.overlaps(inner, outer)
    same = pg.Interval(10, 30)
    assert pg.Interval.overlaps(same, pg.Interval(10, 30))


def test_to_string():
    """str/repr expose the endpoints; repr round-trips the constructor form."""
    pg = _pg()
    interval = pg.Interval(100, 200)
    str_repr = str(interval)
    assert "100" in str_repr
    assert "200" in str_repr
    assert repr(pg.Interval(100, 200)) == "Interval(100, 200)"


def test_getters_and_setters():
    """start/end are read-only; set_range is the only mutation path (pre-insertion).

    Mutating an inserted interval would corrupt B+ tree ordering, so the
    endpoint properties are intentionally read-only.
    """
    pg = _pg()
    interval = pg.Interval(100, 200)
    with pytest.raises(AttributeError):
        interval.start = 150
    with pytest.raises(AttributeError):
        interval.end = 250

    interval.set_range(150, 250)
    assert interval.start == 150
    assert interval.end == 250


def test_large_coordinate_values():
    """Large (genome-scale) coordinates round-trip and overlap correctly."""
    pg = _pg()
    big = pg.Interval(3_000_000_000, 3_000_000_500)
    assert big.start == 3_000_000_000
    assert big.end == 3_000_000_500
    assert pg.Interval.overlaps(big, pg.Interval(3_000_000_250, 3_000_001_000))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])