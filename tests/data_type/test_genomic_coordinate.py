"""
Tests for the GenomicCoordinate value type (gdt::genomic_coordinate).

Mirrors genogrove/tests/data_type/genomic_coordinate_test.cpp over the surface
the Python bindings expose: construction + validation, strand-aware overlap,
coordinate-first ordering, getters/setters, and to_string.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_default_constructor():
    """Default coordinate is strand '.', [0, 0]."""
    pg = _pg()
    c = pg.GenomicCoordinate()
    assert c.strand == "."
    assert c.start == 0
    assert c.end == 0


def test_parameterized_constructor():
    pg = _pg()
    c = pg.GenomicCoordinate("+", 100, 200)
    assert c.strand == "+"
    assert c.start == 100
    assert c.end == 200


def test_constructor_rejects_invalid_strand():
    pg = _pg()
    with pytest.raises(ValueError):
        pg.GenomicCoordinate("x", 100, 200)


def test_constructor_rejects_inverted_interval():
    pg = _pg()
    with pytest.raises(ValueError):
        pg.GenomicCoordinate("+", 200, 100)


def test_constructor_accepts_all_valid_strands():
    pg = _pg()
    for strand in ("+", "-", ".", "*"):
        c = pg.GenomicCoordinate(strand, 10, 20)
        assert c.strand == strand


def test_equality_requires_all_three_components():
    pg = _pg()
    a = pg.GenomicCoordinate("+", 100, 200)
    assert a == pg.GenomicCoordinate("+", 100, 200)
    assert not (a == pg.GenomicCoordinate("-", 100, 200))
    assert not (a == pg.GenomicCoordinate("+", 101, 200))
    assert not (a == pg.GenomicCoordinate("+", 100, 201))


def test_comparison_by_start():
    pg = _pg()
    assert pg.GenomicCoordinate("+", 100, 200) < pg.GenomicCoordinate("+", 150, 200)
    assert pg.GenomicCoordinate("+", 150, 200) > pg.GenomicCoordinate("+", 100, 200)


def test_comparison_by_end():
    """Equal start -> ordered by end."""
    pg = _pg()
    assert pg.GenomicCoordinate("+", 100, 200) < pg.GenomicCoordinate("+", 100, 250)


def test_comparison_by_strand():
    """Equal coordinates -> strand order is * < . < + < -."""
    pg = _pg()
    star = pg.GenomicCoordinate("*", 100, 200)
    dot = pg.GenomicCoordinate(".", 100, 200)
    plus = pg.GenomicCoordinate("+", 100, 200)
    minus = pg.GenomicCoordinate("-", 100, 200)
    assert star < dot < plus < minus


def test_overlap_same_strand_overlapping():
    pg = _pg()
    a = pg.GenomicCoordinate("+", 100, 200)
    b = pg.GenomicCoordinate("+", 150, 250)
    assert pg.GenomicCoordinate.overlaps(a, b)


def test_overlap_same_strand_touching_is_overlap():
    """Closed intervals that touch at a point overlap."""
    pg = _pg()
    a = pg.GenomicCoordinate("+", 100, 200)
    b = pg.GenomicCoordinate("+", 200, 300)
    assert pg.GenomicCoordinate.overlaps(a, b)


def test_overlap_same_strand_disjoint():
    pg = _pg()
    a = pg.GenomicCoordinate("+", 100, 200)
    b = pg.GenomicCoordinate("+", 201, 300)
    assert not pg.GenomicCoordinate.overlaps(a, b)


def test_overlap_different_strands_never_overlap():
    """Coordinate overlap but mismatched strands -> no overlap."""
    pg = _pg()
    a = pg.GenomicCoordinate("+", 100, 200)
    b = pg.GenomicCoordinate("-", 150, 250)
    assert not pg.GenomicCoordinate.overlaps(a, b)


def test_overlap_wildcard_matches_any_strand():
    pg = _pg()
    star = pg.GenomicCoordinate("*", 150, 250)
    assert pg.GenomicCoordinate.overlaps(pg.GenomicCoordinate("+", 100, 200), star)
    assert pg.GenomicCoordinate.overlaps(pg.GenomicCoordinate("-", 100, 200), star)
    assert pg.GenomicCoordinate.overlaps(star, pg.GenomicCoordinate(".", 100, 200))


def test_overlap_unstranded_only_matches_unstranded():
    """'.' is a concrete strand, not a wildcard — it only matches '.'."""
    pg = _pg()
    dot = pg.GenomicCoordinate(".", 150, 250)
    assert pg.GenomicCoordinate.overlaps(pg.GenomicCoordinate(".", 100, 200), dot)
    assert not pg.GenomicCoordinate.overlaps(pg.GenomicCoordinate("+", 100, 200), dot)


def test_to_string():
    """Format is 'strand:start-end'."""
    pg = _pg()
    assert str(pg.GenomicCoordinate("+", 100, 200)) == "+:100-200"
    assert str(pg.GenomicCoordinate(".", 0, 0)) == ".:0-0"


def test_repr_roundtrips_shape():
    pg = _pg()
    assert repr(pg.GenomicCoordinate("-", 5, 9)) == "GenomicCoordinate('-', 5, 9)"


def test_set_range_and_set_strand():
    """Mutators work pre-insertion (not on a stored key)."""
    pg = _pg()
    c = pg.GenomicCoordinate("+", 100, 200)
    c.set_range(300, 400)
    assert c.start == 300 and c.end == 400
    c.set_strand("-")
    assert c.strand == "-"


def test_set_range_rejects_inverted():
    pg = _pg()
    c = pg.GenomicCoordinate("+", 100, 200)
    with pytest.raises(ValueError):
        c.set_range(400, 300)


def test_set_strand_rejects_invalid():
    pg = _pg()
    c = pg.GenomicCoordinate("+", 100, 200)
    with pytest.raises(ValueError):
        c.set_strand("?")


def test_sorting_mixed():
    """Coordinate-first sort: start, then end, then strand (* < . < + < -)."""
    pg = _pg()
    coords = [
        pg.GenomicCoordinate("-", 100, 200),
        pg.GenomicCoordinate("+", 100, 200),
        pg.GenomicCoordinate("+", 50, 200),
        pg.GenomicCoordinate("*", 100, 200),
    ]
    ordered = sorted(coords)
    assert [str(c) for c in ordered] == ["+:50-200", "*:100-200", "+:100-200", "-:100-200"]