"""
Tests for the Numeric key value type (gdt::numeric).

Ports the behaviourally-observable cases from genogrove
tests/data_type/numeric_test.cpp: construction, comparisons, exact-equality
overlap, and to_string. Point semantics — overlap is equality, not range
intersection.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_parameterized_constructor():
    pg = _pg()
    assert pg.Numeric(42).value == 42
    assert pg.Numeric(-100).value == -100
    assert pg.Numeric(0).value == 0


def test_default_constructor_is_int_min():
    pg = _pg()
    # Default is the aggregation sentinel INT_MIN (32-bit).
    assert pg.Numeric().value == -(2**31)


def test_equality():
    pg = _pg()
    assert pg.Numeric(5) == pg.Numeric(5)
    assert not (pg.Numeric(5) == pg.Numeric(6))


def test_comparisons():
    pg = _pg()
    assert pg.Numeric(5) < pg.Numeric(10)
    assert pg.Numeric(10) > pg.Numeric(5)
    assert pg.Numeric(-5) < pg.Numeric(5)
    assert not (pg.Numeric(5) < pg.Numeric(5))


def test_overlap_is_exact_equality():
    pg = _pg()
    assert pg.Numeric.overlaps(pg.Numeric(5), pg.Numeric(5)) is True
    assert pg.Numeric.overlaps(pg.Numeric(5), pg.Numeric(6)) is False


def test_str_and_repr():
    pg = _pg()
    assert str(pg.Numeric(42)) == "42"
    assert repr(pg.Numeric(42)) == "Numeric(42)"
    assert str(pg.Numeric(-7)) == "-7"


def test_set_value_pre_insertion():
    pg = _pg()
    n = pg.Numeric(1)
    n.set_value(99)
    assert n.value == 99


def test_sorting_uses_value_order():
    pg = _pg()
    xs = [pg.Numeric(v) for v in (3, 1, 2, -5)]
    assert [n.value for n in sorted(xs)] == [-5, 1, 2, 3]