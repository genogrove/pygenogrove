"""
Tests for the Key wrapper (key<interval>).

Mirrors genogrove/tests/data_type/key_test.cpp (``keyTest``). The Python Key
wraps a pointer into grove storage, so the binding-specific concerns are its
string form, value-copy semantics, and that it keeps its Grove alive.
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_str():
    """Key has a non-empty string representation."""
    pg = _pg()
    grove = pg.Grove(3)
    key = grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    assert isinstance(str(key), str)
    assert str(key) != ""


def test_value_is_a_copy():
    """key.value returns a copy; mutating it must not affect stored ordering."""
    pg = _pg()
    grove = pg.Grove(3)
    key = grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))

    snapshot = key.value
    snapshot.set_range(0, 5)          # mutate the returned copy
    assert key.value.start == 100     # stored key is unchanged
    assert key.value.end == 200


def test_repr():
    pg = _pg()
    grove = pg.Grove(3)
    key = grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    assert repr(key) == f"Key({key.value})"


def test_equality_and_hash_by_value_ignoring_data():
    """Comparisons delegate to the wrapped value, ignoring .data and pointer
    identity — two distinct Keys (different Groves, different data) with equal
    values must compare equal and hash equal."""
    pg = _pg()
    a = pg.Grove(3).insert("chr1", pg.GenomicCoordinate(".", 100, 200), {"x": 1})
    b = pg.Grove(3).insert("chr1", pg.GenomicCoordinate(".", 100, 200), {"x": 2})
    assert a is not b
    assert a == b
    assert hash(a) == hash(b)


def test_ordering_by_value():
    pg = _pg()
    grove = pg.Grove(3)
    lo = grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    hi = grove.insert("chr1", pg.GenomicCoordinate(".", 300, 400))
    assert lo < hi
    assert hi > lo


def test_keeps_grove_alive():
    """A Key holds its Grove alive: using it after the Grove handle is
    dropped must stay safe (reference_internal => keep_alive<0,1>).
    If this contract regressed it would be a use-after-free (crash),
    not a failed assertion."""
    pg = _pg()
    grove = pg.Grove(3)
    key = grove.insert("chr1", pg.GenomicCoordinate(".", 100, 200))

    del grove        # drop the only Python handle to the Grove
    gc.collect()     # ...and force collection

    # The Grove is kept alive by `key`, so this must not crash.
    assert key.value.start == 100
    assert key.value.end == 200


if __name__ == "__main__":
    pytest.main([__file__, "-v"])