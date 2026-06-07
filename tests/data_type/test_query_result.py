"""
Tests for the QueryResult container (query_result<interval>).

Mirrors genogrove/tests/data_type/query_result_test.cpp (``query_result_test``),
restricted to the void-data surface exposed by the Python bindings.
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_constructor():
    """query/keys properties reflect the originating query and its hits."""
    pg = _pg()
    grove = pg.Grove(100)
    grove.insert("chr1", pg.Interval(100, 200))

    query = pg.Interval(150, 175)
    results = grove.intersect(query, "chr1")

    assert results.query.start == 150
    assert results.query.end == 175
    assert len(results.keys) == 1


def test_key_ordering_preserved():
    """Iterating a result yields Keys, len() matches the iterated count."""
    pg = _pg()
    grove = pg.Grove(100)
    grove.insert("chr1", pg.Interval(100, 200))
    grove.insert("chr1", pg.Interval(150, 250))

    query = pg.Interval(175, 225)
    results = grove.intersect(query, "chr1")

    count = 0
    for key in results:
        assert key is not None
        assert hasattr(key, 'value')
        count += 1

    assert count == len(results)


def test_empty_result():
    """An empty result has length 0 and yields nothing on iteration."""
    pg = _pg()
    grove = pg.Grove(3)
    grove.insert("chr1", pg.Interval(100, 200))

    results = grove.intersect(pg.Interval(500, 600), "chr1")
    assert len(results) == 0
    assert list(results) == []
    assert len(results.keys) == 0


def test_keys_keep_grove_alive():
    """Keys obtained from intersect() keep the grove alive (intersect uses
    keep_alive<0,1>). Using them after the only grove handle is dropped must
    stay safe — a regression here would be a use-after-free (crash), not a
    failed assertion.
    """
    pg = _pg()
    grove = pg.Grove(3)
    grove.insert("chr1", pg.Interval(100, 200))

    keys = list(grove.intersect(pg.Interval(150, 175), "chr1"))
    assert len(keys) == 1

    del grove           # drop the only Python handle to the grove
    gc.collect()

    # The query-result key chain kept the grove alive.
    assert keys[0].value.start == 100
    assert keys[0].value.end == 200


if __name__ == "__main__":
    pytest.main([__file__, "-v"])