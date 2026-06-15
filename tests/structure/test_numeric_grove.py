"""
Behaviour of NumericGrove (grove<numeric, json_value, json_value>): a B+ tree
over integer point keys. Overlap is exact equality, so intersect() is a point
lookup. Exercised dataless (data defaults to None) and with a JSON payload.

Mirrors genogrove/tests/structure/numeric_grove_test.cpp over the bound surface,
plus one split-stress insert and a serialization round-trip.
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_creation_and_insert():
    pg = _pg()
    g = pg.NumericGrove(8)
    assert g.get_order() == 8
    assert g.size() == 0

    key = g.insert("ids", pg.Numeric(42))
    assert g.size() == 1
    assert key.value == pg.Numeric(42)
    assert key.value.value == 42
    assert key.data is None  # dataless insert


def test_intersect_is_exact_point_lookup():
    pg = _pg()
    g = pg.NumericGrove(8)
    for v in (3, 7, 5, 42, -1):
        g.insert("ids", pg.Numeric(v))

    hits = list(g.intersect(pg.Numeric(7), "ids"))
    assert [k.value.value for k in hits] == [7]

    assert len(g.intersect(pg.Numeric(8), "ids")) == 0  # no equal key


def test_json_payload_roundtrips():
    pg = _pg()
    g = pg.NumericGrove(8)
    g.insert("ids", pg.Numeric(10), {"label": "ten"})
    hit = list(g.intersect(pg.Numeric(10), "ids"))[0]
    assert hit.data == {"label": "ten"}


def test_many_inserts_force_splits():
    pg = _pg()
    g = pg.NumericGrove(3)  # small order -> many node splits
    for v in range(100):
        g.insert("ids", pg.Numeric(v))
    assert g.size() == 100
    # Every point is individually findable after rebalancing.
    for v in (0, 1, 50, 98, 99):
        assert [k.value.value for k in g.intersect(pg.Numeric(v), "ids")] == [v]


def test_serialization_roundtrip(tmp_path):
    pg = _pg()
    g = pg.NumericGrove(4)
    for v in (1, 2, 3):
        g.insert("ids", pg.Numeric(v), {"v": v})
    path = str(tmp_path / "numeric.gg")
    g.serialize(path)

    loaded = pg.NumericGrove.deserialize(path)
    assert loaded.size() == 3
    hit = list(loaded.intersect(pg.Numeric(2), "ids"))[0]
    assert hit.value.value == 2
    assert hit.data == {"v": 2}


def test_keys_keep_grove_alive():
    """A NumericKey keeps its grove alive (use-after-free guard)."""
    pg = _pg()
    g = pg.NumericGrove(3)
    key = g.insert("ids", pg.Numeric(5))
    del g
    gc.collect()
    assert key.value.value == 5