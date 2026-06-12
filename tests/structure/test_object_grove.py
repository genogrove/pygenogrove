"""
Tests for the universal Grove's JSON payload: grove<genomic_coordinate, json_value>
stores any JSON-serializable Python object (dict / list / scalar / None) as the
per-key data, round-tripping it transparently (no user-facing json import).

This is the headline behaviour of the redesigned, genomic-coordinate-standard
Grove. The typed BedGrove/GffGrove (tested elsewhere) are the schema'd
alternative for full C++ interop.
"""

import gc as _gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _c(pg, strand, start, end):
    return pg.GenomicCoordinate(strand, start, end)


def test_dict_payload_roundtrips():
    pg = _pg()
    g = pg.Grove()
    payload = {"gene": "FOO", "score": 5, "exons": [1, 2, 3], "meta": {"x": True}}
    k = g.insert("chr1", _c(pg, "+", 100, 200), payload)
    assert k.data == payload
    hit = list(g.intersect(_c(pg, "+", 150, 160), "chr1"))[0]
    assert hit.data == payload


def test_omitted_data_defaults_to_none():
    """insert without a data argument stores None (data is optional)."""
    pg = _pg()
    g = pg.Grove()
    k = g.insert("chr1", _c(pg, "+", 100, 200))  # no data arg
    assert k.data is None
    assert g.size() == 1


def test_explicit_none_payload():
    pg = _pg()
    g = pg.Grove()
    k = g.insert("chr1", _c(pg, "+", 100, 200), None)
    assert k.data is None


def test_heterogeneous_payloads_in_one_grove():
    """Each key may carry a different JSON shape — no schema is enforced."""
    pg = _pg()
    g = pg.Grove()
    g.insert("chr1", _c(pg, "+", 100, 200), {"gene": "FOO"})
    g.insert("chr1", _c(pg, "+", 300, 400), [1, 2, 3])
    g.insert("chr1", _c(pg, "-", 500, 600), "label")
    g.insert("chr1", _c(pg, ".", 700, 800), 42)

    by_start = {k.value.start: k.data for k in g.intersect(_c(pg, "*", 0, 10**6), "chr1")}
    assert by_start[100] == {"gene": "FOO"}
    assert by_start[300] == [1, 2, 3]
    assert by_start[500] == "label"
    assert by_start[700] == 42


def test_non_json_serializable_raises():
    pg = _pg()
    g = pg.Grove()
    with pytest.raises((TypeError, ValueError)):
        g.insert("chr1", _c(pg, "+", 100, 200), {1, 2, 3})  # a set is not JSON


def test_strand_aware_query_with_payload():
    pg = _pg()
    g = pg.Grove()
    g.insert("chr1", _c(pg, "+", 100, 200), {"name": "plus"})
    g.insert("chr1", _c(pg, "-", 100, 200), {"name": "minus"})
    res = g.intersect(_c(pg, "+", 150, 160), "chr1")
    assert [k.data["name"] for k in res] == ["plus"]


def test_serialization_roundtrip_preserves_strand_and_payload(tmp_path):
    pg = _pg()
    g = pg.Grove()
    g.insert("chr1", _c(pg, "+", 100, 200), {"gene": "FOO"})
    g.insert("chr1", _c(pg, "-", 300, 400), {"gene": "BAR", "n": 7})
    g.insert("chr2", _c(pg, ".", 150, 250), None)

    path = str(tmp_path / "obj.gg")
    g.serialize(path)
    loaded = pg.Grove.deserialize(path)

    assert loaded.size() == 3
    res = loaded.intersect(_c(pg, "+", 150, 160), "chr1")
    assert [(k.value.strand, k.data) for k in res] == [("+", {"gene": "FOO"})]


def test_keys_keep_grove_alive_with_payload():
    pg = _pg()
    g = pg.Grove()
    g.insert("chr1", _c(pg, "+", 100, 200), {"gene": "FOO"})
    keys = list(g.intersect(_c(pg, "+", 150, 160), "chr1"))
    del g
    _gc.collect()
    assert keys[0].data == {"gene": "FOO"}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])