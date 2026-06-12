"""
Tests for GffKey (key<interval, gff_entry>) — the value/data accessors and
lifetime semantics of a data-carrying key.

Mirrors tests/data_type/test_key.py (the dataless Key) plus the associated
GffEntry payload. Keys are obtained from a GffGrove because that is the only
way to construct one; the assertions here are about the key, not the grove.
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _key(pg):
    g = pg.GffGrove(3)
    return g, g.insert("chr1", pg.GenomicCoordinate(".", 100, 200),
                       pg.GffEntry("chr1", 100, 200, "gene"))


def test_str():
    """GffKey has a non-empty string representation."""
    pg = _pg()
    _g, key = _key(pg)
    assert isinstance(str(key), str)
    assert str(key) != ""


def test_value_is_a_copy():
    """key.value returns a copy; mutating it can't corrupt tree ordering."""
    pg = _pg()
    _g, key = _key(pg)

    snapshot = key.value
    snapshot.set_range(0, 5)
    assert key.value.start == 100
    assert key.value.end == 200


def test_data_is_a_live_mutable_reference():
    """key.data is a live reference: mutating it (incl. attributes) persists.

    Unlike .value, the data payload is not part of the B+ tree ordering, so
    in-place mutation is safe and intended.
    """
    pg = _pg()
    _g, key = _key(pg)

    key.data.source = "ensembl"
    key.data.attributes = {"gene_id": "ENSG2"}
    assert key.data.source == "ensembl"
    assert key.data.get_gene_id() == "ENSG2"


def test_keeps_grove_alive():
    """A GffKey keeps its GffGrove alive (reference_internal => keep_alive).

    A regression here would be a use-after-free (crash), not a failed assert.
    """
    pg = _pg()
    g = pg.GffGrove(3)
    key = g.insert("chr1", pg.GenomicCoordinate(".", 100, 200),
                   pg.GffEntry("chr1", 100, 200, "gene"))

    del g
    gc.collect()

    assert key.value.start == 100
    assert key.data.type == "gene"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
