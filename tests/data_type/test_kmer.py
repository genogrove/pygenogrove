"""
Tests for the Kmer key value type (gdt::kmer).

Ports the behaviourally-observable cases from genogrove
tests/data_type/kmer_test.cpp: sequence <-> 2-bit encoding round-trip, exact
equality overlap (same bases AND same length), ordering, validation, and the
k <= 32 / canonical-base constraints. Point semantics — overlap is equality.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_sequence_roundtrips_via_to_string():
    pg = _pg()
    for seq in ("ACGT", "AAAA", "TTTT", "GATTACA"):
        assert str(pg.Kmer(seq)) == seq
        assert len(pg.Kmer(seq)) == len(seq)
        assert pg.Kmer(seq).k == len(seq)


def test_lowercase_is_accepted_and_normalized():
    pg = _pg()
    assert str(pg.Kmer("acgt")) == "ACGT"
    assert pg.Kmer("acgt") == pg.Kmer("ACGT")


def test_default_is_empty():
    pg = _pg()
    assert pg.Kmer().k == 0
    assert len(pg.Kmer()) == 0
    assert str(pg.Kmer()) == ""


def test_encoding_constructor_roundtrip():
    pg = _pg()
    km = pg.Kmer("ACGT")
    assert pg.Kmer(km.encoding, km.k) == km
    assert str(pg.Kmer(km.encoding, km.k)) == "ACGT"


def test_overlap_requires_same_sequence_and_length():
    pg = _pg()
    assert pg.Kmer.overlaps(pg.Kmer("ACGT"), pg.Kmer("ACGT")) is True
    assert pg.Kmer.overlaps(pg.Kmer("ACGT"), pg.Kmer("ACGA")) is False
    # Same prefix but different length -> not equal.
    assert pg.Kmer.overlaps(pg.Kmer("ACG"), pg.Kmer("ACGT")) is False


def test_equality_and_repr():
    pg = _pg()
    assert pg.Kmer("ACGT") == pg.Kmer("ACGT")
    assert pg.Kmer("ACGT") != pg.Kmer("TGCA")
    assert repr(pg.Kmer("ACGT")) == "Kmer('ACGT')"


def test_ordering_length_then_encoding():
    pg = _pg()
    # Shorter k sorts first; within equal length, A < C < G < T.
    assert pg.Kmer("AC") < pg.Kmer("ACG")          # length first
    assert pg.Kmer("AAAA") < pg.Kmer("AAAC")        # encoding within equal length


def test_invalid_base_raises():
    pg = _pg()
    with pytest.raises(ValueError):
        pg.Kmer("ACGN")
    with pytest.raises(ValueError):
        pg.Kmer("xyz")


def test_length_over_32_raises():
    pg = _pg()
    with pytest.raises(ValueError):
        pg.Kmer("A" * 33)


def test_is_valid_and_max_k():
    pg = _pg()
    assert pg.Kmer.is_valid("ACGTacgt") is True
    assert pg.Kmer.is_valid("ACGN") is False
    assert pg.Kmer.max_k == 32