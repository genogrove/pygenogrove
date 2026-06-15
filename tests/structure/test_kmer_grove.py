"""
Behaviour of KmerGrove (grove<kmer, json_value, json_value>): a k-mer dictionary
for membership lookups. Overlap is exact equality (same bases AND length), so
intersect() answers "is this k-mer present". Exercised dataless and with a JSON
payload.

Mirrors genogrove/tests/structure/kmer_grove_test.cpp over the bound surface,
plus one split-stress insert and a serialization round-trip.
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_creation_and_insert():
    pg = _pg()
    g = pg.KmerGrove(8)
    assert g.size() == 0

    key = g.insert("seqs", pg.Kmer("ACGT"))
    assert g.size() == 1
    assert str(key.value) == "ACGT"
    assert key.data is None


def test_membership_lookup():
    pg = _pg()
    g = pg.KmerGrove(8)
    for seq in ("ACGT", "AAAA", "TTTT", "GATC"):
        g.insert("seqs", pg.Kmer(seq))

    assert [str(k.value) for k in g.intersect(pg.Kmer("ACGT"), "seqs")] == ["ACGT"]
    assert len(g.intersect(pg.Kmer("ACGA"), "seqs")) == 0   # absent
    assert len(g.intersect(pg.Kmer("ACG"), "seqs")) == 0    # different length


def test_json_payload_counts():
    pg = _pg()
    g = pg.KmerGrove(8)
    g.insert("seqs", pg.Kmer("ACGT"), {"count": 3})
    hit = list(g.intersect(pg.Kmer("ACGT"), "seqs"))[0]
    assert hit.data == {"count": 3}


def test_many_inserts_force_splits():
    pg = _pg()
    g = pg.KmerGrove(3)
    bases = "ACGT"
    seqs = [a + b + c for a in bases for b in bases for c in bases]  # 64 distinct 3-mers
    for s in seqs:
        g.insert("seqs", pg.Kmer(s))
    assert g.size() == 64
    for s in ("AAA", "ACG", "TTT", "GAT"):
        assert [str(k.value) for k in g.intersect(pg.Kmer(s), "seqs")] == [s]


def test_serialization_roundtrip(tmp_path):
    pg = _pg()
    g = pg.KmerGrove(4)
    for s in ("ACGT", "GGGG", "TACG"):
        g.insert("seqs", pg.Kmer(s), {"seq": s})
    path = str(tmp_path / "kmers.gg")
    g.serialize(path)

    loaded = pg.KmerGrove.deserialize(path)
    assert loaded.size() == 3
    hit = list(loaded.intersect(pg.Kmer("TACG"), "seqs"))[0]
    assert str(hit.value) == "TACG"
    assert hit.data == {"seq": "TACG"}


def test_keys_keep_grove_alive():
    """A KmerKey keeps its grove alive (use-after-free guard)."""
    pg = _pg()
    g = pg.KmerGrove(3)
    key = g.insert("seqs", pg.Kmer("ACGT"))
    del g
    gc.collect()
    assert str(key.value) == "ACGT"