"""
Tests for key removal + storage compaction and the vertex/storage counts:
remove_key(), compact(), vertex_count(), external_vertex_count(),
key_storage_size().

Ports the behaviourally-observable cases from genogrove
tests/structure/grove_remove_test.cpp (GroveRemoveTest / GroveCompactTest). The
pure B+ tree rebalancing internals (borrow/merge/root-collapse) are exercised
through the binding by test_stress_remove_with_rebalancing rather than ported
one-to-one — those are genogrove's own unit tests. compact() invalidates
previously-returned indexed Key handles, so every post-compact assertion
re-discovers keys via a fresh query.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


# ---------------------------------------------------------------- remove_key

def test_remove_single_key_from_root_leaf():
    pg = _pg()
    g = pg.Grove(4)
    k = g.insert("chr1", pg.Interval(100, 200))
    assert g.remove_key("chr1", k) is True
    assert g.indexed_vertex_count() == 0
    assert len(g.intersect(pg.Interval(100, 200), "chr1")) == 0


def test_remove_same_key_twice():
    pg = _pg()
    g = pg.Grove(4)
    k = g.insert("chr1", pg.Interval(100, 200))
    assert g.remove_key("chr1", k) is True
    assert g.remove_key("chr1", k) is False


def test_remove_null_key_returns_false():
    """remove_key(None) is a no-op returning False (mirrors RemoveNullKey)."""
    pg = _pg()
    g = pg.Grove(4)
    assert g.remove_key("chr1", None) is False


def test_remove_from_nonexistent_index_returns_false():
    pg = _pg()
    g = pg.Grove(4)
    k = g.insert("chr1", pg.Interval(100, 200))
    assert g.remove_key("chrX", k) is False
    assert g.indexed_vertex_count() == 1


def test_remove_all_keys_empties_index():
    pg = _pg()
    g = pg.Grove(4)
    keys = [g.insert("chr1", pg.Interval(i * 10, i * 10 + 5)) for i in range(3)]
    for k in keys:
        assert g.remove_key("chr1", k) is True
    assert g.indexed_vertex_count() == 0
    assert len(g.intersect(pg.Interval(0, 1000), "chr1")) == 0


def test_remove_multi_index_independence():
    """Removing from one index leaves other indices untouched."""
    pg = _pg()
    g = pg.Grove(4)
    chr1 = [g.insert("chr1", pg.Interval(i * 10, i * 10 + 5)) for i in range(5)]
    [g.insert("chr2", pg.Interval(i * 10, i * 10 + 5)) for i in range(5)]

    assert g.remove_key("chr1", chr1[2]) is True
    assert len(g.intersect(pg.Interval(0, 1000), "chr1")) == 4
    assert len(g.intersect(pg.Interval(0, 1000), "chr2")) == 5


def test_remove_key_removes_incoming_and_outgoing_edges():
    pg = _pg()
    g = pg.Grove(4)
    a = g.insert("chr1", pg.Interval(100, 200))
    b = g.insert("chr1", pg.Interval(300, 400))
    c = g.insert("chr1", pg.Interval(500, 600))
    g.add_edge(a, b)   # incoming to b
    g.add_edge(b, c)   # outgoing from b
    assert g.edge_count() == 2

    g.remove_key("chr1", b)
    assert g.edge_count() == 0  # both edges touching b are gone
    assert g.has_edge(a, b) is False
    assert g.has_edge(b, c) is False


def test_stress_remove_with_rebalancing():
    """Drive splits + borrow/merge/root-collapse through the binding (order 4)."""
    pg = _pg()
    g = pg.Grove(4)
    keys = [g.insert("chr1", pg.Interval(i * 10, i * 10 + 5)) for i in range(40)]
    assert g.indexed_vertex_count() == 40

    # Remove the first 25 (forces underflow / merges / root collapse).
    for k in keys[:25]:
        assert g.remove_key("chr1", k) is True
    assert g.indexed_vertex_count() == 15

    # The remaining 15 are still all discoverable in order.
    res = g.intersect(pg.Interval(0, 100000), "chr1")
    assert len(res) == 15
    assert sorted(k.value.start for k in res) == [i * 10 for i in range(25, 40)]

    # Remove the rest -> empty.
    for k in keys[25:]:
        assert g.remove_key("chr1", k) is True
    assert g.indexed_vertex_count() == 0


# ------------------------------------------------------------------- counts

def test_vertex_counts():
    pg = _pg()
    g = pg.Grove(4)
    g.insert("chr1", pg.Interval(100, 200))
    g.insert("chr1", pg.Interval(300, 400))
    g.add_external_key(pg.Interval(500, 600))

    assert g.indexed_vertex_count() == 2
    assert g.external_vertex_count() == 1
    assert g.vertex_count() == 3  # indexed + external


# ------------------------------------------------------------------ compact

def test_compact_noop_on_empty_grove():
    pg = _pg()
    g = pg.Grove(4)
    assert g.key_storage_size() == 0
    g.compact()
    assert g.key_storage_size() == 0
    assert g.indexed_vertex_count() == 0


def test_compact_preserves_keys_without_removal():
    """compact() with no prior removal keeps storage and all keys intact."""
    pg = _pg()
    g = pg.Grove(4)
    for i in range(20):
        g.insert("chr1", pg.Interval(i * 10, i * 10 + 5))
    before = g.key_storage_size()

    g.compact()
    assert g.key_storage_size() == before
    assert g.indexed_vertex_count() == 20
    for i in range(20):
        assert len(g.intersect(pg.Interval(i * 10, i * 10 + 5), "chr1")) == 1


def test_compact_reclaims_storage_after_removals():
    pg = _pg()
    g = pg.Grove(4)
    keys = [g.insert("chr1", pg.Interval(i * 10, i * 10 + 5)) for i in range(20)]
    for k in keys[:10]:
        assert g.remove_key("chr1", k) is True

    before = g.key_storage_size()
    assert g.indexed_vertex_count() == 10
    assert g.key_storage_size() == before  # dead slots not yet reclaimed

    g.compact()
    assert g.key_storage_size() < before
    assert g.key_storage_size() >= g.indexed_vertex_count()
    assert g.indexed_vertex_count() == 10


def test_compact_queries_work_after_and_requery():
    """After compact(), old handles are invalid; re-querying still finds keys."""
    pg = _pg()
    g = pg.Grove(4)
    keys = [g.insert("chr1", pg.Interval(i * 10, i * 10 + 5)) for i in range(20)]
    for k in keys[:10]:
        g.remove_key("chr1", k)

    g.compact()  # invalidates every handle in `keys` — do not touch them

    res = g.intersect(pg.Interval(0, 100000), "chr1")
    assert len(res) == 10
    assert sorted(k.value.start for k in res) == [i * 10 for i in range(10, 20)]


def test_compact_graph_edges_survive():
    """Edges between surviving indexed keys are remapped, not dropped."""
    pg = _pg()
    g = pg.Grove(4)
    keys = [g.insert("chr1", pg.Interval(i * 10, i * 10 + 5)) for i in range(10)]
    for i in range(0, 10, 2):
        g.add_edge(keys[i], keys[i + 1])
    pre = g.edge_count()

    assert g.remove_key("chr1", keys[0]) is True  # drops 1 edge
    assert g.edge_count() == pre - 1

    g.compact()
    assert g.edge_count() == pre - 1

    # Re-discover a surviving source/target pair and confirm the edge remapped.
    src = list(g.intersect(pg.Interval(20, 25), "chr1"))[0]
    tgt = list(g.intersect(pg.Interval(30, 35), "chr1"))[0]
    assert g.has_edge(src, tgt) is True


def test_compact_external_keys_and_cross_edges_preserved():
    """External keys keep their pointer/data and cross-edges across compact()."""
    pg = _pg()
    g = pg.Grove(4)
    indexed = g.insert("chr1", pg.Interval(0, 5))
    external = g.add_external_key(pg.Interval(1000, 1100))
    g.add_edge(indexed, external)   # indexed -> external
    g.add_edge(external, indexed)   # external -> indexed

    victim = g.insert("chr1", pg.Interval(10, 15))
    assert g.remove_key("chr1", victim) is True

    g.compact()

    # External key handle is still valid after compact().
    assert external.value == pg.Interval(1000, 1100)
    assert g.external_vertex_count() == 1

    # The indexed endpoint moved; re-discover it and check both edge directions.
    new_indexed = list(g.intersect(pg.Interval(0, 5), "chr1"))[0]
    assert g.has_edge(new_indexed, external) is True
    assert g.has_edge(external, new_indexed) is True


def test_compact_roundtrip_serialization(tmp_path):
    """A grove serialized after compact() round-trips correctly."""
    pg = _pg()
    g = pg.Grove(4)
    keys = [g.insert("chr1", pg.Interval(i * 10, i * 10 + 5)) for i in range(20)]
    for k in keys[:10]:
        g.remove_key("chr1", k)
    g.compact()

    path = str(tmp_path / "compacted.gg")
    g.serialize(path)
    loaded = pg.Grove.deserialize(path)

    assert loaded.indexed_vertex_count() == 10
    res = loaded.intersect(pg.Interval(0, 100000), "chr1")
    assert sorted(k.value.start for k in res) == [i * 10 for i in range(10, 20)]


# ----------------------------------------------------- data-carrying grove

def test_remove_key_on_data_grove():
    pg = _pg()
    g = pg.BedGrove(4)
    a = g.insert("chr1", pg.Interval(100, 200), pg.BedEntry("chr1", 100, 201))
    g.insert("chr1", pg.Interval(300, 400), pg.BedEntry("chr1", 300, 401))
    assert g.size() == 2

    assert g.remove_key("chr1", a) is True
    assert g.size() == 1
    assert len(g.intersect(pg.Interval(100, 200), "chr1")) == 0