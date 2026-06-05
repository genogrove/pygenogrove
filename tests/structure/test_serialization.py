"""
Tests for Grove serialization / deserialization (zlib-compressed .gg binary).

Mirrors genogrove/tests/structure/serialization_test.cpp (``SerializationTest``)
plus the round-trip cases from ExternalKeyTest's serialization coverage, over
the surface the Python bindings expose.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_roundtrip_intervals(tmp_path):
    pg = _pg()

    g = pg.Grove(4)
    g.insert("chr1", pg.Interval(100, 200))
    g.insert("chr1", pg.Interval(300, 400))
    g.insert("chr2", pg.Interval(150, 250))

    path = str(tmp_path / "groove.gg")
    g.serialize(path)

    loaded = pg.Grove.deserialize(path)
    assert loaded.size() == 3
    assert loaded.get_order() == 4

    hits = loaded.intersect(pg.Interval(150, 350), "chr1")
    assert len(hits) == 2


def test_roundtrip_preserves_edges(tmp_path):
    pg = _pg()

    g = pg.Grove(3)
    a = g.insert("chr1", pg.Interval(100, 200))
    b = g.insert("chr1", pg.Interval(300, 400))
    g.add_edge(a, b)

    path = str(tmp_path / "with_edges.gg")
    g.serialize(path)

    loaded = pg.Grove.deserialize(path)
    assert loaded.edge_count() == 1

    # recover the source key in the loaded grove and traverse the edge
    src = list(loaded.intersect(pg.Interval(100, 200), "chr1"))[0]
    neighbors = loaded.get_neighbors(src)
    assert len(neighbors) == 1
    assert neighbors[0].value.start == 300


def test_roundtrip_empty_grove(tmp_path):
    pg = _pg()

    g = pg.Grove(3)
    path = str(tmp_path / "empty.gg")
    g.serialize(path)

    loaded = pg.Grove.deserialize(path)
    assert loaded.size() == 0
    assert loaded.edge_count() == 0
    assert loaded.get_order() == 3


def test_roundtrip_preserves_external_keys(tmp_path):
    pg = _pg()

    g = pg.Grove(5)
    exon = g.insert("chr1", pg.Interval(1000, 1200))
    enhancer = g.add_external_key(pg.Interval(5000, 5500))
    g.add_edge(exon, enhancer)

    path = str(tmp_path / "external.gg")
    g.serialize(path)

    loaded = pg.Grove.deserialize(path)
    assert loaded.edge_count() == 1

    # Recover the indexed source and traverse to the external target.
    src = list(loaded.intersect(pg.Interval(1000, 1200), "chr1"))[0]
    neighbors = loaded.get_neighbors(src)
    assert len(neighbors) == 1
    assert neighbors[0].value.start == 5000

    # The external key is still excluded from spatial queries after reload.
    assert len(loaded.intersect(pg.Interval(5000, 5500))) == 0
    assert loaded.size() == 1  # external key is not an indexed vertex


def test_roundtrip_edge_chain_across_splits(tmp_path):
    """Order 3 forces splits with separator keys; the flat edge index must
    still map correctly after a serialize/deserialize round-trip."""
    pg = _pg()

    g = pg.Grove(3)
    e1 = g.insert("chr1", pg.Interval(1000, 1200))
    e2 = g.insert("chr1", pg.Interval(2000, 2200))
    e3 = g.insert("chr1", pg.Interval(3000, 3200))
    ext = g.add_external_key(pg.Interval(5000, 5500))
    g.add_edge(e1, e2)
    g.add_edge(e2, e3)
    g.add_edge(e3, ext)

    path = str(tmp_path / "chain.gg")
    g.serialize(path)

    loaded = pg.Grove.deserialize(path)
    assert loaded.edge_count() == 3

    keys = list(loaded.intersect(pg.Interval(1000, 3200), "chr1"))
    assert len(keys) == 3
    keys.sort(key=lambda k: k.value.start)
    # exon1 -> exon2 -> exon3 -> external enhancer
    assert loaded.get_neighbors(keys[0])[0].value.start == 2000
    assert loaded.get_neighbors(keys[1])[0].value.start == 3000
    assert loaded.get_neighbors(keys[2])[0].value.start == 5000


def test_deserialize_corrupt_file_raises(tmp_path):
    """A file that is not a valid compressed grove fails loudly."""
    pg = _pg()
    bad = tmp_path / "garbage.gg"
    bad.write_bytes(b"this is definitely not a zlib-compressed grove")
    with pytest.raises((RuntimeError, ValueError, IOError, OSError)):
        pg.Grove.deserialize(str(bad))


def test_deserialize_missing_file_raises(tmp_path):
    pg = _pg()
    missing = str(tmp_path / "does_not_exist.gg")
    with pytest.raises((RuntimeError, IOError, OSError)):
        pg.Grove.deserialize(missing)


def test_serialize_to_unwritable_path_raises(tmp_path):
    pg = _pg()
    g = pg.Grove(3)
    g.insert("chr1", pg.Interval(100, 200))
    # a directory component that does not exist -> open fails
    bad = str(tmp_path / "no_such_dir" / "out.gg")
    with pytest.raises((RuntimeError, IOError, OSError)):
        g.serialize(bad)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])