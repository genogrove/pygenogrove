"""
Tests for Grove serialization / deserialization (zlib-compressed .gg binary).
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


class TestSerialize:
    def test_roundtrip_intervals(self, tmp_path):
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

    def test_roundtrip_preserves_edges(self, tmp_path):
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

    def test_deserialize_missing_file_raises(self, tmp_path):
        pg = _pg()
        missing = str(tmp_path / "does_not_exist.gg")
        with pytest.raises((RuntimeError, IOError, OSError)):
            pg.Grove.deserialize(missing)

    def test_serialize_to_unwritable_path_raises(self, tmp_path):
        pg = _pg()
        g = pg.Grove(3)
        g.insert("chr1", pg.Interval(100, 200))
        # a directory component that does not exist -> open fails
        bad = str(tmp_path / "no_such_dir" / "out.gg")
        with pytest.raises((RuntimeError, IOError, OSError)):
            g.serialize(bad)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])