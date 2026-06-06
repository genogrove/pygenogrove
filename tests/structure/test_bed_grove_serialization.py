"""
Serialization round-trip tests for BedGrove (grove<interval, bed_entry>).

The critical check: the associated BedEntry data survives a serialize /
deserialize cycle through the zlib-compressed .gg format (exercises the
out-of-line bed_entry::serialize / deserialize linked from the genogrove lib).
Mirrors tests/structure/test_serialization.py for the dataless grove.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_roundtrip_preserves_bed_data(tmp_path):
    pg = _pg()

    g = pg.BedGrove(4)
    e1 = pg.BedEntry("chr1", 100, 201)
    e1.name = "gene1"
    e1.score = 500
    e1.strand = "+"
    g.insert("chr1", pg.Interval(100, 200), e1)
    g.insert("chr1", pg.Interval(300, 400), pg.BedEntry("chr1", 300, 401))

    path = str(tmp_path / "bed.gg")
    g.serialize(path)

    loaded = pg.BedGrove.deserialize(path)
    assert loaded.size() == 2
    assert loaded.get_order() == 4

    hits = list(loaded.intersect(pg.Interval(100, 200), "chr1"))
    assert len(hits) == 1
    data = hits[0].data
    assert data.chrom == "chr1"
    assert data.start == 100
    assert data.end == 201
    assert data.name == "gene1"
    assert data.score == 500
    assert data.strand == "+"


def test_roundtrip_preserves_block_info(tmp_path):
    pg = _pg()

    g = pg.BedGrove(3)
    entry = pg.BedEntry("chr1", 1000, 1500)
    entry.blocks = pg.BlockInfo(2, [100, 200], [0, 300])
    g.insert("chr1", pg.Interval(1000, 1499), entry)

    path = str(tmp_path / "blocks.gg")
    g.serialize(path)

    loaded = pg.BedGrove.deserialize(path)
    hit = list(loaded.intersect(pg.Interval(1000, 1499), "chr1"))[0]
    assert hit.data.blocks is not None
    assert hit.data.blocks.count == 2
    assert list(hit.data.blocks.sizes) == [100, 200]
    assert list(hit.data.blocks.starts) == [0, 300]


def test_roundtrip_preserves_edges_and_external_data(tmp_path):
    pg = _pg()

    g = pg.BedGrove(3)
    exon = g.insert("chr1", pg.Interval(1000, 1200), pg.BedEntry("chr1", 1000, 1200))
    enhancer = g.add_external_key(pg.Interval(5000, 5500),
                                  pg.BedEntry("chr1", 5000, 5500))
    g.add_edge(exon, enhancer)

    path = str(tmp_path / "bed_edges.gg")
    g.serialize(path)

    loaded = pg.BedGrove.deserialize(path)
    assert loaded.edge_count() == 1

    src = list(loaded.intersect(pg.Interval(1000, 1200), "chr1"))[0]
    neighbors = loaded.get_neighbors(src)
    assert len(neighbors) == 1
    assert neighbors[0].value.start == 5000
    assert neighbors[0].data.start == 5000

    # External key still excluded from spatial queries after reload.
    assert len(loaded.intersect(pg.Interval(5000, 5500))) == 0
    assert loaded.size() == 1


def test_roundtrip_empty_bed_grove(tmp_path):
    pg = _pg()
    g = pg.BedGrove(3)
    path = str(tmp_path / "empty_bed.gg")
    g.serialize(path)

    loaded = pg.BedGrove.deserialize(path)
    assert loaded.size() == 0
    assert loaded.edge_count() == 0
    assert loaded.get_order() == 3


if __name__ == "__main__":
    pytest.main([__file__, "-v"])