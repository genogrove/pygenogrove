"""
Serialization round-trip tests for GffGrove (grove<interval, gff_entry>).

Ports genogrove gfffile-test.cpp::gffEntrySerialization (minimalRoundTrip,
gff3FullRoundTrip, gtfRoundTrip, emptyAttributesMap,
preservesEmbeddedNulInAttributes, groveRoundTrip) to the Python grove surface:
gff_entry serialization is only reachable through the grove, so each case
inserts the entry, serializes the grove, deserializes, and checks the payload.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _roundtrip(pg, tmp_path, grove):
    path = str(tmp_path / "gff.gg")
    grove.serialize(path)
    return pg.GffGrove.deserialize(path)


def _only(grove, lo, hi, index="chr1"):
    import pygenogrove as pg
    hits = list(grove.intersect(pg.Interval(lo, hi), index))
    assert len(hits) == 1
    return hits[0].data


def test_minimal_round_trip(tmp_path):
    """No optionals; empty attributes; UNKNOWN format — all round-trip absent."""
    pg = _pg()
    g = pg.GffGrove(4)
    g.insert("chr1", pg.Interval(100, 200), pg.GffEntry("chr1", 100, 200, "gene"))

    out = _only(_roundtrip(pg, tmp_path, g), 150, 150)
    assert out.seqid == "chr1"
    assert out.start == 100
    assert out.end == 200
    assert out.type == "gene"
    assert out.source == ""
    assert out.score is None
    assert out.strand is None
    assert out.phase is None
    assert dict(out.attributes) == {}
    assert out.format == pg.GffFormat.UNKNOWN


def test_gff3_full_round_trip(tmp_path):
    """All optional columns + GFF3 attributes survive."""
    pg = _pg()
    e = pg.GffEntry("chrX", 1000, 2000, "exon")
    e.source = "ENSEMBL"
    e.score = 42.5
    e.strand = "-"
    e.phase = 1
    e.format = pg.GffFormat.GFF3
    e.attributes = {
        "ID": "exon:ENSE00001",
        "Parent": "transcript:ENST00001",
        "gene_id": "ENSG00001",
    }
    g = pg.GffGrove(4)
    g.insert("chrX", pg.Interval(1000, 2000), e)

    out = _only(_roundtrip(pg, tmp_path, g), 1500, 1500, "chrX")
    assert out.source == "ENSEMBL"
    assert out.score == 42.5
    assert out.strand == "-"
    assert out.phase == 1
    assert out.format == pg.GffFormat.GFF3
    assert dict(out.attributes) == {
        "ID": "exon:ENSE00001",
        "Parent": "transcript:ENST00001",
        "gene_id": "ENSG00001",
    }


def test_gtf_round_trip(tmp_path):
    """GTF format flag + attributes survive; is_gtf/is_gff3 reflect it."""
    pg = _pg()
    e = pg.GffEntry("chr2", 500, 600, "CDS")
    e.format = pg.GffFormat.GTF
    e.phase = 0
    e.attributes = {"gene_id": "G1", "transcript_id": "T1"}
    g = pg.GffGrove(4)
    g.insert("chr2", pg.Interval(500, 600), e)

    out = _only(_roundtrip(pg, tmp_path, g), 550, 550, "chr2")
    assert out.format == pg.GffFormat.GTF
    assert out.is_gtf()
    assert not out.is_gff3()
    assert out.attributes["gene_id"] == "G1"
    assert out.attributes["transcript_id"] == "T1"


def test_empty_attributes_map(tmp_path):
    """An explicit zero-entry attributes map round-trips empty."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1, 10, "region")
    e.format = pg.GffFormat.GFF3
    g = pg.GffGrove(3)
    g.insert("chr1", pg.Interval(1, 10), e)

    out = _only(_roundtrip(pg, tmp_path, g), 5, 5)
    assert dict(out.attributes) == {}
    assert out.format == pg.GffFormat.GFF3


def test_preserves_embedded_nul_in_attributes(tmp_path):
    """Attribute keys/values go through length-prefixed strings, so embedded
    NULs survive. Mirrors gffEntrySerialization::preservesEmbeddedNulInAttributes."""
    pg = _pg()
    e = pg.GffEntry("chr1", 5, 10, "feature")
    e.attributes = {"k\x00ey": "va\x00lue"}
    g = pg.GffGrove(3)
    g.insert("chr1", pg.Interval(5, 10), e)

    out = _only(_roundtrip(pg, tmp_path, g), 7, 7)
    attrs = dict(out.attributes)
    assert attrs == {"k\x00ey": "va\x00lue"}
    assert len(list(attrs.values())[0]) == 6


def test_grove_round_trip(tmp_path):
    """End-to-end multi-entry grove round-trip. Mirrors gffEntrySerialization::groveRoundTrip."""
    pg = _pg()
    g = pg.GffGrove(4)

    e1 = pg.GffEntry("chr1", 100, 200, "gene")
    e1.format = pg.GffFormat.GFF3
    e1.attributes = {"ID": "gene1"}
    e1.strand = "+"

    e2 = pg.GffEntry("chr1", 300, 400, "exon")
    e2.format = pg.GffFormat.GFF3
    e2.attributes = {"Parent": "gene1"}

    g.insert("chr1", pg.Interval(100, 200), e1)
    g.insert("chr1", pg.Interval(300, 400), e2)

    restored = _roundtrip(pg, tmp_path, g)

    d1 = _only(restored, 150, 150)
    assert d1.type == "gene"
    assert d1.attributes["ID"] == "gene1"
    assert d1.strand == "+"

    d2 = _only(restored, 350, 350)
    assert d2.type == "exon"
    assert d2.attributes["Parent"] == "gene1"


def test_roundtrip_preserves_edges_and_external_data(tmp_path):
    """Edges + external GffEntry keys survive the round-trip."""
    pg = _pg()
    g = pg.GffGrove(3)
    exon = g.insert("chr1", pg.Interval(1000, 1200),
                    pg.GffEntry("chr1", 1000, 1200, "exon"))
    enh = g.add_external_key(pg.Interval(5000, 5500),
                             pg.GffEntry("chr1", 5000, 5500, "enhancer"))
    g.add_edge(exon, enh)

    restored = _roundtrip(pg, tmp_path, g)
    assert restored.edge_count() == 1
    src = list(restored.intersect(pg.Interval(1000, 1200), "chr1"))[0]
    neighbors = restored.get_neighbors(src)
    assert len(neighbors) == 1
    assert neighbors[0].data.type == "enhancer"
    assert len(restored.intersect(pg.Interval(5000, 5500))) == 0


def test_roundtrip_empty_gff_grove(tmp_path):
    pg = _pg()
    g = pg.GffGrove(3)
    restored = _roundtrip(pg, tmp_path, g)
    assert restored.size() == 0
    assert restored.edge_count() == 0
    assert restored.get_order() == 3


def test_deserialize_corrupt_file_raises(tmp_path):
    """A file that is not a valid compressed grove fails loudly."""
    pg = _pg()
    bad = tmp_path / "garbage.gg"
    bad.write_bytes(b"this is definitely not a zlib-compressed grove")
    with pytest.raises((RuntimeError, ValueError, IOError, OSError)):
        pg.GffGrove.deserialize(str(bad))


def test_deserialize_missing_file_raises(tmp_path):
    pg = _pg()
    with pytest.raises((RuntimeError, IOError, OSError)):
        pg.GffGrove.deserialize(str(tmp_path / "does_not_exist.gg"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])