"""
Tests for GffReader — the single-pass iterator over GFF3/GTF files.

Ports the applicable cases from genogrove gfffile-test.cpp (readGFF3Format,
readGTFFormat, format detection, GTF helpers, skip/throw on invalid,
file-not-found, validate_gtf, gzip) to the Python iterator surface.
"""

import gzip

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _write(path, lines):
    """Write GFF lines (each a list of 9 columns) and return the path as str."""
    path.write_text("".join("\t".join(map(str, c)) + "\n" for c in lines))
    return str(path)


GFF3 = [
    ["chr1", "HAVANA", "gene", 1000, 2000, ".", "+", ".",
     "ID=gene1;Name=TEST1;biotype=protein_coding"],
    ["chr1", "HAVANA", "exon", 1000, 1500, ".", "+", ".",
     "ID=exon1;Parent=gene1"],
    ["chr2", "HAVANA", "gene", 3000, 4000, 100, "-", ".", "ID=gene2"],
    ["chrX", "HAVANA", "CDS", 5000, 5100, ".", "+", 0, "ID=cds1;Parent=gene2"],
]

GTF = [
    ["chr1", "HAVANA", "gene", 1000, 2000, ".", "+", ".",
     'gene_id "ENSG00000001"; gene_name "TEST1"; gene_biotype "protein_coding";'],
    ["chr1", "HAVANA", "exon", 1000, 1500, ".", "+", ".",
     'gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "1";'],
    ["chr2", "HAVANA", "gene", 3000, 4000, 100, "-", ".",
     'gene_id "ENSG00000002";'],
    ["chrX", "HAVANA", "CDS", 5000, 5100, ".", "+", 0,
     'gene_id "ENSG00000002"; transcript_id "ENST00000002";'],
]


def test_read_gff3(tmp_path):
    """A GFF3 file yields entries with key=value attributes. Mirrors readGFF3Format."""
    pg = _pg()
    entries = list(pg.GffReader(_write(tmp_path / "a.gff3", GFF3)))
    assert len(entries) == 4

    g = entries[0]
    assert g.seqid == "chr1"
    assert g.source == "HAVANA"
    assert g.type == "gene"
    assert g.start == 1000          # native GFF 1-based inclusive
    assert g.end == 2000
    assert g.score is None
    assert g.strand == "+"
    assert g.phase is None
    assert g.format == pg.GffFormat.GFF3
    assert g.attributes["ID"] == "gene1"
    assert g.attributes["Name"] == "TEST1"
    assert g.attributes["biotype"] == "protein_coding"

    assert entries[1].type == "exon"
    assert entries[1].attributes["Parent"] == "gene1"
    assert entries[2].score == 100
    assert entries[2].strand == "-"
    assert entries[3].type == "CDS"
    assert entries[3].phase == 0


def test_read_gtf(tmp_path):
    """A GTF file is detected as GTF; helpers read the GTF attributes. Mirrors readGTFFormat."""
    pg = _pg()
    entries = list(pg.GffReader(_write(tmp_path / "a.gtf", GTF)))
    assert len(entries) == 4

    g = entries[0]
    assert g.format == pg.GffFormat.GTF
    assert g.is_gtf()
    assert not g.is_gff3()
    assert g.get_gene_id() == "ENSG00000001"
    assert g.get_gene_name() == "TEST1"
    assert g.get_gene_biotype() == "protein_coding"

    exon = entries[1]
    assert exon.get_transcript_id() == "ENST00000001"
    assert exon.get_exon_number() == 1
    assert entries[2].score == 100
    assert entries[3].phase == 0


def test_empty_file(tmp_path):
    """An empty GFF file yields no entries."""
    pg = _pg()
    (tmp_path / "e.gff3").write_text("")
    assert list(pg.GffReader(str(tmp_path / "e.gff3"))) == []


def test_comment_lines_skipped(tmp_path):
    """GFF3 '##' header / '#' comment lines are not records."""
    pg = _pg()
    path = tmp_path / "h.gff3"
    path.write_text("##gff-version 3\n# a comment\n"
                    "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g1\n")
    entries = list(pg.GffReader(str(path)))
    assert len(entries) == 1
    assert entries[0].attributes["ID"] == "g1"


def test_invalid_line_raises_by_default(tmp_path):
    """With skip_invalid_lines=False (default), a malformed line raises."""
    pg = _pg()
    bad = [
        ["chr1", "src", "gene", 1000, 2000, ".", "+", ".", "ID=g1"],
        ["chr2", "src", "gene", "NOTANUMBER", 4000, ".", "-", ".", "ID=g2"],
    ]
    with pytest.raises(RuntimeError):
        list(pg.GffReader(_write(tmp_path / "bad.gff3", bad)))


def test_skip_invalid_lines(tmp_path):
    """With skip_invalid_lines=True, malformed lines are skipped."""
    pg = _pg()
    bad = [
        ["chr1", "src", "gene", 1000, 2000, ".", "+", ".", "ID=g1"],
        ["chr2", "src", "gene", "NOTANUMBER", 4000, ".", "-", ".", "ID=g2"],
        ["chr3", "src", "gene", 5000, 6000, ".", "+", ".", "ID=g3"],
    ]
    entries = list(pg.GffReader(_write(tmp_path / "bad.gff3", bad),
                                skip_invalid_lines=True))
    assert [e.seqid for e in entries] == ["chr1", "chr3"]


def test_invalid_first_line_raises_even_when_skipping(tmp_path):
    """The first data record is validated at construction, so a malformed FIRST
    line raises immediately — even with skip_invalid_lines=True."""
    pg = _pg()
    bad = [["chr1", "src", "gene", "NOTANUMBER", 2000, ".", "+", ".", "ID=g1"]]
    path = _write(tmp_path / "badfirst.gff3", bad)
    with pytest.raises(RuntimeError):
        pg.GffReader(path)
    with pytest.raises(RuntimeError):
        pg.GffReader(path, skip_invalid_lines=True)


def test_file_not_found(tmp_path):
    """Opening a missing file raises."""
    pg = _pg()
    with pytest.raises((RuntimeError, IOError, OSError)):
        pg.GffReader(str(tmp_path / "nope.gff3"))


def test_validate_gtf_rejects_missing_gene_id(tmp_path):
    """validate_gtf=True rejects a GTF record lacking the mandatory gene_id."""
    pg = _pg()
    lines = [["chr1", "src", "gene", 1000, 2000, ".", "+", ".",
              'transcript_id "ENST1";']]   # GTF-style but no gene_id
    path = _write(tmp_path / "v.gtf", lines)
    with pytest.raises(RuntimeError):
        list(pg.GffReader(path, validate_gtf=True))
    # off by default: the same file iterates without raising
    assert len(list(pg.GffReader(path))) == 1


def test_start_one_boundary(tmp_path):
    """A feature starting at GFF position 1 reads back as start == 1."""
    pg = _pg()
    path = _write(tmp_path / "b.gff3",
                  [["chr1", "src", "gene", 1, 100, ".", "+", ".", "ID=g1"]])
    e = list(pg.GffReader(path))[0]
    assert e.start == 1
    assert e.end == 100


def test_reads_gzip(tmp_path):
    """A gzip-compressed GFF file is decompressed transparently."""
    pg = _pg()
    path = tmp_path / "a.gff3.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("chr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tID=g1\n")
    entries = list(pg.GffReader(str(path)))
    assert len(entries) == 1
    assert entries[0].attributes["ID"] == "g1"


def test_build_grove_from_reader(tmp_path):
    """The common workflow: read a GFF file into a GffGrove.

    GFF coordinates are 1-based inclusive [start, end]; the grove key is the
    0-based closed interval [start - 1, end - 1].
    """
    pg = _pg()
    g = pg.GffGrove(64)
    for e in pg.GffReader(_write(tmp_path / "a.gff3", GFF3)):
        g.insert(e.seqid, pg.GenomicCoordinate(".", e.start - 1, e.end - 1), e)

    assert g.size() == 4
    hits = list(g.intersect(pg.GenomicCoordinate(".", 1200, 1200), "chr1"))   # inside gene + exon
    types = sorted(h.data.type for h in hits)
    assert types == ["exon", "gene"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])