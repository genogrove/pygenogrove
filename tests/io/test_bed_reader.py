"""
Tests for BedReader — the single-pass iterator over BED files.

Ports the applicable cases from genogrove bedfile-test.cpp (readBED3/BED6,
iteration, skip/throw on invalid, file-not-found, line counter, gzip) to the
Python iterator surface.
"""

import gzip

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _write(path, rows):
    """Write tab-separated rows to a BED file and return its path as str."""
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


def test_read_bed3(tmp_path):
    """A BED3 file yields entries with all optional fields None."""
    pg = _pg()
    path = _write(tmp_path / "a.bed", [
        ("chr1", 1000, 2000),
        ("chr2", 5000, 6000),
        ("chrX", 100, 500),
    ])
    entries = list(pg.BedReader(path))
    assert [(e.chrom, e.start, e.end) for e in entries] == [
        ("chr1", 1000, 2000), ("chr2", 5000, 6000), ("chrX", 100, 500),
    ]
    assert entries[0].name is None
    assert entries[0].score is None
    assert entries[0].strand is None


def test_read_bed6(tmp_path):
    """A BED6 file yields name/score/strand on each entry."""
    pg = _pg()
    path = _write(tmp_path / "a.bed", [
        ("chr1", 1000, 2000, "feature1", 500, "+"),
        ("chr2", 3000, 4000, "feature2", 200, "-"),
        ("chrX", 100, 500, "feature3", 0, "."),
    ])
    entries = list(pg.BedReader(path))
    assert len(entries) == 3
    assert entries[0].name == "feature1"
    assert entries[0].score == 500
    assert entries[0].strand == "+"
    assert entries[1].strand == "-"
    assert entries[2].strand == "."


def test_read_bed12(tmp_path):
    """A BED12 line yields thickness / item_rgb / blocks. Mirrors readBED12Format."""
    pg = _pg()
    path = _write(tmp_path / "b12.bed", [
        ("chr1", 1000, 2000, "item1", 100, "+", 1200, 1800, "255,0,0",
         2, "400,400", "0,600"),
    ])
    e = list(pg.BedReader(path))[0]
    assert e.name == "item1"
    assert e.score == 100
    assert e.strand == "+"
    assert (e.thickness.start, e.thickness.end) == (1200, 1800)
    assert (e.item_rgb.red, e.item_rgb.green, e.item_rgb.blue) == (255, 0, 0)
    assert e.blocks.count == 2
    assert list(e.blocks.sizes) == [400, 400]
    assert list(e.blocks.starts) == [0, 600]


def test_mixed_formats_no_stale_optionals(tmp_path):
    """Optional fields from a previous record must not leak into the next.

    Mirrors mixedBedFormatsNoStaleOptionals: BED12, then BED3, then BED6.
    """
    pg = _pg()
    path = _write(tmp_path / "mixed.bed", [
        ("chr1", 1000, 2000, "item1", 100, "+", 1200, 1800, "255,0,0",
         2, "400,400", "0,600"),
        ("chr2", 5000, 6000),
        ("chr3", 3000, 4000, "feature3", 500, "-"),
    ])
    entries = list(pg.BedReader(path))
    assert len(entries) == 3

    # BED12 — all optionals present
    assert entries[0].name == "item1"
    assert entries[0].thickness is not None
    assert entries[0].item_rgb is not None
    assert entries[0].blocks is not None

    # BED3 — every optional reset to None
    assert entries[1].chrom == "chr2"
    assert entries[1].name is None
    assert entries[1].score is None
    assert entries[1].strand is None
    assert entries[1].thickness is None
    assert entries[1].item_rgb is None
    assert entries[1].blocks is None

    # BED6 — name/score/strand only; the BED12-only fields stay None
    assert entries[2].name == "feature3"
    assert entries[2].strand == "-"
    assert entries[2].thickness is None
    assert entries[2].item_rgb is None
    assert entries[2].blocks is None


def test_only_comments(tmp_path):
    """A comments-only file is valid and yields zero entries. Mirrors validationOnlyComments."""
    pg = _pg()
    path = tmp_path / "c.bed"
    path.write_text("# Comment 1\n# Comment 2\n# Comment 3\n")
    reader = pg.BedReader(str(path))
    assert list(reader) == []
    assert reader.get_error_message() == ""


def test_embedded_nul_in_field_preserved(tmp_path):
    """A NUL embedded in the name column survives into the parsed entry and
    does not truncate the line. Mirrors embeddedNulInFieldPreserved."""
    pg = _pg()
    path = tmp_path / "nul.bed"
    path.write_bytes(b"chr1\t100\t200\tna\x00me\t0\t+\n")
    e = list(pg.BedReader(str(path)))[0]
    assert e.name is not None
    assert len(e.name) == 5          # 'n','a','\0','m','e'
    assert e.name[2] == "\x00"
    # trailing fields still parsed -> line wasn't cut at the NUL
    assert e.score == 0
    assert e.strand == "+"


def test_iterator_is_single_pass(tmp_path):
    """Iterating yields each record once; the reader does not restart."""
    pg = _pg()
    path = _write(tmp_path / "a.bed", [("chr1", 0, 10), ("chr1", 20, 30)])
    reader = pg.BedReader(path)
    assert iter(reader) is reader            # __iter__ returns self
    first = list(reader)
    assert len(first) == 2
    # already drained: a second pass yields nothing (single-pass)
    assert list(reader) == []


def test_empty_file(tmp_path):
    """An empty BED file yields no entries."""
    pg = _pg()
    path = str(tmp_path / "empty.bed")
    (tmp_path / "empty.bed").write_text("")
    assert list(pg.BedReader(path)) == []


def test_invalid_line_raises_by_default(tmp_path):
    """With skip_invalid_lines=False (default), a malformed line raises."""
    pg = _pg()
    path = _write(tmp_path / "bad.bed", [
        ("chr1", 1000, 2000),
        ("chr2", "NOTANUMBER", 3000),    # non-numeric coordinate
    ])
    with pytest.raises(RuntimeError):
        list(pg.BedReader(path))


def test_skip_invalid_lines(tmp_path):
    """With skip_invalid_lines=True, malformed (non-first) lines are skipped."""
    pg = _pg()
    path = _write(tmp_path / "bad.bed", [
        ("chr1", 1000, 2000),
        ("chr2", "NOTANUMBER", 3000),
        ("chr3", 4000, 5000),
    ])
    reader = pg.BedReader(path, skip_invalid_lines=True)
    entries = list(reader)
    assert [e.chrom for e in entries] == ["chr1", "chr3"]
    # after a full pass the error message reflects the last read (clean EOF)
    assert reader.get_error_message() == ""


def test_invalid_first_line_raises_even_when_skipping(tmp_path):
    """The first data record is validated when the reader is constructed, so a
    malformed FIRST line raises immediately — even with skip_invalid_lines=True.
    Mirrors genogrove validationInvalidCoordinates (constructor validates line 1).
    """
    pg = _pg()
    path = _write(tmp_path / "badfirst.bed", [("chr1", "INVALID", 2000)])
    with pytest.raises(RuntimeError):
        pg.BedReader(path)
    with pytest.raises(RuntimeError):
        pg.BedReader(path, skip_invalid_lines=True)


def test_file_not_found(tmp_path):
    """Opening a missing file raises."""
    pg = _pg()
    with pytest.raises((RuntimeError, IOError, OSError)):
        pg.BedReader(str(tmp_path / "does_not_exist.bed"))


def test_current_line_advances(tmp_path):
    """get_current_line tracks input consumed (0 before the first read)."""
    pg = _pg()
    path = _write(tmp_path / "a.bed", [("chr1", 0, 10), ("chr1", 20, 30)])
    reader = pg.BedReader(path)
    assert reader.get_current_line() == 0
    list(reader)
    assert reader.get_current_line() >= 2


def test_reads_gzip(tmp_path):
    """A gzip-compressed BED file is decompressed transparently."""
    pg = _pg()
    path = tmp_path / "a.bed.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("chr1\t1000\t2000\tfeat\t500\t+\n")
    entries = list(pg.BedReader(str(path)))
    assert len(entries) == 1
    assert entries[0].chrom == "chr1"
    assert entries[0].name == "feat"


def test_build_grove_from_reader(tmp_path):
    """The common workflow: read a BED file into a BedGrove.

    BED coordinates are 0-based half-open [start, end); the grove key is the
    closed interval [start, end - 1].
    """
    pg = _pg()
    path = _write(tmp_path / "a.bed", [
        ("chr1", 1000, 2000, "geneA", 100, "+"),
        ("chr1", 3000, 4000, "geneB", 200, "-"),
    ])
    g = pg.BedGrove(64)
    for e in pg.BedReader(path):
        g.insert(e.chrom, pg.Interval(e.start, e.end - 1), e)

    assert g.size() == 2
    hits = list(g.intersect(pg.Interval(1500, 1500), "chr1"))
    assert len(hits) == 1
    assert hits[0].data.name == "geneA"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
