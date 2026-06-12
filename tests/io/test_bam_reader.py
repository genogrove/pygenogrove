"""
Tests for the SAM/BAM alignment reader (BamReader -> SamEntry).

Uses a hand-written SAM file (htslib auto-detects SAM; no index needed). Covers
field/strand/flag access, the filtering options, SamEntry.to_coordinate() /
.to_dict(), and the "load alignments into the universal Grove" flow.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


# read1: forward, primary, mapq 60, POS 101 (1-based) -> 0-based [100, 104)
# read2: reverse (FLAG 16), mapq 30, POS 201 -> [200, 204)
# read3: unmapped (FLAG 4)
_SAM = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr1\tLN:100000\n"
    "read1\t0\tchr1\t101\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
    "read2\t16\tchr1\t201\t30\t4M\t*\t0\t0\tACGT\tIIII\n"
    "read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
)


def _sam(tmp_path):
    p = tmp_path / "reads.sam"
    p.write_text(_SAM)
    return str(p)


def test_read_mapped_alignments(tmp_path):
    pg = _pg()
    alns = list(pg.BamReader(_sam(tmp_path)))  # skip_unmapped defaults True
    assert [a.qname for a in alns] == ["read1", "read2"]
    a1 = alns[0]
    assert a1.chrom == "chr1"
    assert a1.start == 100 and a1.end == 104
    assert a1.mapq == 60
    assert a1.get_strand() == "+"
    assert a1.cigar == "4M"
    assert alns[1].get_strand() == "-"  # FLAG 16 (reverse)


def test_alignment_flags(tmp_path):
    pg = _pg()
    a1, a2 = list(pg.BamReader(_sam(tmp_path)))
    assert a1.flags.is_reverse() is False
    assert a2.flags.is_reverse() is True
    assert a1.flags.has_flag(pg.SamFlags.REVERSE) is False
    assert a2.flags.has_flag(pg.SamFlags.REVERSE) is True
    assert a1.is_primary() and a1.is_mapped()
    assert a1.flags.value() == 0 and a2.flags.value() == 16


def test_to_coordinate_strand_aware(tmp_path):
    pg = _pg()
    a1, a2 = list(pg.BamReader(_sam(tmp_path)))
    assert a1.to_coordinate() == pg.GenomicCoordinate("+", 100, 103)
    assert a2.to_coordinate() == pg.GenomicCoordinate("-", 200, 203)


def test_to_coordinate_raises_for_unmapped(tmp_path):
    pg = _pg()
    alns = list(pg.BamReader(_sam(tmp_path), skip_unmapped=False))
    unmapped = [a for a in alns if not a.is_mapped()]
    assert len(unmapped) == 1
    assert unmapped[0].consumes_reference() is False
    with pytest.raises(ValueError):
        unmapped[0].to_coordinate()


def test_skip_unmapped_option(tmp_path):
    pg = _pg()
    assert len(list(pg.BamReader(_sam(tmp_path)))) == 2  # default skips unmapped
    assert len(list(pg.BamReader(_sam(tmp_path), skip_unmapped=False))) == 3


def test_min_mapq_filter(tmp_path):
    pg = _pg()
    # read1 mapq 60, read2 mapq 30 -> min_mapq=50 keeps only read1
    alns = list(pg.BamReader(_sam(tmp_path), min_mapq=50))
    assert [a.qname for a in alns] == ["read1"]


def test_to_dict(tmp_path):
    pg = _pg()
    d = list(pg.BamReader(_sam(tmp_path)))[0].to_dict()
    assert d["qname"] == "read1"
    assert d["mapq"] == 60
    assert d["strand"] == "+"
    assert d["cigar"] == "4M"
    assert d["is_primary"] is True
    assert d["is_mapped"] is True


def test_load_bam_into_grove(tmp_path):
    """The intended flow: alignments -> universal Grove via to_coordinate/to_dict."""
    pg = _pg()
    g = pg.Grove()
    for aln in pg.BamReader(_sam(tmp_path)):
        g.insert(aln.chrom, aln.to_coordinate(), aln.to_dict())
    assert g.size() == 2

    # strand-aware: a '+' query matches only the forward read
    hits = list(g.intersect(pg.GenomicCoordinate("+", 100, 103), "chr1"))
    assert len(hits) == 1
    assert hits[0].data["qname"] == "read1"
    assert hits[0].data["strand"] == "+"


def test_clean_eof_error_message(tmp_path):
    pg = _pg()
    r = pg.BamReader(_sam(tmp_path))
    list(r)
    assert r.get_error_message() == ""


def test_missing_file_raises():
    pg = _pg()
    with pytest.raises((RuntimeError, IOError, OSError)):
        pg.BamReader("/nonexistent_dir_xyz/reads.bam")