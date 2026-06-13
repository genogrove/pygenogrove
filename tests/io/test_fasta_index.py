"""
Tests for the random-access FASTA index (FastaIndex, htslib faidx).

FastaIndex(path) opens a FASTA and loads/creates its .fai index, then serves
region/whole-sequence fetches and per-sequence metadata. Region coordinates are
0-based half-open [start, end); the index is built on first open (writes a
sibling .fai), so tmp_path must be writable — which it is.

These port the behavioural faidx cases (region fetch, whole-sequence fetch,
metadata, membership, error handling) and the GenomicCoordinate pairing.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


# chr1 = 20 bp, chr2 = 16 bp; single-line records keep faidx happy.
_FASTA = ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nTTTTAAAACCCCGGGG\n"


def _index(pg, tmp_path):
    p = tmp_path / "genome.fa"
    p.write_text(_FASTA)
    return pg.FastaIndex(str(p))


def test_fetch_region(tmp_path):
    pg = _pg()
    idx = _index(pg, tmp_path)
    # 0-based half-open [start, end).
    assert idx.fetch("chr1", 0, 4) == "ACGT"
    assert idx.fetch("chr1", 4, 8) == "ACGT"
    assert idx.fetch("chr2", 0, 4) == "TTTT"


def test_fetch_whole_sequence(tmp_path):
    pg = _pg()
    idx = _index(pg, tmp_path)
    assert idx.fetch("chr1") == "ACGTACGTACGTACGTACGT"
    assert idx.fetch("chr2") == "TTTTAAAACCCCGGGG"


def test_sequence_metadata(tmp_path):
    pg = _pg()
    idx = _index(pg, tmp_path)
    assert idx.sequence_count() == 2
    assert idx.sequence_name(0) == "chr1"
    assert idx.sequence_name(1) == "chr2"
    assert idx.sequence_length("chr1") == 20
    assert idx.sequence_length("chr2") == 16


def test_membership_and_names(tmp_path):
    pg = _pg()
    idx = _index(pg, tmp_path)
    assert idx.has_sequence("chr1") is True
    assert idx.has_sequence("chrX") is False
    # Pythonic surface.
    assert len(idx) == 2
    assert "chr2" in idx
    assert "chrX" not in idx
    assert idx.names() == ["chr1", "chr2"]


def test_fetch_unknown_sequence_raises(tmp_path):
    pg = _pg()
    idx = _index(pg, tmp_path)
    with pytest.raises(IndexError):
        idx.fetch("chrX", 0, 4)
    with pytest.raises(IndexError):
        idx.fetch("chrX")
    with pytest.raises(IndexError):
        idx.sequence_length("chrX")


def test_invalid_region_raises(tmp_path):
    pg = _pg()
    idx = _index(pg, tmp_path)
    with pytest.raises(IndexError):
        idx.fetch("chr1", 8, 4)  # start >= end


def test_sequence_name_out_of_range_raises(tmp_path):
    pg = _pg()
    idx = _index(pg, tmp_path)
    with pytest.raises(IndexError):
        idx.sequence_name(99)


def test_open_missing_file_raises(tmp_path):
    pg = _pg()
    with pytest.raises(RuntimeError):
        pg.FastaIndex(str(tmp_path / "does_not_exist.fa"))


def test_pairs_with_genomic_coordinate(tmp_path):
    """fetch a GenomicCoordinate's bases: gc is 0-based closed, fetch is half-open.

    The chromosome is the grove *index* / FASTA sequence name (passed to fetch);
    a GenomicCoordinate carries (strand, start, end).
    """
    pg = _pg()
    idx = _index(pg, tmp_path)
    # A '+'-strand feature spanning bases [4, 7] of chr1 (closed) -> fetch(name, 4, 8).
    gc = pg.GenomicCoordinate("+", 4, 7)
    assert idx.fetch("chr1", gc.start, gc.end + 1) == "ACGT"