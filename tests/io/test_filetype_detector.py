"""
Tests for the file-type detector (FiletypeDetector + Filetype / CompressionType).

detect_filetype(path) -> (Filetype, CompressionType), inferred from the path's
extension (compression extension stripped first) and magic bytes.
"""

import gzip

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _detect(pg, path):
    return pg.FiletypeDetector().detect_filetype(str(path))


def test_detect_plain_bed(tmp_path):
    pg = _pg()
    p = tmp_path / "peaks.bed"
    p.write_text("chr1\t100\t200\n")
    ft, comp = _detect(pg, p)
    assert ft == pg.Filetype.BED
    assert comp == pg.CompressionType.NONE


def test_detect_gzip_bed(tmp_path):
    pg = _pg()
    p = tmp_path / "peaks.bed.gz"
    with gzip.open(p, "wt") as f:
        f.write("chr1\t100\t200\n")
    ft, comp = _detect(pg, p)
    assert ft == pg.Filetype.BED  # compression extension stripped -> .bed
    assert comp == pg.CompressionType.GZIP


def test_detect_fasta(tmp_path):
    pg = _pg()
    p = tmp_path / "genome.fa"
    p.write_text(">s\nACGT\n")
    ft, comp = _detect(pg, p)
    assert ft == pg.Filetype.FASTA
    assert comp == pg.CompressionType.NONE


def test_detect_gff(tmp_path):
    pg = _pg()
    p = tmp_path / "genes.gff"
    p.write_text("chr1\tx\tgene\t1\t100\t.\t+\t.\tID=g1\n")
    ft, _ = _detect(pg, p)
    assert ft == pg.Filetype.GFF


def test_detect_gg(tmp_path):
    pg = _pg()
    p = tmp_path / "grove.gg"
    p.write_bytes(b"\x00\x01\x02")
    ft, _ = _detect(pg, p)
    assert ft == pg.Filetype.GG


def test_detect_unknown(tmp_path):
    pg = _pg()
    p = tmp_path / "mystery.xyz"
    p.write_text("nope")
    ft, _ = _detect(pg, p)
    assert ft == pg.Filetype.UNKNOWN


def test_returns_two_tuple_of_enums(tmp_path):
    pg = _pg()
    p = tmp_path / "x.bed"
    p.write_text("chr1\t1\t2\n")
    result = _detect(pg, p)
    assert isinstance(result, tuple) and len(result) == 2
    ft, comp = result
    assert isinstance(ft, pg.Filetype)
    assert isinstance(comp, pg.CompressionType)