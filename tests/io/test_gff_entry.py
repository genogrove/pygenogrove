"""
Tests for the GFF/GTF value types (GffEntry + GffFormat).

Mirrors the data carriers from genogrove/include/genogrove/io/gff_reader.hpp.
The gff_reader file iterator itself is not bound yet; these cover only the
entry value type exposed for the data-carrying grove.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_construct_minimal():
    """GffEntry(seqid, start, end, type) stores the core columns."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "gene")
    assert e.seqid == "chr1"
    assert e.start == 1000
    assert e.end == 2000
    assert e.type == "gene"


def test_optionals_default_to_none():
    """Optional GFF columns are None on a fresh entry; attributes is empty."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "gene")
    assert e.score is None
    assert e.strand is None
    assert e.phase is None
    assert dict(e.attributes) == {}


def test_optional_scalar_round_trip():
    """score (float), strand (char), phase (int) read back the assigned value."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "CDS")
    e.score = 3.5
    e.strand = "+"
    e.phase = 2
    assert e.score == 3.5
    assert e.strand == "+"
    assert e.phase == 2


def test_strand_must_be_single_char():
    """strand maps to a C++ char: empty / multi-char strings raise ValueError."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "gene")
    with pytest.raises(ValueError):
        e.strand = ""
    with pytest.raises(ValueError):
        e.strand = "+-"


def test_attributes_dict_round_trip():
    """attributes is a dict[str, str] (column 9)."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "transcript")
    e.attributes = {"gene_id": "ENSG1", "transcript_id": "ENST1"}
    assert dict(e.attributes) == {"gene_id": "ENSG1", "transcript_id": "ENST1"}


def test_gtf_helper_methods():
    """GTF helper getters read the GTF attribute keys.

    Mirrors genogrove gfffile-test.cpp::gtfHelperMethods (gene then exon).
    """
    pg = _pg()
    gene = pg.GffEntry("chr1", 1000, 2000, "gene")
    gene.attributes = {
        "gene_id": "ENSG00000001",
        "gene_name": "TEST1",
        "gene_biotype": "protein_coding",
    }
    assert gene.get_gene_id() == "ENSG00000001"
    assert gene.get_gene_name() == "TEST1"
    assert gene.get_gene_biotype() == "protein_coding"
    assert gene.get_transcript_id() is None   # a gene has no transcript_id

    exon = pg.GffEntry("chr1", 1000, 1100, "exon")
    exon.attributes = {
        "gene_id": "ENSG00000001",
        "transcript_id": "ENST00000001",
        "exon_number": "1",
    }
    assert exon.get_gene_id() == "ENSG00000001"
    assert exon.get_transcript_id() == "ENST00000001"
    assert exon.get_exon_number() == 1        # parsed to int


def test_gff3_helper_methods():
    """GFF3 uses different attribute keys; the helpers fall back accordingly.

    Mirrors genogrove gfffile-test.cpp::gff3HelperMethods.
    """
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "gene")
    e.attributes = {"ID": "gene1", "Name": "TEST1", "biotype": "protein_coding"}
    # GFF3 has no gene_id key.
    assert e.get_gene_id() is None
    # gene_name falls back to the GFF3 "Name" attribute.
    assert e.get_gene_name() == "TEST1"
    # gene_biotype falls back to the GFF3 "biotype" attribute.
    assert e.get_gene_biotype() == "protein_coding"
    # generic getter reaches the raw column-9 keys.
    assert e.get_attribute("ID") == "gene1"


def test_generic_attribute_getter():
    """get_attribute returns the value or None. Mirrors genericAttributeGetter."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "gene")
    e.attributes = {"gene_id": "ENSG00000001"}
    assert e.get_attribute("gene_id") == "ENSG00000001"
    assert e.get_attribute("nonexistent") is None


def test_format_enum_and_predicates():
    """format is a GffFormat enum; is_gtf/is_gff3 reflect it."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "gene")
    assert e.format == pg.GffFormat.UNKNOWN
    assert not e.is_gtf()
    assert not e.is_gff3()

    e.format = pg.GffFormat.GTF
    assert e.is_gtf()
    assert not e.is_gff3()

    e.format = pg.GffFormat.GFF3
    assert e.is_gff3()
    assert not e.is_gtf()


def test_fields_are_mutable():
    """Core fields are read/write (GffEntry is a plain data carrier, not a key)."""
    pg = _pg()
    e = pg.GffEntry("chr1", 1000, 2000, "gene")
    e.seqid = "chr2"
    e.source = "ensembl"
    e.start = 5
    e.end = 50
    assert (e.seqid, e.source, e.start, e.end) == ("chr2", "ensembl", 5, 50)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
