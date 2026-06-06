"""
Tests for the BED value types (BedEntry + BlockInfo / ThickInfo / RgbColor).

Mirrors the data carriers from genogrove/include/genogrove/io/bed_reader.hpp.
The bed_reader file iterator itself is not bound yet; these cover only the
entry value types exposed for the data-carrying grove.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_construct_minimal():
    """BedEntry(chrom, start, end) stores the core BED3 fields."""
    pg = _pg()
    e = pg.BedEntry("chr1", 100, 200)
    assert e.chrom == "chr1"
    assert e.start == 100
    assert e.end == 200


def test_optionals_default_to_none():
    """All optional BED4+ fields are None on a fresh entry."""
    pg = _pg()
    e = pg.BedEntry("chr1", 100, 200)
    assert e.name is None
    assert e.score is None
    assert e.strand is None
    assert e.thickness is None
    assert e.item_rgb is None
    assert e.blocks is None


def test_optional_scalar_round_trip():
    """Setting name/score/strand reads back the assigned value."""
    pg = _pg()
    e = pg.BedEntry("chr1", 100, 200)
    e.name = "gene1"
    e.score = 500
    e.strand = "+"
    assert e.name == "gene1"
    assert e.score == 500
    assert e.strand == "+"


def test_optional_cleared_with_none():
    """Assigning None clears a previously-set optional."""
    pg = _pg()
    e = pg.BedEntry("chr1", 100, 200)
    e.name = "gene1"
    e.name = None
    assert e.name is None


def test_strand_must_be_single_char():
    """strand maps to a C++ char: empty / multi-char strings raise ValueError."""
    pg = _pg()
    e = pg.BedEntry("chr1", 100, 200)
    with pytest.raises(ValueError):
        e.strand = ""
    with pytest.raises(ValueError):
        e.strand = "++"


def test_fields_are_mutable():
    """Core fields are read/write (BedEntry is a plain data carrier, not a key)."""
    pg = _pg()
    e = pg.BedEntry("chr1", 100, 200)
    e.chrom = "chr2"
    e.start = 10
    e.end = 20
    assert (e.chrom, e.start, e.end) == ("chr2", 10, 20)


def test_rgb_color():
    """RgbColor channels are ints in [0, 255] and assignable into item_rgb."""
    pg = _pg()
    c = pg.RgbColor(255, 128, 0)
    assert (c.red, c.green, c.blue) == (255, 128, 0)

    e = pg.BedEntry("chr1", 100, 200)
    e.item_rgb = pg.RgbColor(1, 2, 3)
    assert (e.item_rgb.red, e.item_rgb.green, e.item_rgb.blue) == (1, 2, 3)


def test_thick_info():
    """ThickInfo(start, end) round-trips and assigns into thickness."""
    pg = _pg()
    t = pg.ThickInfo(120, 180)
    assert t.start == 120
    assert t.end == 180

    e = pg.BedEntry("chr1", 100, 200)
    e.thickness = pg.ThickInfo(110, 190)
    assert e.thickness.start == 110
    assert e.thickness.end == 190


def test_block_info():
    """BlockInfo carries count + parallel size/start lists (list[int])."""
    pg = _pg()
    b = pg.BlockInfo(2, [10, 20], [0, 30])
    assert b.count == 2
    assert list(b.sizes) == [10, 20]
    assert list(b.starts) == [0, 30]

    e = pg.BedEntry("chr1", 100, 200)
    e.blocks = pg.BlockInfo(1, [50], [0])
    assert e.blocks.count == 1
    assert list(e.blocks.sizes) == [50]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])