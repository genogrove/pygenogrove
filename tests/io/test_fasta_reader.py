"""
Tests for the FASTA/FASTQ sequence reader (FastaReader -> FastaEntry).

FASTA records are named sequences (not intervals), so this reader is standalone.
Covers FASTA + FASTQ parsing (name / comment / sequence / quality), the
skip_empty_sequences option, the FastaEntry value type, and error handling.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


_FASTA = ">seq1 my description\nACGTACGT\n>seq2\nTTT\n"
_FASTQ = "@read1 a comment\nACGT\n+\nIIII\n"
_FASTA_WITH_EMPTY = ">full\nACGT\n>empty\n>another\nGG\n"


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content)
    return str(p)


def test_read_fasta(tmp_path):
    pg = _pg()
    recs = list(pg.FastaReader(_write(tmp_path, "x.fa", _FASTA)))
    assert [r.name for r in recs] == ["seq1", "seq2"]
    r1 = recs[0]
    assert r1.comment == "my description"
    assert r1.sequence == "ACGTACGT"
    assert len(r1) == 8
    assert r1.quality is None
    assert r1.is_fastq() is False
    assert recs[1].comment == ""
    assert recs[1].sequence == "TTT"


def test_read_fastq(tmp_path):
    pg = _pg()
    recs = list(pg.FastaReader(_write(tmp_path, "x.fq", _FASTQ)))
    assert len(recs) == 1
    r = recs[0]
    assert r.name == "read1"
    assert r.comment == "a comment"
    assert r.sequence == "ACGT"
    assert r.quality == "IIII"
    assert r.is_fastq() is True


def test_skip_empty_sequences(tmp_path):
    pg = _pg()
    path = _write(tmp_path, "e.fa", _FASTA_WITH_EMPTY)
    # default keeps the empty record
    assert [r.name for r in pg.FastaReader(path)] == ["full", "empty", "another"]
    # skip it
    kept = [r.name for r in pg.FastaReader(path, skip_empty_sequences=True)]
    assert kept == ["full", "another"]


def test_fasta_entry_construct():
    pg = _pg()
    e = pg.FastaEntry("seq1", "ACGT")
    assert e.name == "seq1"
    assert e.sequence == "ACGT"
    assert len(e) == 4
    assert e.quality is None


def test_repr(tmp_path):
    pg = _pg()
    reader = pg.FastaReader(_write(tmp_path, "x.fa", _FASTA))
    assert repr(reader) == "FastaReader(line=0)"
    list(reader)
    assert repr(reader) == f"FastaReader(line={reader.get_current_line()})"


def test_concurrent_next_raises(tmp_path):
    """__next__ guards against concurrent calls on the same reader (issue #84).
    Pokes the internal flag directly rather than racing real threads — the
    guard's check-then-set only needs to be correct, not literally raced."""
    pg = _pg()
    reader = pg.FastaReader(_write(tmp_path, "x.fa", _FASTA))
    reader._pygg_reading = True
    with pytest.raises(RuntimeError):
        next(reader)
    # The guard doesn't get stuck: a normal call after the flag is cleared works.
    reader._pygg_reading = False
    assert next(reader) is not None


def test_clean_eof_error_message(tmp_path):
    pg = _pg()
    r = pg.FastaReader(_write(tmp_path, "x.fa", _FASTA))
    list(r)
    assert r.get_error_message() == ""


def test_missing_file_raises():
    pg = _pg()
    with pytest.raises((RuntimeError, IOError, OSError)):
        pg.FastaReader("/nonexistent_dir_xyz/genome.fa")