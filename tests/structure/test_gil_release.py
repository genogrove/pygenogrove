"""
Concurrency smoke tests for the GIL-releasing I/O bindings (issue #31, #83).

`serialize` / `deserialize`, bulk insert, the reader `__next__`, and
`FastaIndex` now release the GIL around their pure-C++ / htslib work. `intersect`,
`flanking`, `remove_key`, `compact`, and `Registry.serialize` do too (#83). These
tests don't try to prove overlap (timing-dependent); they verify those paths
stay correct when driven concurrently from many Python threads — i.e. releasing
the GIL didn't corrupt anything or break the return / exception paths. A
misplaced release on something that touches Python would crash or corrupt here.

Each thread uses its OWN objects: a reader/grove isn't shared-thread-safe; the
realistic pattern is one driver per object, overlapping only the C++ work.
`Registry` is a process-wide singleton, so its test instead overlaps read-only
(mutex-protected) `serialize()` calls against one shared instance.
"""

from concurrent.futures import ThreadPoolExecutor

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def test_concurrent_serialize_deserialize(tmp_path):
    pg = _pg()

    def roundtrip(i):
        g = pg.Grove(4)
        for j in range(50):
            g.insert("chr1", pg.GenomicCoordinate(".", j * 10, j * 10 + 5),
                     {"i": i, "j": j})
        path = str(tmp_path / f"g{i}.gg")
        g.serialize(path)
        loaded = pg.Grove.deserialize(path)
        assert loaded.size() == 50
        # Assert the exact hit count — a corrupt/cross-contaminated result from a
        # mis-placed GIL release would surface as an extra or missing hit here.
        hits = list(loaded.intersect(pg.GenomicCoordinate(".", 20, 21), "chr1"))
        assert len(hits) == 1
        return hits[0].data

    with ThreadPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(roundtrip, range(8)))

    assert [r["j"] for r in results] == [2] * 8          # [20,25] is the only hit
    assert sorted(r["i"] for r in results) == list(range(8))


def test_concurrent_bulk_insert():
    pg = _pg()

    def build(_):
        g = pg.Grove(3)  # small order -> many splits in the released region
        items = [(pg.GenomicCoordinate(".", k * 10, k * 10 + 5), {"k": k})
                 for k in range(100)]
        keys = g.insert_bulk("chr1", items)
        assert len(keys) == 100
        return g.size()

    with ThreadPoolExecutor(max_workers=8) as ex:
        sizes = list(ex.map(build, range(8)))

    assert sizes == [100] * 8


def test_concurrent_reader(tmp_path):
    pg = _pg()
    bed = tmp_path / "x.bed"
    bed.write_text("".join(f"chr1\t{j * 10}\t{j * 10 + 5}\tf{j}\n" for j in range(100)))
    expected = [f"f{j}" for j in range(100)]

    def read(_):
        # Each thread drives its own reader (readers are single-pass, not shared).
        return [e.name for e in pg.BedReader(str(bed))]

    with ThreadPoolExecutor(max_workers=8) as ex:
        all_names = list(ex.map(read, range(8)))

    assert all(names == expected for names in all_names)


def test_concurrent_flanking_remove_compact():
    pg = _pg()

    def build_and_shrink(_):
        g = pg.Grove(3)
        for j in range(50):
            g.insert("chr1", pg.GenomicCoordinate(".", j * 10, j * 10 + 5),
                     {"j": j})
        # [96,97] falls in the gap between j=9's [90,95] and j=10's [100,105] —
        # overlaps neither, so both flank it cleanly.
        res = g.flanking(pg.GenomicCoordinate(".", 96, 97), "chr1")
        assert res.predecessor.data["j"] == 9
        assert res.successor.data["j"] == 10

        removed = g.remove_key("chr1", res.predecessor)
        assert removed
        g.compact()
        return g.size()

    with ThreadPoolExecutor(max_workers=8) as ex:
        sizes = list(ex.map(build_and_shrink, range(8)))

    assert sizes == [49] * 8


def test_concurrent_registry_serialize(tmp_path):
    pg = _pg()
    reg = pg.Registry.instance()
    reg.clear()
    for name in ("chr1", "chr2", "chr3"):
        reg.intern(name)

    def dump(i):
        path = str(tmp_path / f"r{i}.bin")
        reg.serialize(path)
        return (tmp_path / f"r{i}.bin").read_bytes()

    with ThreadPoolExecutor(max_workers=8) as ex:
        blobs = list(ex.map(dump, range(8)))

    assert all(b == blobs[0] for b in blobs)
    reg.clear()


def test_concurrent_fasta_fetch(tmp_path):
    pg = _pg()
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\nACGTACGTACGT\n>chr2\nTTTTGGGGCCCC\n")
    # Build the .fai once up front so the threads don't race to create it; each
    # then opens its own handle and only reads (faidx fetch is the released path).
    pg.FastaIndex(str(fa))

    def fetch(_):
        idx = pg.FastaIndex(str(fa))
        return idx.fetch("chr1", 0, 4), idx.fetch("chr2")

    with ThreadPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(fetch, range(8)))

    assert all(r == ("ACGT", "TTTTGGGGCCCC") for r in results)
