"""
SIF (Simple Interaction Format) export — Grove.to_sif(path) (#34).

A thin pass-through over genogrove's node-less grove_to_sif(ostream) (v0.24.7).
It writes the grove's B+ tree structure (nodelink / leaflink lines) and the
graph-overlay edges (keylink lines) as tab-separated interactions for graph
visualization (e.g. Cytoscape).

SIF line/index order is NOT stable across runs (hash-map index iteration), so
these assert on sets of lines, never exact order. The keylink token for a key
is its value's to_string — i.e. exactly `str(key.value)` — so expected lines are
derived from the keys themselves rather than hard-coding a format.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _read_lines(path):
    with open(path) as f:
        return [ln for ln in f.read().splitlines() if ln.strip()]


def test_to_sif_writes_graph_edges(tmp_path):
    pg = _pg()
    g = pg.Grove(8)  # order high enough that 3 keys stay in one leaf (no splits)
    a = g.insert("chr1", pg.GenomicCoordinate(".", 100, 200))
    b = g.insert("chr1", pg.GenomicCoordinate(".", 300, 400))
    c = g.insert("chr1", pg.GenomicCoordinate(".", 500, 600))
    g.add_edge(a, b)
    g.add_edge(a, c)

    path = str(tmp_path / "graph.sif")
    g.to_sif(path)

    keylinks = {ln for ln in _read_lines(path) if "\tkeylink\t" in ln}
    assert keylinks == {
        f"{a.value}\tkeylink\t{b.value}",
        f"{a.value}\tkeylink\t{c.value}",
    }


def test_to_sif_empty_grove_writes_empty_file(tmp_path):
    pg = _pg()
    g = pg.Grove(3)
    path = str(tmp_path / "empty.sif")
    g.to_sif(path)
    with open(path) as f:
        assert f.read() == ""


def test_to_sif_includes_tree_structure_after_splits(tmp_path):
    """Enough keys to split → SIF carries the B+ tree structure links too."""
    pg = _pg()
    g = pg.Grove(3)  # small order forces node splits
    for i in range(20):
        g.insert("chr1", pg.GenomicCoordinate(".", i * 10, i * 10 + 5))

    path = str(tmp_path / "tree.sif")
    g.to_sif(path)
    text = open(path).read()
    assert "nodelink" in text   # internal-node → child links
    assert "leaflink" in text   # leaf-chain links


def test_to_sif_is_generic_over_key_type(tmp_path):
    """to_sif works on any grove — here a NumericGrove (numeric key to_string)."""
    pg = _pg()
    g = pg.NumericGrove(8)
    a = g.insert("ids", pg.Numeric(1))
    b = g.insert("ids", pg.Numeric(2))
    g.add_edge(a, b)

    path = str(tmp_path / "num.sif")
    g.to_sif(path)
    keylinks = {ln for ln in _read_lines(path) if "\tkeylink\t" in ln}
    assert keylinks == {f"{a.value}\tkeylink\t{b.value}"}


def test_to_sif_unwritable_path_raises(tmp_path):
    pg = _pg()
    g = pg.Grove(3)
    g.insert("chr1", pg.GenomicCoordinate(".", 10, 20))
    bad = str(tmp_path / "no_such_dir" / "out.sif")
    with pytest.raises((RuntimeError, IOError, OSError)):
        g.to_sif(bad)
