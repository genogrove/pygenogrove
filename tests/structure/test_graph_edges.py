"""
Tests for the completed graph overlay on the universal Grove:

  * labelled edges (#1) — grove<gc, json_value, json_value>: add_edge(s, t, data),
    get_edges, get_neighbors_if, link_with.
  * edge removal / bulk linking (#2, available on every grove): remove_edges_from,
    remove_edges_to, remove_all_edges, clear_graph, graph_empty, link_if.

Mirrors genogrove graph_overlay_test.cpp's metadata + link_if cases. The basic
unlabelled edge surface (add_edge/remove_edge/has_edge/get_neighbors/…) is
covered in test_graph_overlay.py.
"""

import gc

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _chain(g, *spans, index="chr1"):
    """Insert ascending (start, end) spans, return the Keys in insertion order."""
    import pygenogrove as pg
    return [g.insert(index, pg.GenomicCoordinate(".", s, e)) for s, e in spans]


# --------------------------------------------------------------------------- #
# #1 Labelled edges (JSON edge metadata)
# --------------------------------------------------------------------------- #

def test_add_edge_with_metadata():
    pg = _pg()
    g = pg.Grove(3)
    a, b = _chain(g, (100, 200), (300, 400))
    g.add_edge(a, b, {"type": "exon", "weight": 7})
    assert g.get_edges(a) == [{"type": "exon", "weight": 7}]
    assert g.has_edge(a, b)


def test_two_arg_add_edge_yields_none_metadata():
    """The unlabelled overload still works on the JSON grove — payload is None."""
    pg = _pg()
    g = pg.Grove(3)
    a, b = _chain(g, (100, 200), (300, 400))
    g.add_edge(a, b)  # no metadata
    assert g.get_edges(a) == [None]


def test_get_edges_parallel_to_neighbors():
    """get_edges(source) is in the same order as get_neighbors(source)."""
    pg = _pg()
    g = pg.Grove(5)
    a, b, c = _chain(g, (100, 200), (300, 400), (500, 600))
    g.add_edge(a, b, {"to": 300})
    g.add_edge(a, c, {"to": 500})
    paired = {n.value.start: e for n, e in zip(g.get_neighbors(a), g.get_edges(a))}
    assert paired == {300: {"to": 300}, 500: {"to": 500}}


def test_get_neighbors_if_filters_by_metadata():
    pg = _pg()
    g = pg.Grove(5)
    a, b, c, d = _chain(g, (100, 200), (300, 400), (500, 600), (700, 800))
    g.add_edge(a, b, {"w": 1})
    g.add_edge(a, c, {"w": 10})
    g.add_edge(a, d, {"w": 20})
    strong = g.get_neighbors_if(a, lambda meta: meta["w"] >= 10)
    assert sorted(n.value.start for n in strong) == [500, 700]


def test_edge_metadata_accepts_all_json_shapes():
    """dict / list / scalar / None all round-trip as edge metadata."""
    pg = _pg()
    g = pg.Grove(7)
    src = g.insert("chr1", pg.GenomicCoordinate(".", 0, 10))
    payloads = [{"k": "v"}, [1, 2, 3], 42, "text", 3.5, True, None]
    targets = _chain(g, *[(100 * i + 100, 100 * i + 150) for i in range(len(payloads))])
    for t, p in zip(targets, payloads):
        g.add_edge(src, t, p)
    assert g.get_edges(src) == payloads


def test_link_with_attaches_metadata_over_a_run():
    pg = _pg()
    g = pg.Grove(5)
    keys = _chain(g, (100, 200), (300, 400), (500, 600))
    # Link each adjacent pair, labelling the edge with the intron gap.
    g.link_with(keys, lambda a, b: {"gap": b.value.start - a.value.end - 1})
    assert g.edge_count() == 2
    assert g.get_edges(keys[0]) == [{"gap": 99}]   # 300 - 200 - 1
    assert g.get_edges(keys[1]) == [{"gap": 99}]   # 500 - 400 - 1


def test_link_with_none_skips_the_edge():
    pg = _pg()
    g = pg.Grove(5)
    keys = _chain(g, (100, 200), (300, 400), (500, 600))
    # Only link the first pair (return None for the rest).
    g.link_with(keys, lambda a, b: {"x": 1} if a.value.start == 100 else None)
    assert g.edge_count() == 1
    assert g.has_edge(keys[0], keys[1])
    assert not g.has_edge(keys[1], keys[2])


def test_serialization_preserves_edge_metadata(tmp_path):
    pg = _pg()
    g = pg.Grove(4)
    a, b, c = _chain(g, (100, 200), (300, 400), (500, 600))
    g.add_edge(a, b, {"type": "exon", "n": 1})
    g.add_edge(a, c, {"type": "intron"})

    path = str(tmp_path / "edges.gg")
    g.serialize(path)
    restored = pg.Grove.deserialize(path)

    # Re-discover the source key by query, then read its labelled edges back.
    src = list(restored.intersect(pg.GenomicCoordinate(".", 100, 200), "chr1"))[0]
    assert restored.edge_count() == 2
    edges = restored.get_edges(src)
    assert {"type": "exon", "n": 1} in edges
    assert {"type": "intron"} in edges


# --------------------------------------------------------------------------- #
# #2 Edge removal / bulk linking (available on every grove)
# --------------------------------------------------------------------------- #

def test_remove_edges_from():
    pg = _pg()
    g = pg.Grove(5)
    a, b, c = _chain(g, (10, 20), (30, 40), (50, 60))
    g.add_edge(a, b)
    g.add_edge(a, c)
    assert g.remove_edges_from(a) == 2
    assert g.out_degree(a) == 0
    assert g.edge_count() == 0


def test_remove_edges_to():
    pg = _pg()
    g = pg.Grove(5)
    a, b, c = _chain(g, (10, 20), (30, 40), (50, 60))
    g.add_edge(a, c)
    g.add_edge(b, c)
    assert g.remove_edges_to(c) == 2
    assert g.edge_count() == 0


def test_remove_all_edges_in_and_out():
    pg = _pg()
    g = pg.Grove(5)
    a, b, c = _chain(g, (10, 20), (30, 40), (50, 60))
    g.add_edge(a, b)   # out of a
    g.add_edge(a, c)   # out of a
    g.add_edge(b, a)   # into a
    assert g.remove_all_edges(a) == 3
    assert g.edge_count() == 0


def test_clear_graph_and_graph_empty():
    pg = _pg()
    g = pg.Grove(3)
    a, b = _chain(g, (10, 20), (30, 40))
    assert g.graph_empty() is True
    g.add_edge(a, b)
    assert g.graph_empty() is False
    g.clear_graph()
    assert g.graph_empty() is True
    assert g.edge_count() == 0
    # keys themselves survive clearing the graph
    assert g.size() == 2


def test_link_if_chains_adjacent_keys():
    pg = _pg()
    g = pg.Grove(5)
    keys = _chain(g, (100, 200), (300, 400), (500, 600), (700, 800))
    g.link_if(keys, lambda a, b: True)
    assert g.edge_count() == len(keys) - 1
    # walk the chain
    node = keys[0]
    seen = [node.value.start]
    for _ in range(3):
        node = g.get_neighbors(node)[0]
        seen.append(node.value.start)
    assert seen == [100, 300, 500, 700]


def test_link_if_predicate_can_skip():
    pg = _pg()
    g = pg.Grove(5)
    keys = _chain(g, (100, 200), (300, 400), (500, 600))
    # Only link pairs whose gap is exactly 99 (all of them here) but flip one off.
    g.link_if(keys, lambda a, b: a.value.start != 300)
    assert g.has_edge(keys[0], keys[1])
    assert not g.has_edge(keys[1], keys[2])
    assert g.edge_count() == 1


# --------------------------------------------------------------------------- #
# remove_edges_if — predicate-filtered edge removal (#33)
# --------------------------------------------------------------------------- #

def test_remove_edges_if_by_metadata():
    """JSON Grove: predicate(target, metadata) removes edges whose metadata matches."""
    pg = _pg()
    g = pg.Grove(5)
    a, b, c, d = _chain(g, (100, 200), (300, 400), (500, 600), (700, 800))
    g.add_edge(a, b, {"w": 1})
    g.add_edge(a, c, {"w": 10})
    g.add_edge(a, d, {"w": 20})

    removed = g.remove_edges_if(lambda target, meta: meta["w"] < 10)
    assert removed == 1
    assert sorted(n.value.start for n in g.get_neighbors(a)) == [500, 700]


def test_remove_edges_if_by_target():
    """The predicate can also inspect the target Key, not just the metadata."""
    pg = _pg()
    g = pg.Grove(5)
    a, b, c = _chain(g, (10, 20), (30, 40), (50, 60))
    g.add_edge(a, b, {"w": 1})
    g.add_edge(a, c, {"w": 2})

    removed = g.remove_edges_if(lambda target, meta: target.value.start == 30)
    assert removed == 1
    assert [n.value.start for n in g.get_neighbors(a)] == [50]


def test_remove_edges_if_no_match_returns_zero():
    pg = _pg()
    g = pg.Grove(3)
    a, b = _chain(g, (10, 20), (30, 40))
    g.add_edge(a, b, {"w": 5})
    assert g.remove_edges_if(lambda t, m: m["w"] > 100) == 0
    assert g.edge_count() == 1


def test_remove_edges_if_predicate_exception_propagates():
    pg = _pg()
    g = pg.Grove(3)
    a, b = _chain(g, (10, 20), (30, 40))
    g.add_edge(a, b, {"w": 5})
    with pytest.raises(KeyError):
        g.remove_edges_if(lambda t, m: m["nonexistent"])


def test_remove_edges_if_void_grove_takes_target_only():
    """Void-edge groves (BedGrove) expose remove_edges_if with a 1-arg predicate."""
    pg = _pg()
    g = pg.BedGrove(5)

    def ins(s, e):
        return g.insert("chr1", pg.GenomicCoordinate(".", s, e),
                        pg.BedEntry("chr1", s, e))

    a, b, c = ins(10, 20), ins(30, 40), ins(50, 60)
    g.add_edge(a, b)
    g.add_edge(a, c)

    removed = g.remove_edges_if(lambda target: target.value.start == 30)
    assert removed == 1
    assert [n.value.start for n in g.get_neighbors(a)] == [50]


# --------------------------------------------------------------------------- #
# Typed groves keep void edges (C++ binary interop) — no labelled-edge API
# --------------------------------------------------------------------------- #

def test_typed_grove_has_unlabelled_edges_only():
    """BedGrove keeps void edge metadata: the basic + cleanup edge API is present,
    but the labelled-edge methods are not."""
    pg = _pg()
    g = pg.BedGrove(3)
    a = g.insert("chr1", pg.GenomicCoordinate(".", 100, 200), pg.BedEntry("chr1", 100, 200))
    b = g.insert("chr1", pg.GenomicCoordinate(".", 300, 400), pg.BedEntry("chr1", 300, 400))
    g.add_edge(a, b)                 # unlabelled edge works
    assert g.remove_all_edges(a) == 1  # cleanup methods work
    # labelled-edge methods are gated out for void-edge groves
    assert not hasattr(g, "get_edges")
    assert not hasattr(g, "get_neighbors_if")
    assert not hasattr(g, "link_with")


def test_labelled_edge_key_keeps_grove_alive():
    """A neighbour Key from get_neighbors_if keeps the Grove alive (UAF guard)."""
    pg = _pg()
    g = pg.Grove(3)
    a, b = _chain(g, (100, 200), (300, 400))
    g.add_edge(a, b, {"w": 99})
    nbr = g.get_neighbors_if(a, lambda m: m["w"] > 0)[0]
    del g, a, b
    gc.collect()
    assert nbr.value.start == 300


if __name__ == "__main__":
    pytest.main([__file__, "-v"])