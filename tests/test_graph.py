"""
Tests for the graph overlay (directed edges between keys) exposed on Grove.

Covers add_edge / remove_edge / has_edge / get_neighbors / out_degree /
edge_count / vertex_count_with_edges / add_external_key.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


class TestGraphOverlay:
    def test_add_and_query_edges(self):
        pg = _pg()

        g = pg.Grove(3)
        a = g.insert("chr1", pg.Interval(100, 200))
        b = g.insert("chr1", pg.Interval(300, 400))
        c = g.insert("chr1", pg.Interval(500, 600))

        g.add_edge(a, b)
        g.add_edge(a, c)

        assert g.out_degree(a) == 2
        assert g.out_degree(b) == 0
        assert g.has_edge(a, b)
        assert g.has_edge(a, c)
        assert not g.has_edge(b, a)  # directed

        neighbors = g.get_neighbors(a)
        starts = sorted(n.value.start for n in neighbors)
        assert starts == [300, 500]

        assert g.edge_count() == 2
        assert g.vertex_count_with_edges() == 1

    def test_remove_edge(self):
        pg = _pg()

        g = pg.Grove(3)
        a = g.insert("chr1", pg.Interval(100, 200))
        b = g.insert("chr1", pg.Interval(300, 400))
        g.add_edge(a, b)

        assert g.has_edge(a, b)
        assert g.remove_edge(a, b) is True
        assert g.has_edge(a, b) is False
        # removing a non-existent edge returns False
        assert g.remove_edge(a, b) is False

    def test_edge_from_query_result_key(self):
        """Keys yielded by intersect() are the real stored keys, so edges built
        from them are visible when traversing from an insert()-returned key."""
        pg = _pg()

        g = pg.Grove(3)
        a = g.insert("chr1", pg.Interval(100, 200))
        b = g.insert("chr1", pg.Interval(300, 400))

        # recover `a` via a query and use it as the edge source
        found = list(g.intersect(pg.Interval(100, 200), "chr1"))
        assert len(found) == 1
        g.add_edge(found[0], b)

        assert g.has_edge(a, b)
        assert g.out_degree(a) == 1

    def test_external_key(self):
        pg = _pg()

        g = pg.Grove(3)
        exon = g.insert("chr1", pg.Interval(1000, 1200))
        enhancer = g.add_external_key(pg.Interval(5000, 5500))
        g.add_edge(exon, enhancer)

        # external key is not part of the index
        assert g.size() == 1
        assert g.has_edge(exon, enhancer)

        neighbors = g.get_neighbors(exon)
        assert len(neighbors) == 1
        assert neighbors[0].value.start == 5000

        # external keys are not returned by intersect()
        hits = g.intersect(pg.Interval(5000, 5500))
        assert len(hits) == 0

    def test_add_edge_none_raises(self):
        pg = _pg()

        g = pg.Grove(3)
        a = g.insert("chr1", pg.Interval(100, 200))
        with pytest.raises((ValueError, TypeError)):
            g.add_edge(a, None)

    def test_get_neighbors_none_raises(self):
        pg = _pg()

        g = pg.Grove(3)
        with pytest.raises((ValueError, TypeError)):
            g.get_neighbors(None)

    def test_remove_edge_leaves_other_edges(self):
        """Removing one edge from a multi-edge source keeps the rest."""
        pg = _pg()

        g = pg.Grove(5)
        a = g.insert("chr1", pg.Interval(10, 20))
        b = g.insert("chr1", pg.Interval(25, 35))
        c = g.insert("chr1", pg.Interval(40, 50))
        g.add_edge(a, b)
        g.add_edge(a, c)
        assert g.out_degree(a) == 2

        assert g.remove_edge(a, b) is True
        assert g.out_degree(a) == 1
        neighbors = g.get_neighbors(a)
        assert len(neighbors) == 1
        assert neighbors[0].value.start == 40  # c remains

    def test_no_edges_queries(self):
        """A key with no outgoing edges: out_degree 0, empty neighbor list."""
        pg = _pg()

        g = pg.Grove(3)
        a = g.insert("chr1", pg.Interval(10, 20))
        assert g.out_degree(a) == 0
        assert g.get_neighbors(a) == []
        assert g.edge_count() == 0
        assert g.vertex_count_with_edges() == 0

    def test_branching_path_traversal(self):
        """A source with two targets reports out_degree 2 and both neighbors."""
        pg = _pg()

        g = pg.Grove(5)
        root = g.insert("chr1", pg.Interval(100, 200))
        left = g.insert("chr1", pg.Interval(300, 400))
        right = g.insert("chr1", pg.Interval(500, 600))
        g.add_edge(root, left)
        g.add_edge(root, right)

        assert g.out_degree(root) == 2
        starts = sorted(n.value.start for n in g.get_neighbors(root))
        assert starts == [300, 500]
        assert g.vertex_count_with_edges() == 1

    def test_linear_multi_hop_traversal(self):
        """Walk a 3-hop chain a -> b -> c -> d via get_neighbors."""
        pg = _pg()

        g = pg.Grove(3)
        a = g.insert("chr1", pg.Interval(100, 200))
        b = g.insert("chr1", pg.Interval(300, 400))
        c = g.insert("chr1", pg.Interval(500, 600))
        d = g.insert("chr1", pg.Interval(700, 800))
        g.add_edge(a, b)
        g.add_edge(b, c)
        g.add_edge(c, d)

        node = a
        visited = [node.value.start]
        for _ in range(3):
            nbrs = g.get_neighbors(node)
            assert len(nbrs) == 1
            node = nbrs[0]
            visited.append(node.value.start)
        assert visited == [100, 300, 500, 700]
        assert g.edge_count() == 3
        assert g.vertex_count_with_edges() == 3  # a, b, c each have one out-edge

    def test_edge_directions_with_external_keys(self):
        """Edges in all three directions: indexed->external, external->indexed,
        external->external."""
        pg = _pg()

        g = pg.Grove(5)
        exon = g.insert("chr1", pg.Interval(1000, 1200))
        enh = g.add_external_key(pg.Interval(5000, 5500))
        prom = g.add_external_key(pg.Interval(500, 600))

        g.add_edge(exon, enh)   # indexed -> external
        g.add_edge(prom, exon)  # external -> indexed
        g.add_edge(prom, enh)   # external -> external

        assert g.edge_count() == 3
        assert g.has_edge(exon, enh)
        assert g.has_edge(prom, exon)
        assert g.has_edge(prom, enh)
        assert g.out_degree(prom) == 2
        assert g.out_degree(exon) == 1
        assert g.out_degree(enh) == 0
        assert g.vertex_count_with_edges() == 2  # exon and prom are sources

    def test_mixed_graph_multi_hop(self):
        """Traverse external -> external -> indexed (enhancer -> TF -> gene)."""
        pg = _pg()

        g = pg.Grove(5)
        gene = g.insert("chr1", pg.Interval(1000, 1500))
        tf = g.add_external_key(pg.Interval(0, 0))
        enhancer = g.add_external_key(pg.Interval(50000, 50500))

        g.add_edge(enhancer, tf)
        g.add_edge(tf, gene)

        step1 = g.get_neighbors(enhancer)
        assert len(step1) == 1
        step2 = g.get_neighbors(step1[0])
        assert len(step2) == 1
        assert step2[0].value.start == 1000  # reached the indexed gene

    def test_pointer_stability_across_splits(self):
        """A Key captured before many splits stays valid as an edge endpoint."""
        pg = _pg()

        g = pg.Grove(3)  # small order -> frequent splits
        first = g.insert("chr1", pg.Interval(0, 50))
        for i in range(1, 60):
            g.insert("chr1", pg.Interval(i * 100, i * 100 + 50))
        last = g.insert("chr1", pg.Interval(100_000, 100_050))

        # The pre-split key is still a usable edge endpoint.
        g.add_edge(first, last)
        assert g.has_edge(first, last)
        assert g.out_degree(first) == 1
        assert g.get_neighbors(first)[0].value.start == 100_000

    def test_external_key_pointer_stability(self):
        """Many external keys keep stable identity for edge construction."""
        pg = _pg()

        g = pg.Grove(3)
        ext = [g.add_external_key(pg.Interval(i * 1000, i * 1000 + 500))
               for i in range(100)]
        # Early-created keys must still be valid endpoints after later additions.
        for i in range(99):
            g.add_edge(ext[i], ext[i + 1])
        assert g.edge_count() == 99
        assert g.out_degree(ext[0]) == 1
        assert g.get_neighbors(ext[0])[0].value.start == 1000


if __name__ == "__main__":
    pytest.main([__file__, "-v"])