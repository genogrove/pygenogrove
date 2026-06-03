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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])