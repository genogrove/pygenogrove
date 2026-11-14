"""
insert_sorted / insert_bulk append at the rightmost leaf and never re-check order,
so out-of-order input would silently corrupt the tree. The bindings reject it with
ValueError and leave the grove untouched.
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _gc(pg, start, end=None):
    return pg.GenomicCoordinate(".", start, start + 5 if end is None else end)


def _grove_with(pg, starts, index="chr1"):
    g = pg.Grove()
    for s in starts:
        g.insert(index, _gc(pg, s), None)
    return g


def test_insert_sorted_rejects_value_below_max():
    pg = _pg()
    g = _grove_with(pg, [100, 200])
    with pytest.raises(ValueError, match="chr1"):
        g.insert_sorted("chr1", _gc(pg, 50), None)
    assert g.size() == 2


def test_insert_sorted_accepts_equal_and_greater():
    """Equal values are ties, not violations; ascending appends keep working."""
    pg = _pg()
    g = _grove_with(pg, [100])
    g.insert_sorted("chr1", _gc(pg, 100), None)
    g.insert_sorted("chr1", _gc(pg, 300), None)
    assert g.size() == 3
    assert len(list(g.intersect(_gc(pg, 0, 1000), "chr1"))) == 3


def test_insert_sorted_checks_each_index_separately():
    pg = _pg()
    g = _grove_with(pg, [500], index="chr1")
    g.insert_sorted("chr2", _gc(pg, 1), None)        # chr2 is empty: any value is fine
    g.insert_sorted("chrNew", _gc(pg, 1), None)      # unknown index too
    assert g.size() == 3


@pytest.mark.parametrize("presorted", [False, True])
def test_insert_bulk_rejects_values_below_existing_max(presorted):
    pg = _pg()
    g = _grove_with(pg, [1000])
    items = [(_gc(pg, 10), None), (_gc(pg, 20), None)]
    with pytest.raises(ValueError, match="chr1"):
        g.insert_bulk("chr1", items, presorted=presorted)
    assert g.size() == 1


def test_insert_bulk_presorted_rejects_unsorted_input():
    pg = _pg()
    g = pg.Grove()
    items = [(_gc(pg, 200), None), (_gc(pg, 100), None)]
    with pytest.raises(ValueError, match="presorted"):
        g.insert_bulk("chr1", items, presorted=True)
    assert g.size() == 0


def test_insert_bulk_unsorted_input_still_sorted_when_not_presorted():
    pg = _pg()
    g = pg.Grove()
    items = [(_gc(pg, 200), None), (_gc(pg, 100), None), (_gc(pg, 300), None)]
    keys = g.insert_bulk("chr1", items)
    assert [k.value.start for k in keys] == [100, 200, 300]


def test_insert_bulk_appends_above_existing_max():
    pg = _pg()
    g = _grove_with(pg, [100])
    g.insert_bulk("chr1", [(_gc(pg, 200), None), (_gc(pg, 300), None)], presorted=True)
    assert g.size() == 3
    assert len(list(g.intersect(_gc(pg, 0, 1000), "chr1"))) == 3


def test_insert_bulk_empty_list_is_fine():
    pg = _pg()
    g = _grove_with(pg, [100])
    assert g.insert_bulk("chr1", [], presorted=True) == []


@pytest.mark.parametrize("presorted", [False, True])
def test_entry_overload_rejects_order_violations(presorted):
    """The entry-deriving insert_bulk checks the derived coordinates."""
    pg = _pg()
    g = pg.BedGrove()
    g.insert("chr1", pg.BedEntry("chr1", 1000, 1010))
    with pytest.raises(ValueError, match="chr1"):
        g.insert_bulk("chr1", [pg.BedEntry("chr1", 10, 20)], presorted=presorted)
    assert g.size() == 1


def test_entry_overload_presorted_rejects_unsorted_entries():
    pg = _pg()
    g = pg.BedGrove()
    entries = [pg.BedEntry("chr1", 200, 210), pg.BedEntry("chr1", 100, 110)]
    with pytest.raises(ValueError, match="presorted"):
        g.insert_bulk("chr1", entries, presorted=True)
    assert g.size() == 0


def test_numeric_grove_checked_too():
    pg = _pg()
    n = pg.NumericGrove()
    n.insert("x", pg.Numeric(10))
    with pytest.raises(ValueError):
        n.insert_bulk("x", [(pg.Numeric(1), None)], presorted=True)
    with pytest.raises(ValueError):
        n.insert_sorted("x", pg.Numeric(1), None)
    assert n.size() == 1


def test_rejected_call_leaves_tree_queryable_and_appendable():
    """After a rejection the tree is intact: queries are right and a valid append still works."""
    pg = _pg()
    g = _grove_with(pg, [100, 200])
    with pytest.raises(ValueError):
        g.insert_sorted("chr1", _gc(pg, 50), None)
    with pytest.raises(ValueError):
        g.insert_bulk("chr1", [(_gc(pg, 60), None)])
    assert len(list(g.intersect(_gc(pg, 0, 1000), "chr1"))) == 2
    g.insert_sorted("chr1", _gc(pg, 300), None)
    assert len(list(g.intersect(_gc(pg, 0, 1000), "chr1"))) == 3


def test_presorted_batch_with_ties_is_accepted():
    pg = _pg()
    g = pg.Grove()
    items = [(_gc(pg, 100), None), (_gc(pg, 100), None), (_gc(pg, 200), None)]
    assert len(g.insert_bulk("chr1", items, presorted=True)) == 3


def test_check_follows_the_max_after_remove_key():
    """Removing the rightmost key lowers the index maximum the check compares against."""
    pg = _pg()
    g = pg.Grove()
    keys = [g.insert("chr1", _gc(pg, s), None) for s in (100, 200, 300)]
    assert g.remove_key("chr1", keys[2])
    g.insert_sorted("chr1", _gc(pg, 250), None)          # above 200, below the removed 300
    with pytest.raises(ValueError):
        g.insert_sorted("chr1", _gc(pg, 150), None)

    for k in g.intersect(_gc(pg, 0, 1000), "chr1"):      # empty the index entirely
        g.remove_key("chr1", k)
    g.insert_sorted("chr1", _gc(pg, 1), None)            # gone index: any value is fine
    assert g.size() == 1


def test_check_after_deserialize(tmp_path):
    pg = _pg()
    g = _grove_with(pg, [100, 200, 300])
    path = str(tmp_path / "g.gg")
    g.serialize(path)
    h = pg.Grove.deserialize(path)
    with pytest.raises(ValueError, match="chr1"):
        h.insert_sorted("chr1", _gc(pg, 50), None)
    h.insert_sorted("chr1", _gc(pg, 400), None)
    assert h.size() == 4


def test_gff_entry_overload_checked():
    pg = _pg()
    g = pg.GffGrove()
    g.insert("chr1", pg.GffEntry("chr1", 1000, 1100, "gene"))
    with pytest.raises(ValueError, match="chr1"):
        g.insert_bulk("chr1", [pg.GffEntry("chr1", 10, 20, "gene")])
    g.insert_bulk("chr1", [pg.GffEntry("chr1", 2000, 2100, "gene")], presorted=True)
    assert g.size() == 2
