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