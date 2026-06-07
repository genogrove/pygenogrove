"""
Build sanity check: the compiled extension imports and exposes its public types.

This has no genogrove C++ counterpart — it guards the Python packaging/build
itself. It is a hard failure (not a skip) on purpose: if the compiled extension
is not installable/importable, the rest of the suite would otherwise pass by
silently skipping every test, hiding a broken build.
"""


def test_imports():
    import pygenogrove
    assert hasattr(pygenogrove, 'Interval')
    assert hasattr(pygenogrove, 'Grove')
    assert hasattr(pygenogrove, 'Key')
    assert hasattr(pygenogrove, 'QueryResult')


def test_imports_bed_data_grove():
    """The data-carrying grove for bed_entry and its value types are exposed."""
    import pygenogrove
    assert hasattr(pygenogrove, 'BedGrove')
    assert hasattr(pygenogrove, 'BedKey')
    assert hasattr(pygenogrove, 'BedQueryResult')
    assert hasattr(pygenogrove, 'BedEntry')
    assert hasattr(pygenogrove, 'BlockInfo')
    assert hasattr(pygenogrove, 'ThickInfo')
    assert hasattr(pygenogrove, 'RgbColor')


def test_imports_gff_data_grove():
    """The data-carrying grove for gff_entry and its value types are exposed."""
    import pygenogrove
    assert hasattr(pygenogrove, 'GffGrove')
    assert hasattr(pygenogrove, 'GffKey')
    assert hasattr(pygenogrove, 'GffQueryResult')
    assert hasattr(pygenogrove, 'GffEntry')
    assert hasattr(pygenogrove, 'GffFormat')


def test_imports_file_readers():
    """The BED/GFF file readers are exposed."""
    import pygenogrove
    assert hasattr(pygenogrove, 'BedReader')
    assert hasattr(pygenogrove, 'GffReader')