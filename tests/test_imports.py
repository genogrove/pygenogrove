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