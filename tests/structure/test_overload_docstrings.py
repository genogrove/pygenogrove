"""
Overloaded insert/insert_bulk docstrings must be dedented: pybind11 3.x leaves an
overload's body as written, and Sphinx rejects indented numpy section titles
("Unexpected section title").
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


@pytest.mark.parametrize("cls", ["BedGrove", "GffGrove"])
@pytest.mark.parametrize("method", ["insert", "insert_bulk"])
def test_overloaded_docstring_is_dedented(cls, method):
    """Section titles and their underlines start at column 0."""
    pg = _pg()
    doc = getattr(getattr(pg, cls), method).__doc__
    assert "Overloaded function." in doc
    assert "\nParameters\n----------\n" in doc
    assert "\nReturns\n-------\n" in doc
    assert not [line for line in doc.splitlines() if line.startswith((" ", "\t")) and set(line.strip()) == {"-"}]