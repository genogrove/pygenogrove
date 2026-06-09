"""
The module reports its own version and the genogrove it was built against.

No genogrove C++ counterpart — this guards the version wiring (``__version__``
single-sourced from pyproject.toml, ``__genogrove_version__`` from genogrove's
generated ``config/version.hpp``). Both follow independent SemVer.
"""

import re

SEMVER_CORE = re.compile(r"^\d+\.\d+\.\d+")


def test_has_version_attrs():
    import pygenogrove
    assert hasattr(pygenogrove, "__version__")
    assert hasattr(pygenogrove, "__genogrove_version__")


def test_versions_are_nonempty_strings():
    import pygenogrove
    assert isinstance(pygenogrove.__version__, str)
    assert isinstance(pygenogrove.__genogrove_version__, str)
    assert pygenogrove.__version__
    assert pygenogrove.__genogrove_version__


def test_genogrove_version_is_semver():
    """The bound genogrove version comes from real MAJOR.MINOR.PATCH macros."""
    import pygenogrove
    assert SEMVER_CORE.match(pygenogrove.__genogrove_version__)