"""
pytest configuration for pygenogrove tests.
"""

import sys
from pathlib import Path


def pytest_configure(config):
    """
    Configure pytest to add the build directory to the Python path.

    This allows importing the compiled pygenogrove module from the build directory.
    """
    # Add build directory to sys.path to find the compiled module
    repo_root = Path(__file__).parent.parent
    build_dir = repo_root / "build"

    if build_dir.exists():
        sys.path.insert(0, str(build_dir))
