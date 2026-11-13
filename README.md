# pygenogrove

Python bindings for the [genogrove](https://github.com/genogrove/genogrove) C++ library - a specialized B+ tree data structure optimized for genomic interval storage and querying.

**Documentation: https://genogrove.readthedocs.io**

## Installation

```bash
pip install pygenogrove
```

To build from source (C++20 compiler, CMake 3.15+, Python 3.10+):

```bash
git clone --recursive https://github.com/genogrove/pygenogrove.git
cd pygenogrove
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
```

A GCC 13 dev image with htslib and Python headers: `docker/gcc13-dev.Dockerfile`.

## Quick Start

```python
import pygenogrove as pg

grove = pg.Grove()
grove.insert("chr1", pg.GenomicCoordinate("+", 100, 200), {"gene": "FOO"})

for key in grove.intersect(pg.GenomicCoordinate("+", 150, 160), "chr1"):
    print(key.value, key.data)
```

See the [documentation](https://genogrove.readthedocs.io) for the full API and usage guides.

## License

GPL-3.0-or-later, inherited from the genogrove C++ library (see `LICENSE`).

## Citation

If you use pygenogrove in your research, please cite the original genogrove library.
