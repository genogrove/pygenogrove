# pygenogrove

Python bindings for the [genogrove](https://github.com/genogrove/genogrove) C++ library - a specialized B+ tree data structure optimized for genomic interval storage and querying.

**Documentation: https://genogrove.readthedocs.io**

## Installation

```bash
pip install pygenogrove
```

### Building from source

Requires a C++20 compiler, CMake 3.15+ and Python 3.10+.

```bash
git clone --recursive https://github.com/genogrove/pygenogrove.git
cd pygenogrove
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build        # module lands in build/pygenogrove.so (.pyd on Windows)
```

A GCC 13 dev image with htslib and Python headers: `docker/gcc13-dev.Dockerfile`.

## Quick Start

The standard key is a **`GenomicCoordinate`** (stranded, 0-based closed `[start, end]`), and
the standard **`Grove`** stores any JSON-serializable payload per key:

```python
import pygenogrove as pg

grove = pg.Grove()
grove.insert("chr1", pg.GenomicCoordinate("+", 100, 200), {"gene": "FOO", "score": 5})
grove.insert("chr1", pg.GenomicCoordinate("-", 100, 200), {"gene": "BAR"})
grove.insert("chr1", pg.GenomicCoordinate(".", 300, 400))   # data defaults to None

for key in grove.intersect(pg.GenomicCoordinate("+", 150, 160), "chr1"):
    print(key.value, key.data)   # GenomicCoordinate('+', 100, 200) {'gene': 'FOO', 'score': 5}

grove.serialize("out.gg")        # compressed .gg file
```

Queries are strand-aware: `'+'` / `'-'` match only their own strand, `'.'` is a concrete
unstranded value (matches only `'.'`), and `'*'` is a wildcard that matches any strand.

**Do not mutate an inserted coordinate.** `start`, `end` and `strand` are read-only;
`set_range()` / `set_strand()` are for coordinates you have not inserted yet. Mutating a
stored key silently corrupts the B+ tree ordering.

## What's included

- `Grove` (JSON payloads) and typed `BedGrove` / `GffGrove`, plus `Numeric` / `Kmer` point-key groves
- Strand-aware `intersect` and `flanking` queries, bulk and sorted insertion
- Graph overlay (directed, optionally labelled edges), `.gg` serialization and random-access `GroveView`
- Readers for BED, GFF/GTF, SAM/BAM, FASTA/FASTQ and VCF/BCF, plus `FastaIndex`, `FiletypeDetector` and `Registry`

See the [documentation](https://genogrove.readthedocs.io) for the full API and usage guides.
Roadmap and unexposed features are tracked in [#1](https://github.com/genogrove/pygenogrove/issues/1).

## License

GPL-3.0-or-later, inherited from the genogrove C++ library (see `LICENSE`).

## Citation

If you use pygenogrove in your research, please cite the original genogrove library.