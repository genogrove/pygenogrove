# pygenogrove

Python bindings for the [genogrove](https://github.com/genogrove/genogrove) C++ library - a specialized B+ tree data structure optimized for genomic interval storage and querying.

## Installation

### Building from Source

Requirements:
- C++20 compatible compiler
- CMake 3.15+
- Python 3.8+

```bash
# Clone with submodules
git clone --recursive https://github.com/genogrove/pygenogrove.git
cd pygenogrove

# Build using CMake
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# The built module will be in build/pygenogrove.so (or .pyd on Windows)
```

### Using pip 

```bash
pip install pygenogrove
```

### Using conda/mamba

## Quick Start

```python
import pygenogrove as pg

# Create a grove with order 100 (max 99 keys per node)
grove = pg.Grove(100)

# Create intervals — coordinates are closed [start, end] (both inclusive)
interval1 = pg.Interval(100, 200)
interval2 = pg.Interval(150, 250)
interval3 = pg.Interval(300, 400)

# Insert intervals into different chromosomes
grove.insert("chr1", interval1)
grove.insert("chr1", interval2)
grove.insert("chr2", interval3)

print(f"Total intervals: {len(grove)}")  # Output: Total intervals: 3

# Query for overlapping intervals
query = pg.Interval(175, 225)
results = grove.intersect(query, "chr1")

print(f"Found {len(results)} overlapping intervals")
for key in results:
    interval = key.value
    print(f"  {interval.start}-{interval.end}")
```

## Usage Examples

### Basic Operations

```python
import pygenogrove as pg

# Create a grove (default order is 3; minimum is 3)
grove = pg.Grove()

# Create and insert intervals
interval = pg.Interval(1000, 2000)
key = grove.insert("chr1", interval)

# Access interval properties (read-only)
print(f"Start: {interval.start}")  # Output: Start: 1000
print(f"End: {interval.end}")      # Output: End: 2000
```

**Important — do not mutate an inserted interval.** `Interval.start` and
`Interval.end` are intentionally read-only, and `Interval.set_range(start, end)`
must only be used on intervals you have NOT yet inserted (e.g. a query
interval you want to reuse). Mutating a stored key silently corrupts B+ tree
ordering — overlap queries will start returning wrong answers with no error.

### Querying Intervals

```python
import pygenogrove as pg

grove = pg.Grove(100)

# Insert some intervals
grove.insert("chr1", pg.Interval(100, 200))
grove.insert("chr1", pg.Interval(300, 400))
grove.insert("chr2", pg.Interval(100, 200))

# Query specific chromosome
query = pg.Interval(150, 350)
results = grove.intersect(query, "chr1")

print(f"Found {len(results)} overlaps in chr1")
for key in results:
    print(f"  Interval: {key.value}")

# Query across all chromosomes
all_results = grove.intersect(query)
print(f"Found {len(all_results)} overlaps across all chromosomes")
```

### Overlap Detection

```python
import pygenogrove as pg

# Static method for checking overlap
interval1 = pg.Interval(100, 200)
interval2 = pg.Interval(150, 250)
interval3 = pg.Interval(300, 400)

print(pg.Interval.overlaps(interval1, interval2))  # True (they overlap)
print(pg.Interval.overlaps(interval1, interval3))  # False (no overlap)
```

## API Reference

### Interval

```python
Interval(start: int, end: int)
```

A genomic interval with closed `[start, end]` coordinates (0-based, both inclusive).

**Attributes** (read-only):
- `start`: Start position (inclusive)
- `end`: End position (inclusive)

**Methods**:
- `set_range(start, end)`: Atomically set both endpoints. Only safe on intervals not yet inserted into a Grove (mutating a stored key corrupts B+ tree ordering).
- `Interval.overlaps(a, b)`: Static method to check if two intervals overlap

### Grove

```python
Grove(order: int = 3)
```

A B+ tree container for genomic intervals with multi-index support.

**Parameters**:
- `order`: Maximum branching factor (max keys per node = order - 1). Minimum 3.

**Methods**:
- `len(grove)` / `size()` / `indexed_vertex_count()`: Number of indexed intervals across all indices
- `get_order()`: Get the order (branching factor) of the tree
- `insert(index: str, interval: Interval) -> Key`: Insert an interval at the specified index
- `intersect(query: Interval) -> QueryResult`: Find overlapping intervals across all indices
- `intersect(query: Interval, index: str) -> QueryResult`: Find overlapping intervals in specific index

**Graph overlay** (directed edges between keys):
- `add_edge(source: Key, target: Key)`: Add a directed edge (raises `ValueError` if a key is `None`)
- `remove_edge(source: Key, target: Key) -> bool`: Remove an edge; `True` if one was removed
- `has_edge(source: Key, target: Key) -> bool`: Test whether an edge exists
- `get_neighbors(source: Key) -> list[Key]`: Keys directly reachable from `source`
- `out_degree(source: Key) -> int`: Number of outgoing edges from `source`
- `edge_count() -> int`: Total number of edges in the overlay
- `vertex_count_with_edges() -> int`: Number of keys with at least one outgoing edge
- `add_external_key(interval: Interval) -> Key`: Add a key outside the index that can still participate in the graph (not returned by `intersect`)

**Serialization** (zlib-compressed `.gg` binary):
- `serialize(path: str)`: Write the grove (intervals + graph overlay) to `path`
- `deserialize(path: str) -> Grove` *(static)*: Load a grove written by `serialize`

### Key

Wrapper object for intervals stored in the grove. Returned by insert operations.

**Attributes**:
- `value`: The interval value of this key

### QueryResult

Result object containing matching intervals from a query.

**Attributes**:
- `query`: The query interval used for the search
- `keys`: List of matching keys

**Methods**:
- `__len__()`: Number of results
- `__iter__()`: Iterate over matching keys

### BedGrove (interval grove with BED data)

`BedGrove` is the data-carrying counterpart of `Grove`: each indexed interval
also carries an associated `BedEntry` payload, so prebuilt `.gg` files that
store BED records can be loaded, queried, and traversed from Python.

```python
import pygenogrove as pg

g = pg.BedGrove(100)

# insert(index, interval, data) — the interval is the key, BedEntry is the payload
entry = pg.BedEntry("chr1", 1000, 2000)   # BED-native coords (0-based, half-open)
entry.name = "BRCA1"
entry.score = 900
entry.strand = "+"
key = g.insert("chr1", pg.Interval(1000, 1999), entry)

# the returned key exposes both the interval value and the BED payload
print(key.value.start, key.data.name)     # 1000 BRCA1

for hit in g.intersect(pg.Interval(1500, 1600), "chr1"):
    print(hit.data.name, hit.data.score)

# serialize/deserialize preserves the BedEntry data
g.serialize("genes.gg")
reloaded = pg.BedGrove.deserialize("genes.gg")
```

`BedGrove` exposes the same surface as `Grove` (multi-index `insert`/`intersect`,
the graph overlay, and `serialize`/`deserialize`), with two differences:
- `insert(index: str, interval: Interval, data: BedEntry) -> BedKey` takes the BED payload.
- `add_external_key(interval: Interval, data: BedEntry) -> BedKey` takes the payload too.

**BedKey** is like `Key` but adds a `data` attribute:
- `value`: the interval (returned by copy; do not rely on mutating it)
- `data`: the associated `BedEntry` — a **live, mutable** reference (unlike `value`,
  the payload is not part of the B+ tree ordering, so editing it in place is safe)

**BedQueryResult** is the `BedGrove` analog of `QueryResult` (its keys are `BedKey`s).

### BedEntry

A single BED record. Coordinates are BED-native: 0-based, half-open `[start, end)`
(distinct from the closed `[start, end]` of `Interval` used as the grove key).

```python
BedEntry(chrom: str, start: int, end: int)
```

**Attributes** (read/write):
- `chrom` (str), `start` (int), `end` (int)
- `name`: `Optional[str]` (BED4+)
- `score`: `Optional[int]` (BED5+)
- `strand`: `Optional[str]` — a single character (`'+'`, `'-'`, `'.'`); assigning an
  empty or multi-character string raises `ValueError`, `None` clears it (BED6+)
- `thickness`: `Optional[ThickInfo]` (BED7+)
- `item_rgb`: `Optional[RgbColor]` (BED9+)
- `blocks`: `Optional[BlockInfo]` (BED12)

`ThickInfo(start, end)`, `RgbColor(red, green, blue)` (channels 0–255), and
`BlockInfo(count, sizes, starts)` (with `list[int]` `sizes`/`starts`) are the
supporting value types. List fields are returned/assigned by copy.

### GffGrove (interval grove with GFF/GTF data)

`GffGrove` is the same data-carrying grove for **GFF3/GTF** records — identical
surface to `BedGrove`, with a `GffEntry` payload instead of `BedEntry`:

```python
import pygenogrove as pg

g = pg.GffGrove(100)

entry = pg.GffEntry("chr1", 1000, 2000, "gene")   # GFF-native coords (1-based, inclusive)
entry.source = "ensembl"
entry.strand = "+"
entry.attributes = {"gene_id": "ENSG1", "gene_name": "BRCA1"}
key = g.insert("chr1", pg.Interval(999, 1999), entry)

print(key.data.type, key.data.get_gene_id())      # gene ENSG1

for hit in g.intersect(pg.Interval(1500, 1600), "chr1"):
    print(hit.data.type, dict(hit.data.attributes))

g.serialize("genes.gg")
reloaded = pg.GffGrove.deserialize("genes.gg")
```

`GffKey` mirrors `BedKey` (`value` is a copy, `data` is a live mutable `GffEntry`
reference); `GffQueryResult` is the `GffGrove` analog of `QueryResult`.

### GffEntry

A single GFF3/GTF record. Coordinates are GFF-native: **1-based, both endpoints
inclusive** (distinct from `Interval`'s 0-based closed and `BedEntry`'s 0-based
half-open).

```python
GffEntry(seqid: str, start: int, end: int, type: str)
```

**Attributes** (read/write):
- `seqid` (str), `source` (str), `type` (str), `start` (int), `end` (int)
- `score`: `Optional[float]`
- `strand`: `Optional[str]` — a single character (`'+'`, `'-'`, `'.'`, `'?'`); empty
  or multi-character raises `ValueError`, `None` clears it
- `phase`: `Optional[int]` (CDS phase 0/1/2)
- `attributes`: `dict[str, str]` — the column-9 key/value pairs (returned/assigned by copy)
- `format`: a `GffFormat` enum (`GFF3` / `GTF` / `UNKNOWN`)

**Methods**: `is_gtf()`, `is_gff3()`, `get_attribute(key)`, and the GTF helpers
`get_gene_id()`, `get_transcript_id()`, `get_exon_number()`, `get_gene_name()`,
`get_gene_biotype()` (each returns `None` when the attribute is absent).

### BedReader / GffReader (file iterators)

`BedReader` and `GffReader` are single-pass iterators over BED and GFF3/GTF
files. Iterate them to get `BedEntry` / `GffEntry` records. Plain and
gzip/BGZF-compressed (`.gz`) files are both accepted (auto-detected).

```python
import pygenogrove as pg

# read records one at a time
for entry in pg.BedReader("peaks.bed"):
    print(entry.chrom, entry.start, entry.end, entry.name)

# the common workflow: load a file into a grove (converting to the grove's
# closed [start, end] interval key)
g = pg.BedGrove(256)
for e in pg.BedReader("peaks.bed"):           # BED is 0-based half-open
    g.insert(e.chrom, pg.Interval(e.start, e.end - 1), e)

gff = pg.GffGrove(256)
for e in pg.GffReader("genes.gff3"):          # GFF is 1-based inclusive
    gff.insert(e.seqid, pg.Interval(e.start - 1, e.end - 1), e)
```

```python
BedReader(path: str, skip_invalid_lines: bool = False)
GffReader(path: str, skip_invalid_lines: bool = False, validate_gtf: bool = False)
```

- A missing/unreadable `path` raises on construction.
- With `skip_invalid_lines=False` (default) a malformed line raises `RuntimeError`
  mid-iteration; with `True` such lines are skipped.
- `GffReader(..., validate_gtf=True)` enforces the mandatory GTF2 attributes
  (`gene_id`, `transcript_id`).
- Both expose `get_error_message()` and `get_current_line()` for diagnostics.
- The readers are **single-pass** — they own an htslib file handle and cannot be
  restarted or iterated twice.

> **Coordinate systems** — `Interval` is 0-based closed `[start, end]`; `BedEntry`
> is 0-based half-open `[start, end)`; `GffEntry` is 1-based inclusive `[start, end]`.
> Convert deliberately when building grove keys, as shown above.

## Current Status

This is an early development version. Currently exposed features:

- Basic grove and interval operations
- Insert and query functionality
- Multi-index support (per chromosome)
- Graph overlay (directed edges, external keys)
- Serialization / deserialization to compressed `.gg` files
- Associated data: the `BedEntry` / `GffEntry` value types and the data-carrying
  groves `BedGrove` (`grove<interval, bed_entry>`) and `GffGrove`
  (`grove<interval, gff_entry>`)
- File readers: `BedReader` and `GffReader` (single-pass iterators over BED /
  GFF3 / GTF files, including `.gz`)

**Not yet exposed** (tracked in [#1](https://github.com/genogrove/pygenogrove/issues/1)):
- Genomic coordinates with strand information, and other key types — numeric, kmer
  ([#7](https://github.com/genogrove/pygenogrove/issues/7))
- BAM/SAM and FASTA readers
- Bulk / sorted insertion
- Edge metadata, `get_neighbors_if` / `link_if` (require a metadata-carrying grove)

## Performance Tips

1. **Choose appropriate order**: Higher order (e.g., 100-500) reduces tree height for large datasets
2. **Separate by chromosome**: Use the index parameter to maintain separate trees per chromosome
3. **Query specific indices**: Query specific chromosomes instead of all indices when possible

## License

This project inherits the license from the genogrove C++ library and is therefore licensed under the GPLv3 license. 

## Related Projects

- [genogrove](https://github.com/genogrove/genogrove): The underlying C++ library

## Citation

If you use pygenogrove in your research, please cite the original genogrove library
