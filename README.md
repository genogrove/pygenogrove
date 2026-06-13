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

The standard key is a **`GenomicCoordinate`** (stranded, 0-based closed
`[start, end]`), and the standard **`Grove`** stores any JSON-serializable
payload (dict / list / scalar / `None`) per key:

```python
import pygenogrove as pg

grove = pg.Grove()

# Insert stranded coordinates with arbitrary metadata (or no data at all)
grove.insert("chr1", pg.GenomicCoordinate("+", 100, 200), {"gene": "FOO", "score": 5})
grove.insert("chr1", pg.GenomicCoordinate("-", 100, 200), {"gene": "BAR"})
grove.insert("chr1", pg.GenomicCoordinate(".", 300, 400))   # data defaults to None

# Query is strand-aware: a '+' query matches only '+' (and '*' wildcards)
for key in grove.intersect(pg.GenomicCoordinate("+", 150, 160), "chr1"):
    print(key.value, key.data)        # GenomicCoordinate('+', 100, 200) {'gene': 'FOO', 'score': 5}

# '*' matches any strand; '.' is a concrete unstranded value (matches only '.')
len(grove.intersect(pg.GenomicCoordinate("*", 150, 160), "chr1"))   # 2

grove.serialize("out.gg")             # JSON-text payload; a C++ grove<gc, string> can read it
```

The payload round-trips transparently (no `json` import needed), and each key
may carry a **different** shape — no schema is enforced.

**Important — do not mutate an inserted coordinate.** `GenomicCoordinate.start`,
`.end`, and `.strand` are read-only; `set_range()` / `set_strand()` must only be
used on coordinates you have NOT yet inserted (e.g. a query you want to reuse).
Mutating a stored key silently corrupts B+ tree ordering.

### Strand semantics
- `'+'` / `'-'` — forward / reverse strand
- `'.'` — a concrete *unstranded* value (matches only `'.'`)
- `'*'` — wildcard query strand (matches any strand)

So plain unstranded intervals are just `GenomicCoordinate('.', start, end)`.

### Typed BED/GFF groves (for C++ interop)

The schemaless `Grove` is the everyday tool. When you need a guaranteed BED/GFF
structure and full interop with typed C++ `.gg` files, use the typed groves
(`BedGrove` / `GffGrove`, also genomic-coordinate keyed):

```python
g = pg.BedGrove()
g.insert("chr1", bed_entry)           # entry-deriving: strand taken from the BED6 column
```

## API Reference

### GenomicCoordinate

```python
GenomicCoordinate(strand: str, start: int, end: int)
```

A stranded genomic coordinate with closed `[start, end]` (0-based, both
inclusive). `strand` is one of `'+'`, `'-'`, `'.'`, `'*'`. Overlap requires both
coordinate overlap AND strand compatibility (`'*'` matches any).

**Attributes** (read-only): `strand`, `start`, `end`

**Methods**:
- `set_range(start, end)` / `set_strand(strand)`: pre-insertion only (mutating a stored key corrupts B+ tree ordering).
- `GenomicCoordinate.overlaps(a, b)`: static strand-aware overlap check.

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
- `insert(index: str, key: GenomicCoordinate, data=None) -> Key`: Insert a coordinate (with an optional JSON-serializable payload) at the specified index
- `intersect(query: GenomicCoordinate) -> QueryResult`: Find strand-aware overlaps across all indices
- `intersect(query: GenomicCoordinate, index: str) -> QueryResult`: Find strand-aware overlaps in a specific index
- `flanking(query: GenomicCoordinate, index: str) -> FlankingResult`: Find the nearest **non-overlapping** keys on either side of the query (predecessor / successor). Also `flanking(query, index, is_compatible)` filters candidates by a `bool(candidate, query)` predicate (e.g. same strand)

**FlankingResult** (returned by `flanking`):
- `predecessor`: the closest key entirely before the query (a `Key`), or `None`
- `successor`: the closest key entirely after the query (a `Key`), or `None`

Keys overlapping the query are excluded; for nested intervals the predecessor is
the one with the largest end (smallest gap). Compute the gap distance from the
returned key, e.g. `query.start - result.predecessor.value.end - 1` (closed
coordinates). `BedGrove`/`GffGrove` expose `flanking` too (their results' keys
carry `.data`).

**Graph overlay** (directed edges between keys):
- `add_edge(source: Key, target: Key)`: Add a directed edge (raises `ValueError` if a key is `None`)
- `remove_edge(source: Key, target: Key) -> bool`: Remove an edge; `True` if one was removed
- `has_edge(source: Key, target: Key) -> bool`: Test whether an edge exists
- `get_neighbors(source: Key) -> list[Key]`: Keys directly reachable from `source`
- `out_degree(source: Key) -> int`: Number of outgoing edges from `source`
- `edge_count() -> int`: Total number of edges in the overlay
- `vertex_count_with_edges() -> int`: Number of keys with at least one outgoing edge
- `add_external_key(key: GenomicCoordinate, data=None) -> Key`: Add a key outside the index that can still participate in the graph (not returned by `intersect`)

**Labelled edges** — on the universal `Grove`, edges carry a JSON-serializable payload (the typed `BedGrove`/`GffGrove` keep unlabelled edges for binary interop, so these methods are absent there):
- `add_edge(source: Key, target: Key, data)`: Add an edge with a metadata payload. The 2-argument `add_edge` attaches `None`
- `get_edges(source: Key) -> list`: The edge payloads of `source`'s outgoing edges, parallel to `get_neighbors(source)`
- `get_neighbors_if(source: Key, predicate) -> list[Key]`: Target keys whose edge metadata satisfies `predicate(metadata)` — the predicate receives the **decoded** payload (edges added without a payload yield `None`, so guard for it when mixing labelled and unlabelled edges)
- `link_with(keys: list[Key], predicate)`: Label adjacent pairs — `predicate(k1, k2)` returns the edge payload to attach, or `None` to skip

**Edge removal / bulk linking** (on every grove):
- `remove_edges_from(source: Key) -> int` / `remove_edges_to(target: Key) -> int` / `remove_all_edges(key: Key) -> int`: Remove outgoing / incoming / all touching edges; each returns the count removed
- `clear_graph()`: Remove all edges (keys are left intact); `graph_empty() -> bool`
- `link_if(keys: list[Key], predicate)`: Add an unlabelled edge between each adjacent pair `(keys[i], keys[i+1])` for which `predicate(k1, k2)` returns `True` (typically over the keys returned by a bulk insert)

```python
import pygenogrove as pg

g = pg.Grove()
a = g.insert("chr1", pg.GenomicCoordinate("+", 100, 200))
b = g.insert("chr1", pg.GenomicCoordinate("+", 300, 400))
g.add_edge(a, b, {"type": "exon->transcript", "weight": 7})
g.get_edges(a)                                    # [{"type": ..., "weight": 7}]
g.get_neighbors_if(a, lambda m: m["weight"] > 5)  # [b]
```

**Serialization** (zlib-compressed `.gg` binary):
- `serialize(path: str)`: Write the grove (coordinates + payloads + graph overlay) to `path`
- `deserialize(path: str) -> Grove` *(static)*: Load a grove written by `serialize`

**Removal / storage**:
- `remove_key(index: str, key: Key) -> bool`: Remove a key (and its graph edges); `True` if found. `None`/unknown index → `False`
- `compact()`: Reclaim dead slots left by `remove_key()`. ⚠️ Invalidates every previously-returned indexed `Key` — re-discover via a fresh query afterward
- `vertex_count()` / `external_vertex_count()` / `key_storage_size()`: counts (indexed + external; external-only; total storage slots incl. dead)

### Key

Wrapper object for a coordinate stored in the grove. Returned by insert operations
and yielded by query results. Keeps its `Grove` alive.

**Attributes**:
- `value`: the `GenomicCoordinate` (returned by copy — mutating it cannot corrupt ordering)
- `data`: the payload. On the universal `Grove` this is the JSON value you stored
  (dict / list / scalar / `None`), returned as a freshly decoded copy each access.
  On the typed `BedKey`/`GffKey` it is a live, mutable reference to the record.

### QueryResult

Result object containing matching intervals from a query.

**Attributes**:
- `query`: The query interval used for the search
- `keys`: List of matching keys

**Methods**:
- `__len__()`: Number of results
- `__iter__()`: Iterate over matching keys

### BedGrove (typed BED grove)

`BedGrove` (`grove<genomic_coordinate, bed_entry>`) is the **typed** alternative
to the schemaless `Grove`: instead of a JSON payload, each key carries a
structured `BedEntry`. Use it when you want a guaranteed BED schema and full
interop with typed C++ `.gg` files (prebuilt BED groves load/save with their
records intact, and the GTF-style helpers are available on `GffGrove`).

```python
import pygenogrove as pg

g = pg.BedGrove(100)

# insert(index, coord, data) — the GenomicCoordinate is the key, BedEntry is the payload
entry = pg.BedEntry("chr1", 1000, 2000)   # BED-native coords (0-based, half-open)
entry.name = "BRCA1"
entry.score = 900
entry.strand = "+"
key = g.insert("chr1", pg.GenomicCoordinate(".", 1000, 1999), entry)

# the returned key exposes both the interval value and the BED payload
print(key.value.start, key.data.name)     # 1000 BRCA1

for hit in g.intersect(pg.GenomicCoordinate(".", 1500, 1600), "chr1"):
    print(hit.data.name, hit.data.score)

# serialize/deserialize preserves the BedEntry data
g.serialize("genes.gg")
reloaded = pg.BedGrove.deserialize("genes.gg")
```

`BedGrove` exposes the same surface as `Grove` (multi-index `insert`/`intersect`,
the graph overlay, and `serialize`/`deserialize`), with these differences:
- `insert(index: str, key: GenomicCoordinate, data: BedEntry) -> BedKey` takes the BED payload.
- `add_external_key(key: GenomicCoordinate, data: BedEntry) -> BedKey` takes the payload too.
- **Entry-deriving inserts** (no hand-conversion of coordinates):
  - `insert(index, entry) -> BedKey` — a 2-argument overload: pass a bare
    `BedEntry` and the GenomicCoordinate key is derived from its native coordinates
    (BED's half-open `[s, e)` → closed `[s, e-1]`; GFF's 1-based `[s, e]` →
    `[s-1, e-1]`). This is the foolproof way to load records from a reader.
  - `insert_bulk(index, entries, presorted=False) -> list[BedKey]` — same idea
    for a whole list of bare entries.
- **Fast-path inserts** (data-carrying groves only):
  - `insert_sorted(index, interval, data) -> BedKey` — single insert on the
    rightmost-append path (skips tree traversal).
  - `insert_bulk(index, items, presorted=False) -> list[BedKey]` — insert many
    explicit `(GenomicCoordinate, BedEntry)` records at once (10–20× faster for large
    datasets; an empty index is built bottom-up in O(n)). `presorted=True`
    assumes the records are already sorted by interval (skips the internal sort).
  - **Precondition:** sorted/bulk inserts require ascending intervals, and when
    appending to a non-empty index every new interval must be greater than all
    existing ones. Violating this corrupts B+ tree ordering — use plain `insert`
    if unsure. (`GffGrove` has all the same methods.)

**BedKey** is like `Key` but adds a `data` attribute:
- `value`: the interval (returned by copy; do not rely on mutating it)
- `data`: the associated `BedEntry` — a **live, mutable** reference (unlike `value`,
  the payload is not part of the B+ tree ordering, so editing it in place is safe)

**BedQueryResult** is the `BedGrove` analog of `QueryResult` (its keys are `BedKey`s).

### BedEntry

A single BED record. Coordinates are BED-native: 0-based, half-open `[start, end)`
(distinct from the closed [start, end] of GenomicCoordinate used as the grove key).

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

### GffGrove (typed GFF/GTF grove)

`GffGrove` (`grove<genomic_coordinate, gff_entry>`) is the same typed grove for
**GFF3/GTF** records — identical surface to `BedGrove`, with a `GffEntry` payload
instead of `BedEntry`:

```python
import pygenogrove as pg

g = pg.GffGrove(100)

entry = pg.GffEntry("chr1", 1000, 2000, "gene")   # GFF-native coords (1-based, inclusive)
entry.source = "ensembl"
entry.strand = "+"
entry.attributes = {"gene_id": "ENSG1", "gene_name": "BRCA1"}
key = g.insert("chr1", pg.GenomicCoordinate(".", 999, 1999), entry)

print(key.data.type, key.data.get_gene_id())      # gene ENSG1

for hit in g.intersect(pg.GenomicCoordinate(".", 1500, 1600), "chr1"):
    print(hit.data.type, dict(hit.data.attributes))

g.serialize("genes.gg")
reloaded = pg.GffGrove.deserialize("genes.gg")
```

`GffKey` mirrors `BedKey` (`value` is a copy, `data` is a live mutable `GffEntry`
reference); `GffQueryResult` is the `GffGrove` analog of `QueryResult`.

### GffEntry

A single GFF3/GTF record. Coordinates are GFF-native: **1-based, both endpoints
inclusive** (distinct from GenomicCoordinate's 0-based closed and `BedEntry`'s 0-based
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

# the common workflow: load a file into a grove. The 2-argument insert derives
# the grove's 0-based closed GenomicCoordinate key from each entry's native coordinates,
# so you don't hand-convert (BED half-open, GFF 1-based) yourself.
g = pg.BedGrove(256)
for e in pg.BedReader("peaks.bed"):
    g.insert(e.chrom, e)

gff = pg.GffGrove(256)
for e in pg.GffReader("genes.gff3"):
    gff.insert(e.seqid, e)

# bulk-load one chromosome at a time (insert_bulk is per-index):
g2 = pg.BedGrove(256)
g2.insert_bulk("chr1", [e for e in pg.BedReader("peaks.bed") if e.chrom == "chr1"])
```

```python
BedReader(path: str, skip_invalid_lines: bool = False)
GffReader(path: str, skip_invalid_lines: bool = False, validate_gtf: bool = False)
```

- A missing/unreadable `path` raises on construction.
- With `skip_invalid_lines=False` (default) a malformed line raises `RuntimeError`
  mid-iteration; with `True` such lines are skipped. The **first** data record is
  validated when the reader is constructed, so a malformed first record raises
  immediately regardless of this flag.
- `GffReader(..., validate_gtf=True)` enforces the mandatory GTF2 attributes
  (`gene_id`, `transcript_id`).
- Both expose `get_error_message()` and `get_current_line()` for diagnostics.
- The readers are **single-pass** — they own an htslib file handle and cannot be
  restarted or iterated twice.

> **Coordinate systems** — `GenomicCoordinate` is 0-based closed `[start, end]`; `BedEntry`
> is 0-based half-open `[start, end)`; `GffEntry` is 1-based inclusive `[start, end]`.
> Convert deliberately when building grove keys, as shown above.

### BamReader (SAM/BAM alignments)

`BamReader` is a single-pass iterator over SAM/BAM files (htslib auto-detects the
format) yielding `SamEntry` records, with filtering options applied during
iteration.

```python
import pygenogrove as pg

for aln in pg.BamReader("reads.bam", min_mapq=30):
    print(aln.qname, aln.chrom, aln.start, aln.end, aln.get_strand())

# load alignments into the universal Grove (sam_entry isn't serializable, so
# there's no typed BamGrove — route through to_coordinate() + to_dict())
g = pg.Grove(256)
for aln in pg.BamReader("reads.bam"):
    if aln.is_mapped():
        g.insert(aln.chrom, aln.to_coordinate(), aln.to_dict())
```

```python
BamReader(path, skip_unmapped=True, skip_secondary=False,
          skip_supplementary=False, skip_qc_fail=False,
          skip_duplicates=False, min_mapq=0)
```

- **`SamEntry`** fields: `qname`, `chrom`, `start`, `end` (0-based half-open),
  `mapq`, `sequence`, `quality`, `cigar` (string form), `flags` (an
  `AlignmentFlags`). Helpers: `get_strand()`, `is_primary()` / `is_mapped()` /
  `is_reverse()` / `is_secondary()` / `is_supplementary()` / `is_duplicate()` /
  `is_paired()` / … , `consumes_reference()`, `has_flag(flag)`,
  `to_coordinate()` (strand-aware key) and `to_dict()` (JSON payload).
- **`SamFlags`** exposes the standard FLAG bit constants; **`AlignmentFlags`**
  (the `.flags` object) has `value()` plus the same `is_*()` predicates.
- CIGAR element detail, mate info, and aux tags are not yet exposed.

### FastaReader (FASTA/FASTQ sequences)

`FastaReader` is a single-pass iterator over FASTA/FASTQ files (auto-detected;
`.gz` accepted) yielding `FastaEntry` records. Sequences are named records, not
intervals, so this reader is standalone (no grove integration).

```python
import pygenogrove as pg

for rec in pg.FastaReader("genome.fa"):
    print(rec.name, rec.comment, len(rec), rec.is_fastq())
```

```python
FastaReader(path, skip_empty_sequences=False)
```

- **`FastaEntry`** fields: `name`, `comment`, `sequence`, `quality`
  (`Optional[str]`, FASTQ only); `is_fastq()`, `len(entry)`.

### FastaIndex (random-access FASTA)

`FastaIndex` provides random-access region fetches over a FASTA file, backed by an
`.fai` index (built on first open — the directory must be writable then).

```python
import pygenogrove as pg

fa = pg.FastaIndex("genome.fa")
fa.fetch("chr1", 1000, 2000)   # bases of the 0-based half-open region [1000, 2000)
fa.fetch("chrM")               # the whole sequence
fa.sequence_length("chr1")     # length in bases
list(fa.names()), "chr1" in fa, len(fa)

# fetch a feature's bases: GenomicCoordinate is closed, fetch is half-open
gc = pg.GenomicCoordinate("+", 4, 7)
fa.fetch("chr1", gc.start, gc.end + 1)
```

- Methods: `fetch(name, start, end)` / `fetch(name)`, `sequence_count()`,
  `sequence_name(i)`, `sequence_length(name)`, `has_sequence(name)`, plus the
  Pythonic `len()` / `in` / `names()`. Unknown name / invalid region raise
  `IndexError`.

### FiletypeDetector (format detection)

`FiletypeDetector` infers a file's format and compression from its extension
(compression extension stripped first) and magic bytes.

```python
import pygenogrove as pg

ftype, comp = pg.FiletypeDetector().detect_filetype("peaks.bed.gz")
# (Filetype.BED, CompressionType.GZIP)
```

- `Filetype`: `BED` / `BEDGRAPH` / `GFF` / `GTF` / `VCF` / `SAM` / `BAM` /
  `FASTA` / `FASTQ` / `GG` / `UNKNOWN`.
- `CompressionType`: `NONE` / `GZIP` / `BZIP2` / `XZ` / `ZSTD` / `LZ4` /
  `UNKNOWN`.

### StringRegistry

A process-wide singleton that interns strings into small, stable integer ids
(deduplicated) — handy for chromosome names, sources, gene ids, etc.

```python
import pygenogrove as pg
r = pg.StringRegistry.instance()
a = r.intern("chr1")     # 0
r.intern("chr1")         # 0  (deduplicated)
r.get(a)                 # "chr1"
r.find("chr2")           # None
r.serialize("names.gg")  # also: StringRegistry.deserialize(path), reset(), null_id
```

## Current Status

Currently exposed features:

- **Strand-aware coordinates** — `GenomicCoordinate` is the standard key (`'+'` / `'-'` / `'.'` / `'*'`); overlap and flanking are strand-aware
- **Universal `Grove`** (`grove<genomic_coordinate, json>`) storing arbitrary JSON payloads (dict / list / scalar / `None`), or no payload at all
- Insert / query, multi-index support (per chromosome)
- Graph overlay (directed edges, external keys), including **labelled edges** on the universal `Grove` — `add_edge(s, t, data)` / `get_edges` / `get_neighbors_if` / `link_with` — and edge cleanup / bulk linking on every grove (`remove_edges_from`/`to`, `remove_all_edges`, `clear_graph`, `graph_empty`, `link_if`)
- Key removal + storage compaction: `remove_key()`, `compact()`, `vertex_count()` / `external_vertex_count()` / `key_storage_size()`
- Serialization / deserialization to compressed `.gg` files (an edgeless JSON Grove `.gg` is readable by a C++ `grove<genomic_coordinate, std::string>`; with labelled edges, `grove<genomic_coordinate, std::string, std::string>`)
- Nearest-neighbour queries: `flanking()` (predecessor / successor), incl. a predicate-filtered overload (e.g. same-strand neighbours)
- **Typed** data groves for C++ interop: `BedGrove` (`grove<genomic_coordinate, bed_entry>`) and `GffGrove` (`grove<genomic_coordinate, gff_entry>`), with the `BedEntry` / `GffEntry` value types
- File readers: `BedReader`, `GffReader`, `BamReader` (SAM/BAM), `FastaReader` (FASTA/FASTQ), plus `FastaIndex` (random-access) and `FiletypeDetector` (format detection)
- Fast-path inserts on the typed groves: `insert_sorted` / `insert_bulk`, plus entry-deriving `insert(index, entry)` / `insert_bulk(index, entries)` that derive a **stranded** key from a BED/GFF record's native coordinates
- `StringRegistry` — string interning singleton

**Not yet exposed** (tracked in [#1](https://github.com/genogrove/pygenogrove/issues/1)):
- Other key types — `numeric`, `kmer` ([#7](https://github.com/genogrove/pygenogrove/issues/7))
- `remove_edges_if` ([#33](https://github.com/genogrove/pygenogrove/issues/33)) and SIF export `grove_to_sif` ([#34](https://github.com/genogrove/pygenogrove/issues/34))
- BAM CIGAR-element detail, mate info, and aux tags

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
