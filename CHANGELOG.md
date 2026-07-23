# Changelog

All notable changes to pygenogrove will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.3] - 2026-07-23

### Added

- **`GroveView.get_edge_list(source)`.** The outgoing edges from `source` as
  `(target Key, metadata)` pairs — the zip of `get_neighbors(source)` and
  `get_edges(source)` — read through a partial (random-access) view, paging in
  each target's block on demand. Mirrors the mutable `Grove.get_edge_list`;
  edges added without a payload yield `None` metadata, and a `None` source raises.
  Gated out on the void-edge `BedGroveView` / `GffGroveView`
  ([#1](https://github.com/genogrove/pygenogrove/issues/1),
  [#67](https://github.com/genogrove/pygenogrove/pull/67)).
- **`GroveView.get_order()` and `GroveView.get_index_names()`.** Read the B+ tree
  order and the list of index (chromosome) names straight from the directory
  loaded by `open()`, without paging in any extra blocks — `get_order` mirrors
  `Grove.get_order()`, and `get_index_names` lets a caller discover what
  `intersect` / `flanking` can run against
  ([#1](https://github.com/genogrove/pygenogrove/issues/1),
  [#67](https://github.com/genogrove/pygenogrove/pull/67)).

### Changed

- Bumped the bundled genogrove to
  [v0.25.5](https://github.com/genogrove/genogrove/releases/tag/v0.25.5) (from
  v0.25.3): adds `grove_view::get_edge_list` (parity with `graph_overlay`,
  genogrove [#505](https://github.com/genogrove/genogrove/pull/505)) and
  `grove_view::get_order` / `get_index_names` (genogrove
  [#510](https://github.com/genogrove/genogrove/pull/510)), plus hardening of the
  serialized-`.gg` reader against file-controlled allocation counts and bulk
  builds over `INT_MAX` keys (genogrove
  [#511](https://github.com/genogrove/genogrove/pull/511),
  [#512](https://github.com/genogrove/genogrove/pull/512)).

## [0.7.2] - 2026-07-16

### Added

- **`GroveView.flanking(query, index[, is_compatible])`.** Nearest
  non-overlapping predecessor/successor queries through a partial (random-access)
  view, paging in only the descent-path blocks — the same `FlankingResult` the
  eager `Grove.flanking()` returns, including the predicate-filtered overload for
  same-strand neighbours. This completes the view's read surface (it now mirrors
  every `Grove` query method) ([#61](https://github.com/genogrove/pygenogrove/issues/61)).

### Changed

- Bumped the bundled genogrove to
  [v0.25.3](https://github.com/genogrove/genogrove/releases/tag/v0.25.3) (from
  v0.25.2): adds `grove_view::flanking` (genogrove
  [#483](https://github.com/genogrove/genogrove/pull/483)) and makes the BED/GFF
  readers honor `skip_invalid_lines` for the first record (genogrove
  [#497](https://github.com/genogrove/genogrove/pull/497)) — a malformed leading
  line is now skipped during iteration instead of throwing at construction when
  the option is set ([#62](https://github.com/genogrove/pygenogrove/pull/62)).

## [0.7.1] - 2026-07-13

### Added

- **`GroveView.get_edges(source)` and `get_neighbors_if(source, predicate)`.**
  Read edge metadata payloads through a partial (random-access) view without a
  full `deserialize()`, mirroring the same methods on the mutable `Grove`:
  `get_edges` returns the outgoing edges' payloads in `get_neighbors` order, and
  `get_neighbors_if` returns the target Keys whose decoded payload satisfies a
  Python predicate (paging in each surviving target's block on demand). Both are
  gated out on the void-edge `BedGroveView` / `GffGroveView`
  ([#1](https://github.com/genogrove/pygenogrove/issues/1),
  [#59](https://github.com/genogrove/pygenogrove/pull/59)).

### Changed

- Bumped the bundled genogrove to
  [v0.25.2](https://github.com/genogrove/genogrove/releases/tag/v0.25.2) (from
  v0.25.1): exposes edge-payload reads (`get_edges` / `get_neighbors_if`) on
  `grove_view` (genogrove [#480](https://github.com/genogrove/genogrove/pull/480))
  ([#59](https://github.com/genogrove/pygenogrove/pull/59)).

## [0.7.0] - 2026-07-10

> **⚠️ BREAKING** — the `.gg` serialization format changed to block-structured
> **format 0.2**. Files written by pygenogrove ≤ 0.6.x cannot be read by 0.7.0,
> and 0.7.0 files cannot be read by older versions. Re-serialize existing indexes
> (load with a 0.6.x install, `serialize()` again under 0.7.0). The new format is
> what makes `GroveView`'s partial reads possible.

### Added

- **`GroveView` — partial (random-access) reader over a serialized `.gg`.** A
  read-only view that opens a file written by `Grove.serialize()` and pages in
  only the blocks a query touches, instead of loading the whole grove like
  `Grove.deserialize()` — query a large on-disk index without materializing it.
  Exposes `GroveView.open(path, data_offset=0)`, `intersect(query[, index])`,
  `get_neighbors(key)`, and the `blocks_loaded()` / `block_count()` partial-load
  counters. Bound for every grove flavour (`GroveView` / `NumericGroveView` /
  `KmerGroveView` / `BedGroveView` / `GffGroveView`) via one
  `bind_grove_view<KeyT, DataT, EdgeT>` template mirroring `bind_grove`, reusing
  each grove's `Key` / `QueryResult`. Queries hold the GIL (the view mutates its
  block cache and is not thread-safe — use one view per thread)
  ([#56](https://github.com/genogrove/pygenogrove/issues/56),
  [#57](https://github.com/genogrove/pygenogrove/pull/57)).

### Changed

- Bumped the bundled genogrove to
  [v0.25.1](https://github.com/genogrove/genogrove/releases/tag/v0.25.1) (from
  v0.24.8, skipping v0.25.0): block-structured `.gg` serialization format 0.2
  (which `GroveView` reads), cross-type BED↔GFF intersect in the CLI, and
  `gff_reader` now rejecting a GFF `start < 1`
  ([#57](https://github.com/genogrove/pygenogrove/pull/57)).
- Switched the deprecated `scikit-build-core` config key
  `cmake.minimum-version` to `cmake.version = ">=3.15"` and raised the build
  backend floor to `scikit-build-core>=0.8`; the old key errors on modern
  scikit-build-core, which would break `pip install` from source and CI wheels
  ([#57](https://github.com/genogrove/pygenogrove/pull/57)).

## [0.6.3] - 2026-06-30

### Added

- **Region-based random access for `BedReader`, `GffReader`, and `VcfReader`.**
  Each reader constructor now takes an optional `region` string (e.g.
  `"chr1"`, `"chr1:1000-2000"`) in htslib/tabix query coordinates (1-based,
  inclusive). When set, only records overlapping the region are yielded; this
  requires a bgzip-compressed, tabix/CSI-indexed input (or a BCF for
  `VcfReader`). The default empty string streams the whole file as before
  (genogrove [#456](https://github.com/genogrove/genogrove/pull/456),
  [#458](https://github.com/genogrove/genogrove/pull/458);
  [#55](https://github.com/genogrove/pygenogrove/pull/55)).

### Changed

- Bumped the bundled genogrove to
  [v0.24.8](https://github.com/genogrove/genogrove/releases/tag/v0.24.8).

## [0.6.2] - 2026-06-27

### Fixed

- **`FastaIndex(path)` construction is no longer a data race.** htslib's
  `fai_load()` is not thread-safe; opening several `FastaIndex` handles
  concurrently could abort the interpreter (`SIGABRT`). Construction now takes a
  process-wide lock, so concurrent opens are serialized. The GIL is still
  released during the (potentially long) index build, and `fetch()` on separate
  handles remains fully concurrent ([#50](https://github.com/genogrove/pygenogrove/issues/50)).

## [0.6.1] - 2026-06-26

### Added

- **CPython 3.13 support.** Added `cp313-*` to the cibuildwheel matrix, `3.13`
  to the CI test matrix, and the `Python :: 3.13` classifier. No build changes
  needed.

## [0.6.0] - 2026-06-25

### Added

- **Point key types `Numeric` and `Kmer`.** Two non-interval key types whose overlap is **exact equality** (not range intersection), so their groves act as point dictionaries. `Numeric` wraps an integer (ids / timestamps): `Numeric(value)`, read-only `value` (+ `set_value` for pre-insertion reuse), `overlaps(a, b)`, comparisons, `str`/`repr`. `Kmer` is a 2-bit-encoded DNA k-mer (k ≤ 32, A/C/G/T case-insensitive): `Kmer(sequence)` or `Kmer(encoding, k)`, `encoding` / `k` / `len()`, `overlaps(a, b)`, static `is_valid(sequence)` and `max_k` (invalid bases or `k > 32` raise `ValueError`). Each ships its own `NumericGrove` / `KmerGrove` (plus `*Key` / `*QueryResult` / `*FlankingResult`), each `grove<K, json_value, json_value>` with the same universal surface as `Grove` — optional JSON payload, labelled edges, `.gg` serialization. Reuses the generic `bind_grove<KeyT, DataT>` machinery; no new build dependency ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#7](https://github.com/genogrove/pygenogrove/issues/7), [#39](https://github.com/genogrove/pygenogrove/pull/39)).
- **Value-based `__hash__` on `GenomicCoordinate`, `Numeric`, and `Kmer`.** These key value types now hash by value (consistent with `==` — `GenomicCoordinate` on `(strand, start, end)`, `Numeric` on its int, `Kmer` on `(encoding, k)`), so they work as `set` / `dict` keys instead of hashing by identity ([#39](https://github.com/genogrove/pygenogrove/pull/39)).
- **`remove_edges_if(predicate) -> int`** — predicate-filtered edge removal on every grove, completing the graph-overlay cleanup surface. On the universal `Grove` the predicate is `predicate(target: Key, metadata: object) -> bool` (sees both the target Key and the decoded edge metadata); on void-edge `BedGrove`/`GffGrove` it is `predicate(target: Key) -> bool`. Returns the number of edges removed; the GIL is held since the predicate re-enters Python ([#33](https://github.com/genogrove/pygenogrove/issues/33), [#41](https://github.com/genogrove/pygenogrove/pull/41)).
- **`to_sif(path)` — SIF export on every grove.** Writes the grove's B+ tree structure (`nodelink` / `leaflink` lines) and graph-overlay edges (`keylink` lines) as a tab-separated SIF (Simple Interaction Format) text file for visualization (e.g. Cytoscape). A thin pass-through over genogrove's node-less `grove_to_sif(ostream)` (v0.24.7); GIL released (pure C++). Line/index order is not stable across runs ([#34](https://github.com/genogrove/pygenogrove/issues/34), [#45](https://github.com/genogrove/pygenogrove/pull/45)).
- **`VcfReader` — VCF/BCF variant reader.** A single-pass iterator over VCF/BCF files (plain / bgzip / binary BCF; htslib auto-detects) yielding `VcfEntry`, with `parse_info` / `parse_samples` / `skip_filtered` options and `get_header` / `get_sample_names` / `get_contigs` / `get_error_message` / `get_current_line`. `VcfEntry` exposes the record fields, `info` (htslib-typed INFO → `bool` / `list[int]` / `list[float]` / `str`), per-sample `SampleGenotype` (decoded GT + `gt_string` / `is_hom_ref` + FORMAT fields), and `is_snp` / `is_indel` / `passed_filter` predicates. No typed `VcfGrove` (not serializable) — `to_coordinate()` + `to_dict()` load variants into the universal `Grove`. The GIL is released around the htslib read. Rounds out the file-reader set (BED/GFF/BAM/FASTA → +VCF) ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#46](https://github.com/genogrove/pygenogrove/pull/46)).
- **`get_edge_list(source)` — outgoing edges as `(target, metadata)` pairs.** On the metadata-carrying groves (universal `Grove` / `NumericGrove` / `KmerGrove`), returns `list[tuple[Key, object]]` — the zip of `get_neighbors` and `get_edges`; payload-less edges yield `None` metadata and each target Key keeps the `Grove` alive. Gated out on void-edge `BedGrove`/`GffGrove`, like `get_edges` ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#47](https://github.com/genogrove/pygenogrove/pull/47)).

> **⚠️ BREAKING** — `StringRegistry` is renamed to `Registry`. Single-arg string
> interning is unchanged (`r.intern("chr1")` returns an id; `r.get(id)` returns
> the string); callers must update the class name.

### Changed

- **`StringRegistry` → `Registry`, now mapping a string identity to any JSON payload.** The bound registry is `registry<std::string, void, json_value>`: `intern(key, payload)` interns a string key against any JSON-serializable payload (dict / list / scalar / `None`), deduplicating on the key with **first-write-wins** on the payload, and `get(id)` returns the decoded payload. A single-arg `intern(value)` sugar interns a string as its own payload (plain string interning — `get(id)` returns the string), preserving the previous `StringRegistry` behaviour under the new name. `find` / `contains` / `size` / `__len__` / `empty` / `clear` / static `reset` / `null_id` / `serialize(path)` / static `deserialize(path)` are unchanged; serialization round-trips keys and their JSON payloads. Collapses the per-payload-type registry zoo into one class, mirroring how the universal `Grove` uses `json_value` instead of a typed-per-instantiation template ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#36](https://github.com/genogrove/pygenogrove/pull/36)).
- **Release the GIL on I/O-bound bindings.** The four reader `__next__` methods (`BedReader`/`GffReader`/`BamReader`/`FastaReader`), `FastaIndex` construction + both `fetch` overloads, and grove `serialize` / `deserialize` now drop the GIL around their pure-C++ / htslib work (`py::call_guard<py::gil_scoped_release>()`); `insert_bulk` releases it around just the C++ build loop. Other Python threads run during the native I/O — a throughput win for threaded pipelines with no API or single-threaded behaviour change. The predicate-callback methods (`flanking`/`get_neighbors_if`/`link_*`) keep the GIL since their C++ re-enters Python ([#31](https://github.com/genogrove/pygenogrove/issues/31), [#40](https://github.com/genogrove/pygenogrove/pull/40)).
- **Binding hardening: reject `None` keys on graph methods; cache the JSON caster.** Every graph method taking a `Key` handle (`add_edge`, `remove_edge`, `has_edge`, `get_neighbors`, `out_degree`, `remove_edges_from`/`to`, `remove_all_edges`, `get_edges`, `get_neighbors_if`) now raises `TypeError` on a `None` key — previously the read/remove ones silently returned `False`/`0`, masking caller bugs (`add_edge`'s `None`-key error changes `ValueError` → `TypeError` to match). `remove_key` keeps its documented `None` → `False`. Separately, the `json_value` payload caster caches `json.dumps` / `json.loads` once instead of importing on every conversion ([#43](https://github.com/genogrove/pygenogrove/pull/43)).
- **Bump the bundled genogrove to v0.24.7.** `__genogrove_version__` now reports `0.24.7`. Brings the node-less `grove_to_sif(ostream)` overload (enables an upcoming SIF-export binding), a new htslib VCF/BCF reader in the C++ library, and a bulk-insert contract correction ([#44](https://github.com/genogrove/pygenogrove/pull/44)).

### Fixed

- **Dangling `Key` from vector-of-keys returns.** `get_neighbors`, `get_neighbors_if`, `QueryResult.keys`, and both `insert_bulk` overloads returned `std::vector<key_t*>` with `reference_internal`, which pinned only the resulting list to its parent — extracted `Key` elements got no keep-alive, so a key pulled from the list and outliving it (and every other handle to the owning `Grove`) dangled (use-after-free reachable from pure Python). Each `Key` is now individually pinned to its parent (the `Grove`, or the `QueryResult` which keeps its `Grove` alive); the list stays indexable / `len()`-able, so the API is unchanged ([#37](https://github.com/genogrove/pygenogrove/issues/37), [#38](https://github.com/genogrove/pygenogrove/pull/38)).

## [0.5.0] - 2026-06-13

### Added

- **Completed graph overlay — labelled edges + edge cleanup** — the universal `Grove` is now `grove<genomic_coordinate, json_value, json_value>`, so graph edges carry an arbitrary JSON-serializable payload: `add_edge(source, target, data)`, `get_edges(source)`, `get_neighbors_if(source, predicate)` (predicate on the decoded metadata), and `link_with(keys, predicate)` (optional-returning predicate). Edge cleanup / bulk linking is exposed on every grove: `remove_edges_from`, `remove_edges_to`, `remove_all_edges`, `clear_graph`, `graph_empty`, `link_if`. Typed `BedGrove`/`GffGrove` keep void edge metadata (binary `.gg` C++ interop preserved); the universal `Grove`'s `.gg` now stores per-edge JSON metadata (edgeless files unchanged) ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#32](https://github.com/genogrove/pygenogrove/pull/32)).
- **Random-access FASTA index** — `FastaIndex(path)`, a reader over htslib faidx that builds a missing `.fai` on open. `fetch(name, start, end)` returns the bases of a 0-based half-open region and `fetch(name)` the whole sequence; `sequence_count()` / `sequence_name(i)` / `sequence_length(name)` / `has_sequence(name)` expose per-sequence metadata, with a Pythonic `len()` / `in` / `names()` surface. Pairs with `GenomicCoordinate` (closed `[start, end]`) via `fetch(name, gc.start, gc.end + 1)` ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#30](https://github.com/genogrove/pygenogrove/pull/30)).
- **File-type detector** — `FiletypeDetector.detect_filetype(path) -> (Filetype, CompressionType)`, inferring a file's format and compression from its extension (compression extension stripped first) and magic bytes. The `Filetype` (`BED`/`BEDGRAPH`/`GFF`/`GTF`/`VCF`/`SAM`/`BAM`/`FASTA`/`FASTQ`/`GG`/`UNKNOWN`) and `CompressionType` (`NONE`/`GZIP`/`BZIP2`/`XZ`/`ZSTD`/`LZ4`/`UNKNOWN`) enums are bound too ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#29](https://github.com/genogrove/pygenogrove/pull/29)).
- **FASTA/FASTQ sequence reader** — `FastaReader`, a single-pass iterator over FASTA/FASTQ files (auto-detected; `.gz` accepted) yielding `FastaEntry`, with a `skip_empty_sequences` option and `get_error_message()` / `get_current_line()`. `FastaEntry` exposes `name` / `comment` / `sequence` / `quality` (`Optional[str]`, FASTQ only), `is_fastq()`, and `len()`. Standalone (sequences aren't intervals — no grove integration); random-access (`fasta_index`) deferred ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#28](https://github.com/genogrove/pygenogrove/pull/28)).
- **SAM/BAM alignment reader** — `BamReader`, a single-pass iterator over SAM/BAM files (htslib auto-detects the format) yielding `SamEntry`, with filtering options (`skip_unmapped`/`skip_secondary`/`skip_supplementary`/`skip_qc_fail`/`skip_duplicates`, `min_mapq`). `SamEntry` exposes the core fields + `get_strand()` (from the FLAG), `is_primary()`/`is_mapped()`/`consumes_reference()`, and `cigar` (string form); `AlignmentFlags` (the `.flags` object) + `SamFlags` constants are also bound. Since `sam_entry` isn't serializable there is no typed grove — `SamEntry.to_coordinate()` (strand-aware key) and `.to_dict()` (JSON payload) load alignments into the universal `Grove`. CIGAR element detail, mate info, and aux tags are deferred ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#27](https://github.com/genogrove/pygenogrove/pull/27)).

## [0.4.0] - 2026-06-11

> **⚠️ BREAKING** — `genomic_coordinate` is now the standard key type and `Grove`
> stores arbitrary JSON payloads. The `Interval` value type and the
> interval-keyed groves are removed. Code that used `pg.Interval(...)` /
> interval-keyed `Grove`/`BedGrove`/`GffGrove` must migrate to
> `pg.GenomicCoordinate(strand, start, end)` (plain intervals are
> `GenomicCoordinate('.', start, end)`).

### Changed

- **`Grove` is now `grove<genomic_coordinate, json>` — the universal, strand-aware grove.** Keys are `GenomicCoordinate`; the payload is any JSON-serializable Python object (dict / list / scalar / `None`), stored per key with `insert(index, coord, data=None)` and returned transparently via `key.data` (no user-facing `json` import). Each key may carry a different shape — no schema. It serializes to a `.gg` whose payload is JSON text, so a C++ `grove<genomic_coordinate, std::string>` can still read the file ([#1](https://github.com/genogrove/pygenogrove/issues/1)).
- **`BedGrove` / `GffGrove` are now genomic-coordinate keyed** (`grove<genomic_coordinate, bed_entry/gff_entry>`), kept for typed C++ `.gg` interop + the GTF helper accessors. The entry-deriving `insert(index, entry)` now derives a **stranded** coordinate from the BED6/GFF strand column (absent strand → `'.'`).

### Removed

- **The `Interval` value type and all interval-keyed groves**, plus the verbose `GenomicCoordinate*Grove` / `GenomicCoordinateGrove` names — `genomic_coordinate` is the one key type and the names collapse to `Grove`/`BedGrove`/`GffGrove`. Unstranded intervals are `GenomicCoordinate('.', start, end)`.

## [0.3.0] - 2026-06-11

### Added

- **`StringRegistry`** — genogrove's `registry<std::string>` exposed as a process-wide string-interning singleton: `instance()`, `intern(value) -> int` (idempotent/deduplicated), `find(value) -> int | None`, `get(id) -> str` (`IndexError` on invalid id), `contains`, `size`/`__len__`, `empty`, `clear`, static `reset`, `null_id`, and `serialize(path)` / static `deserialize(path)`. Bound via a generic `bind_registry<Key, Tag, Payload>` template (mirrors genogrove's typed-per-instantiation registry, each with its own id space); only the string instantiation is exposed for now ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#22](https://github.com/genogrove/pygenogrove/pull/22)).
- **Key removal, compaction, and vertex/storage counts** — `remove_key(index, key) -> bool` (unlinks a key from the B+ tree and drops its graph edges; `None`/unknown index returns `False`; the slot stays dead until compaction), `compact()` (reclaims dead slots — **invalidates every previously-returned indexed `Key`**; re-discover via a fresh query; external keys are unaffected), and the counts `vertex_count()` / `external_vertex_count()` / `key_storage_size()`. On every grove ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#21](https://github.com/genogrove/pygenogrove/pull/21)).
- **Predicate-filtered `flanking()`** — a `flanking(query, index, is_compatible)` overload that takes a Python `bool(candidate, query)` callable applied at each leaf candidate, so only matching keys are considered as neighbours. Available on every grove; the canonical use is strand-aware nearest neighbours on a `GenomicCoordinateGrove` (e.g. `g.flanking(q, "chr1", lambda c, q: c.strand == q.strand)` for the nearest same-strand key). The GIL is held across the query (the predicate calls back into Python), and predicate exceptions propagate to Python ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#20](https://github.com/genogrove/pygenogrove/pull/20)).

## [0.2.0] - 2026-06-10

### Added

- **Stranded `GenomicCoordinate` key type and `GenomicCoordinateGrove`** — `grove<genomic_coordinate>` exposed as `GenomicCoordinateGrove` / `GenomicCoordinateKey` / `GenomicCoordinateQueryResult` / `GenomicCoordinateFlankingResult`, plus the `GenomicCoordinate` value type (ctor `(strand, start, end)`; `strand`/`start`/`end`; `set_range`/`set_strand`; `overlaps`; comparisons). Overlap is **strand-aware** — equal strands match, `'*'` is a wildcard matching any strand, `'.'` is a concrete unstranded value — and sorting is coordinate-first (start → end → strand, `* < . < + < -`). Dataless grove only; reuses the generic `bind_grove<KeyT, DataT>` template, so insert / strand-aware intersect / flanking / graph overlay / serialization all come for free ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#15](https://github.com/genogrove/pygenogrove/pull/15)).

### Changed

- **PyPI release pipeline** — tagged releases now build manylinux x86_64 + macOS arm64 (Apple Silicon) wheels for CPython 3.9–3.12 and an sdist via cibuildwheel, then publish to PyPI through OIDC Trusted Publishing (no API tokens). htslib is provisioned and bundled into each wheel (built from source on manylinux, Homebrew on macOS), so `pip install pygenogrove` works without a system htslib. **macOS wheels are arm64-only and require macOS 14.0+** (the arm64 Homebrew bottles set the 14.0 floor; genogrove's `std::format` needs libc++ ≥ 13.3 regardless). **Intel Macs install from the sdist** — GitHub's hosted Intel runners are being retired, so x86_64 macOS wheels are not built ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#14](https://github.com/genogrove/pygenogrove/pull/14), [#16](https://github.com/genogrove/pygenogrove/pull/16), [#17](https://github.com/genogrove/pygenogrove/pull/17), [#18](https://github.com/genogrove/pygenogrove/pull/18)).

## [0.1.0] - 2026-06-09

### Added

- **Version introspection** — `pygenogrove.__version__` (single-sourced from `pyproject.toml` at build time, no longer hardcoded in `bindings.cpp`) and `pygenogrove.__genogrove_version__` (the genogrove C++ library version the wheel was built against, read from genogrove's generated `config/version.hpp`). The two version lines follow independent SemVer cadences ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#13](https://github.com/genogrove/pygenogrove/pull/13)).
- **Graph overlay on `grove<interval>`** — directed-edge operations (`add_edge`, `remove_edge`, `has_edge`, `get_neighbors`, `out_degree`, `edge_count`, `vertex_count_with_edges`) and `add_external_key` for graph-only keys outside the B+ tree index ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#2](https://github.com/genogrove/pygenogrove/pull/2)).
- **Grove serialization** — `serialize(path)` and static `deserialize(path) -> Grove` over the zlib-compressed `.gg` binary format, with open/write failures surfaced as Python exceptions ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#2](https://github.com/genogrove/pygenogrove/pull/2)).
- **CI, packaging, and tests** — GitHub Actions matrix build (Linux/macOS × Python 3.9–3.12), scikit-build-core `pyproject.toml`, and pytest suites for the graph overlay and serialization ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#2](https://github.com/genogrove/pygenogrove/pull/2)).
- **Expanded test coverage** — ported the missing cases from the genogrove C++ suites for the bound surface (interval edge cases, node splits, graph edge directions / traversal / pointer stability, serialization round-trips); test count 25 → 55. The suite now mirrors the genogrove layout: `tests/data_type/` (`test_interval`, `test_key`, `test_query_result`) and `tests/structure/` (`test_dataless_grove`, `test_graph_overlay`, `test_serialization`), with snake_case test names matching the corresponding C++ cases ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#3](https://github.com/genogrove/pygenogrove/pull/3)).
- **Data-carrying `BedGrove`** — `grove<interval, bed_entry>` exposed as `BedGrove` / `BedKey` / `BedQueryResult`, so prebuilt `.gg` files carrying BED records can be loaded, queried, and traversed from Python. `insert(index, interval, data)` and `add_external_key(interval, data)` take a `BedEntry` payload; `BedKey.data` is a live mutable reference (the payload is not part of B+ tree ordering) while `BedKey.value` is still returned by copy; `serialize`/`deserialize` round-trip the BED data. Also exposes the `BedEntry` value type (+ `BlockInfo` / `ThickInfo` / `RgbColor`). The binding sources were reorganized to mirror the genogrove tree (`src/data_type/`, `src/io/`, `src/structure/`), with the dataless and data-carrying groves sharing one `bind_grove<DataT>` template (the existing `Grove` surface is unchanged) ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#4](https://github.com/genogrove/pygenogrove/pull/4)).
- **Data-carrying `GffGrove`** — `grove<interval, gff_entry>` exposed as `GffGrove` / `GffKey` / `GffQueryResult` (same surface as `BedGrove`), with the `GffEntry` value type, the `GffFormat` enum, the column-9 `attributes` as a `dict[str, str]`, and the GTF helper accessors (`get_gene_id`, `get_transcript_id`, `get_exon_number`, `get_gene_name`, `get_gene_biotype`, `is_gtf`/`is_gff3`). Reuses the `bind_interval_grove<DataT>` template; tests port the applicable genogrove `gfffile-test.cpp` cases. The data-carrying key value/data/lifetime tests were split into `tests/data_type/test_{bed,gff}_key.py` to match genogrove's `key_type_tests/` layout ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#8](https://github.com/genogrove/pygenogrove/pull/8)).

- **BED/GFF file readers** — `BedReader` and `GffReader`, single-pass Python iterators over BED and GFF3/GTF files yielding `BedEntry` / `GffEntry`. Options are keyword args (`skip_invalid_lines`; `validate_gtf` for GFF); a missing file raises on construction, a malformed line raises `RuntimeError` (unless skipped), and plain/gzip/BGZF (`.gz`) inputs are auto-detected. Both expose `get_error_message()` / `get_current_line()`. Note: the first data record is validated at construction, so a malformed first record raises regardless of `skip_invalid_lines` ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#9](https://github.com/genogrove/pygenogrove/pull/9)).
- **Sorted / bulk insertion** — `insert_sorted(index, interval, data)` (rightmost-append fast path) and `insert_bulk(index, items, presorted=False)` (insert many `(Interval, data)` records at once; 10–20× faster for large datasets, empty index built bottom-up in O(n)) on the data-carrying groves. Plus **entry-deriving overloads** `insert(index, entry)` and `insert_bulk(index, entries)` that compute the `Interval` key from a BED/GFF record's native coordinates (BED 0-based half-open, GFF 1-based), so callers never hand-convert — the conversion lives in one place (`src/io/entry_interval.hpp`). These require associated data, so they're absent on the dataless `Grove` ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#10](https://github.com/genogrove/pygenogrove/pull/10)).
- **Nearest-neighbour queries** — `Grove.flanking(query, index) -> FlankingResult` returns the closest non-overlapping keys on either side of a query (`.predecessor` / `.successor`, each a `Key` or `None`); keys overlapping the query are excluded and, for nested intervals, the predecessor is the one with the largest end (smallest gap). Exposed on all groves (`FlankingResult` / `BedFlankingResult` / `GffFlankingResult`); the result and its keys keep the grove alive ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#11](https://github.com/genogrove/pygenogrove/pull/11)).

### Fixed

- **`intersect()` result lifetime** — the `QueryResult` returned by `Grove.intersect()` (and the keys it yields) now keep the grove alive via `keep_alive<0,1>`. Previously a key or `.data`/`.value` reference materialized from an intersect result could outlive the grove, a latent use-after-free affecting every grove instantiation ([#8](https://github.com/genogrove/pygenogrove/pull/8)).

### Refactored

- **Generalized the grove binding templates over the key type** — `bind_grove<KeyT, DataT>` (and `bind_key` / `bind_query_result` / `bind_flanking_query_result`) replace the interval-hardcoded versions, so additional key types can be added by instantiation rather than near-duplicate headers; the now-generic `interval_key.hpp` / `interval_grove.hpp` were renamed to `key.hpp` / `grove.hpp`. No change to the Python API ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#12](https://github.com/genogrove/pygenogrove/pull/12)).
