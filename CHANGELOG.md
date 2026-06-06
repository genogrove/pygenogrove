# Changelog

All notable changes to pygenogrove will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Graph overlay on `grove<interval>`** — directed-edge operations (`add_edge`, `remove_edge`, `has_edge`, `get_neighbors`, `out_degree`, `edge_count`, `vertex_count_with_edges`) and `add_external_key` for graph-only keys outside the B+ tree index ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#2](https://github.com/genogrove/pygenogrove/pull/2)).
- **Grove serialization** — `serialize(path)` and static `deserialize(path) -> Grove` over the zlib-compressed `.gg` binary format, with open/write failures surfaced as Python exceptions ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#2](https://github.com/genogrove/pygenogrove/pull/2)).
- **CI, packaging, and tests** — GitHub Actions matrix build (Linux/macOS × Python 3.9–3.12), scikit-build-core `pyproject.toml`, and pytest suites for the graph overlay and serialization ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#2](https://github.com/genogrove/pygenogrove/pull/2)).
- **Expanded test coverage** — ported the missing cases from the genogrove C++ suites for the bound surface (interval edge cases, node splits, graph edge directions / traversal / pointer stability, serialization round-trips); test count 25 → 55. The suite now mirrors the genogrove layout: `tests/data_type/` (`test_interval`, `test_key`, `test_query_result`) and `tests/structure/` (`test_dataless_grove`, `test_graph_overlay`, `test_serialization`), with snake_case test names matching the corresponding C++ cases ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#3](https://github.com/genogrove/pygenogrove/pull/3)).
- **Data-carrying `BedGrove`** — `grove<interval, bed_entry>` exposed as `BedGrove` / `BedKey` / `BedQueryResult`, so prebuilt `.gg` files carrying BED records can be loaded, queried, and traversed from Python. `insert(index, interval, data)` and `add_external_key(interval, data)` take a `BedEntry` payload; `BedKey.data` is a live mutable reference (the payload is not part of B+ tree ordering) while `BedKey.value` is still returned by copy; `serialize`/`deserialize` round-trip the BED data. Also exposes the `BedEntry` value type (+ `BlockInfo` / `ThickInfo` / `RgbColor`). The binding sources were reorganized to mirror the genogrove tree (`src/data_type/`, `src/io/`, `src/structure/`), with the dataless and data-carrying groves sharing one `bind_grove<DataT>` template (the existing `Grove` surface is unchanged) ([#1](https://github.com/genogrove/pygenogrove/issues/1), [#4](https://github.com/genogrove/pygenogrove/pull/4)).
