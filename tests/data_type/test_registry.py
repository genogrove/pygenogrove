"""
Tests for Registry — registry<std::string, void, json_value>.

The universal interning singleton: intern(key, payload) maps a string identity
(gene_id / chrom / transcript_id) to any JSON-serializable Python object,
deduplicating on the key. A single-arg intern(value) sugar interns a string as
its own payload, so get(id) returns the string back (plain string interning).

Ports the behaviourally-observable cases from genogrove
tests/data_type/registry_test.cpp (both the key == payload and key -> payload
forms), plus the JSON round-trip specific to this binding. Registry is a
process-wide singleton, so each test resets it first (autouse fixture) for
isolation. Tagged / concurrency cases are out of scope (not bound / not
Python-testable).
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


@pytest.fixture(autouse=True)
def _reset_registry():
    pg = pytest.importorskip("pygenogrove")
    pg.Registry.reset()
    yield
    pg.Registry.reset()


# --- plain string interning (single-arg sugar) ---------------------------------


def test_singleton_shares_state():
    """instance() returns the one global pool — state is shared."""
    pg = _pg()
    r1 = pg.Registry.instance()
    r2 = pg.Registry.instance()
    i = r1.intern("chr1")
    assert r2.get(i) == "chr1"
    assert r2.find("chr1") == i


def test_single_arg_interns_string_as_its_own_payload():
    pg = _pg()
    r = pg.Registry.instance()
    i = r.intern("chr1")
    assert i == 0
    assert r.get(i) == "chr1"  # string is recoverable from the id


def test_distinct_interns_allocate_sequential_ids():
    pg = _pg()
    r = pg.Registry.instance()
    assert r.intern("chr1") == 0
    assert r.intern("chr2") == 1
    assert r.intern("chr3") == 2
    assert r.size() == 3


def test_intern_is_idempotent():
    pg = _pg()
    r = pg.Registry.instance()
    a = r.intern("chr1")
    b = r.intern("chr1")
    assert a == b
    assert r.size() == 1


# --- key -> JSON payload interning (two-arg form) ------------------------------


def test_intern_returns_id_and_get_returns_payload():
    pg = _pg()
    r = pg.Registry.instance()
    i = r.intern("ENSG001", {"name": "BRCA2", "biotype": "protein_coding"})
    assert i == 0
    assert r.get(i) == {"name": "BRCA2", "biotype": "protein_coding"}


def test_dedup_is_on_the_key_not_the_payload():
    pg = _pg()
    r = pg.Registry.instance()
    a = r.intern("ENSG001", {"name": "BRCA2"})
    b = r.intern("ENSG001", {"name": "ignored"})
    c = r.intern("ENSG002", {"name": "TP53"})
    assert a == b == 0
    assert c == 1
    assert r.size() == 2


def test_first_write_wins_on_payload():
    """Re-interning an existing key keeps the original payload, drops the new one."""
    pg = _pg()
    r = pg.Registry.instance()
    i = r.intern("ENSG001", {"name": "BRCA2", "biotype": "protein_coding"})
    j = r.intern("ENSG001", {"name": "placeholder", "biotype": ""})
    assert i == j
    assert r.get(i) == {"name": "BRCA2", "biotype": "protein_coding"}


def test_heterogeneous_payload_shapes():
    """Payload is arbitrary JSON — dict, list, scalar, None all round-trip."""
    pg = _pg()
    r = pg.Registry.instance()
    assert r.get(r.intern("k_dict", {"a": 1})) == {"a": 1}
    assert r.get(r.intern("k_list", [1, 2, 3])) == [1, 2, 3]
    assert r.get(r.intern("k_str", "label")) == "label"
    assert r.get(r.intern("k_int", 42)) == 42
    assert r.get(r.intern("k_none", None)) is None


def test_single_arg_and_two_arg_share_one_id_space():
    """Both intern forms dedup on the same key in the same pool."""
    pg = _pg()
    r = pg.Registry.instance()
    i = r.intern("ENSG001")  # string-as-payload
    j = r.intern("ENSG001", {"name": "late"})  # first-write-wins -> payload unchanged
    assert i == j == 0
    assert r.get(i) == "ENSG001"
    assert r.size() == 1


# --- lookup / bounds / state ---------------------------------------------------


def test_find_returns_id_for_existing():
    pg = _pg()
    r = pg.Registry.instance()
    i = r.intern("ENSG001", {"name": "BRCA2"})
    assert r.find("ENSG001") == i


def test_find_returns_none_for_missing():
    pg = _pg()
    r = pg.Registry.instance()
    assert r.find("nope") is None


def test_find_does_not_insert():
    pg = _pg()
    r = pg.Registry.instance()
    r.find("ghost")
    assert r.size() == 0
    assert r.empty()


def test_invalid_id_raises_index_error():
    pg = _pg()
    r = pg.Registry.instance()
    r.intern("ENSG001", {"name": "BRCA2"})  # id 0 valid
    with pytest.raises(IndexError):
        r.get(5)
    with pytest.raises(IndexError):
        r.get(pg.Registry.null_id)


def test_contains():
    pg = _pg()
    r = pg.Registry.instance()
    i = r.intern("ENSG001", {"name": "BRCA2"})
    assert r.contains(i) is True
    assert r.contains(i + 1) is False


def test_size_and_empty_and_len():
    pg = _pg()
    r = pg.Registry.instance()
    assert r.empty() is True
    assert len(r) == 0
    r.intern("ENSG001", {"name": "BRCA2"})
    r.intern("ENSG002", {"name": "TP53"})
    assert r.empty() is False
    assert r.size() == 2
    assert len(r) == 2


def test_clear_and_ids_restart_from_zero():
    pg = _pg()
    r = pg.Registry.instance()
    r.intern("ENSG001", {"name": "BRCA2"})
    r.intern("ENSG002", {"name": "TP53"})
    r.clear()
    assert r.empty()
    assert r.intern("ENSG999", {"name": "X"}) == 0  # ids restart after clear


def test_reset_classmethod_clears_singleton():
    pg = _pg()
    pg.Registry.instance().intern("ENSG001", {"name": "BRCA2"})
    pg.Registry.reset()
    assert pg.Registry.instance().empty()


def test_null_id_constant():
    pg = _pg()
    assert pg.Registry.null_id == 2**32 - 1


# --- serialization -------------------------------------------------------------


def test_serialize_deserialize_roundtrips_json_payloads(tmp_path):
    """The headline: keys AND their JSON payloads survive a binary round-trip."""
    pg = _pg()
    r = pg.Registry.instance()
    r.intern("ENSG001", {"name": "BRCA2", "biotype": "protein_coding"})
    r.intern("ENSG002", [1, 2, 3])
    r.intern("ENSG003", "lincRNA")
    path = str(tmp_path / "genes.gg")
    r.serialize(path)

    r.clear()
    assert r.empty()

    pg.Registry.deserialize(path)
    assert r.size() == 3
    assert r.get(0) == {"name": "BRCA2", "biotype": "protein_coding"}
    assert r.get(1) == [1, 2, 3]
    assert r.get(2) == "lincRNA"
    assert r.find("ENSG002") == 1


def test_serialize_deserialize_empty(tmp_path):
    pg = _pg()
    r = pg.Registry.instance()
    path = str(tmp_path / "empty.gg")
    r.serialize(path)
    r.intern("stale", {"x": 1})
    pg.Registry.deserialize(path)
    assert r.empty()


def test_deserialize_replaces_existing_data(tmp_path):
    pg = _pg()
    r = pg.Registry.instance()
    r.intern("ENSG001", {"name": "BRCA2"})
    path = str(tmp_path / "one.gg")
    r.serialize(path)

    r.clear()
    r.intern("OTHER", {"name": "Z"})  # different data
    assert r.get(0) == {"name": "Z"}

    pg.Registry.deserialize(path)  # replaces, not merges
    assert r.size() == 1
    assert r.get(0) == {"name": "BRCA2"}
    assert r.find("OTHER") is None


def test_serialize_open_failure_raises():
    pg = _pg()
    r = pg.Registry.instance()
    r.intern("ENSG001", {"name": "BRCA2"})
    with pytest.raises(RuntimeError):
        r.serialize("/nonexistent_dir_xyz/reg.gg")