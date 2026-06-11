"""
Tests for StringRegistry — the bound registry<std::string> singleton.

Ports the behaviourally-observable cases from genogrove
tests/data_type/registry_test.cpp (RegistryTest, the key == payload form). The
registry is a process-wide singleton, so each test resets it first (autouse
fixture) for isolation. Tagged / key->payload / concurrency cases are out of
scope (not bound / not Python-testable).
"""

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


@pytest.fixture(autouse=True)
def _reset_registry():
    pg = pytest.importorskip("pygenogrove")
    pg.StringRegistry.reset()
    yield
    pg.StringRegistry.reset()


def test_singleton_shares_state():
    """instance() returns the one global pool — state is shared."""
    pg = _pg()
    r1 = pg.StringRegistry.instance()
    r2 = pg.StringRegistry.instance()
    i = r1.intern("chr1")
    assert r2.get(i) == "chr1"
    assert r2.find("chr1") == i


def test_basic_intern_and_get():
    pg = _pg()
    r = pg.StringRegistry.instance()
    i = r.intern("chr1")
    assert i == 0
    assert r.get(i) == "chr1"


def test_distinct_interns_allocate_sequential_ids():
    pg = _pg()
    r = pg.StringRegistry.instance()
    assert r.intern("chr1") == 0
    assert r.intern("chr2") == 1
    assert r.intern("chr3") == 2
    assert r.size() == 3


def test_intern_is_idempotent():
    pg = _pg()
    r = pg.StringRegistry.instance()
    a = r.intern("chr1")
    b = r.intern("chr1")
    assert a == b
    assert r.size() == 1


def test_find_returns_id_for_existing():
    pg = _pg()
    r = pg.StringRegistry.instance()
    i = r.intern("chr1")
    assert r.find("chr1") == i


def test_find_returns_none_for_missing():
    pg = _pg()
    r = pg.StringRegistry.instance()
    assert r.find("nope") is None


def test_find_does_not_insert():
    pg = _pg()
    r = pg.StringRegistry.instance()
    r.find("ghost")
    assert r.size() == 0
    assert r.empty()


def test_invalid_id_raises_index_error():
    pg = _pg()
    r = pg.StringRegistry.instance()
    r.intern("chr1")  # id 0 valid
    with pytest.raises(IndexError):
        r.get(5)
    with pytest.raises(IndexError):
        r.get(pg.StringRegistry.null_id)


def test_contains():
    pg = _pg()
    r = pg.StringRegistry.instance()
    i = r.intern("chr1")
    assert r.contains(i) is True
    assert r.contains(i + 1) is False


def test_size_and_empty_and_len():
    pg = _pg()
    r = pg.StringRegistry.instance()
    assert r.empty() is True
    assert len(r) == 0
    r.intern("chr1")
    r.intern("chr2")
    assert r.empty() is False
    assert r.size() == 2
    assert len(r) == 2


def test_clear_and_ids_restart_from_zero():
    pg = _pg()
    r = pg.StringRegistry.instance()
    r.intern("chr1")
    r.intern("chr2")
    r.clear()
    assert r.empty()
    assert r.intern("chrX") == 0  # ids restart after clear


def test_reset_classmethod_clears_singleton():
    pg = _pg()
    pg.StringRegistry.instance().intern("chr1")
    pg.StringRegistry.reset()
    assert pg.StringRegistry.instance().empty()


def test_null_id_constant():
    pg = _pg()
    assert pg.StringRegistry.null_id == 2**32 - 1


def test_serialize_deserialize_roundtrip(tmp_path):
    pg = _pg()
    r = pg.StringRegistry.instance()
    for s in ("chr1", "chr2", "chrX"):
        r.intern(s)
    path = str(tmp_path / "reg.gg")
    r.serialize(path)

    r.clear()
    assert r.empty()

    pg.StringRegistry.deserialize(path)
    assert r.size() == 3
    assert r.get(0) == "chr1"
    assert r.get(1) == "chr2"
    assert r.get(2) == "chrX"
    assert r.find("chr2") == 1


def test_serialize_deserialize_empty(tmp_path):
    pg = _pg()
    r = pg.StringRegistry.instance()
    path = str(tmp_path / "empty.gg")
    r.serialize(path)
    r.intern("stale")
    pg.StringRegistry.deserialize(path)
    assert r.empty()


def test_deserialize_replaces_existing_data(tmp_path):
    pg = _pg()
    r = pg.StringRegistry.instance()
    r.intern("a")
    r.intern("b")
    path = str(tmp_path / "ab.gg")
    r.serialize(path)

    r.clear()
    r.intern("z")  # different data
    assert r.get(0) == "z"

    pg.StringRegistry.deserialize(path)  # replaces, not merges
    assert r.size() == 2
    assert r.get(0) == "a"
    assert r.get(1) == "b"
    assert r.find("z") is None


def test_serialize_open_failure_raises():
    pg = _pg()
    r = pg.StringRegistry.instance()
    r.intern("chr1")
    with pytest.raises(RuntimeError):
        r.serialize("/nonexistent_dir_xyz/reg.gg")