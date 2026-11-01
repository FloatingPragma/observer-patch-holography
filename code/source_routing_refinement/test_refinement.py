"""Adversarial lifetime, target execution, hierarchy, and receipt controls."""
from __future__ import annotations

import importlib.util
import hashlib
import json
from pathlib import Path
import re
import struct
import sys
from types import ModuleType

import pytest
import yaml

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("m1_refinement_verifier", HERE/"verify.py")
verify = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = verify
spec.loader.exec_module(verify)
storage, hierarchy = verify.storage, verify.hierarchy


@pytest.fixture(scope="module")
def baseline():
    folder = verify.ROOT/"code/source_read_routing/controls"
    packet = json.loads((folder/"q3_baseline.json").read_text(encoding="utf-8"))
    return list(verify.oracle.rows_for(packet, folder)), packet


def test_full_receipt_reproduces():
    actual = verify.build_receipt()
    expected = json.loads((HERE/"receipt.json").read_text(encoding="utf-8"))
    assert actual == expected
    assert sum(c["events"] for c in actual["controls"]) == 91272
    assert all(c["recycled_allocations"] > 0 for c in actual["controls"])


def test_live_zero_slot_cannot_be_reused_for_relay(baseline):
    rows, _ = baseline
    # The protected zero has pending reads when the relay wants its slot.
    with pytest.raises(ValueError, match="live|static allocation alias"):
        storage.replay(rows, 1280, 27, slot_override=lambda eid,op,c,l,p: 23*c+12 if op == storage.CAPTURE else p)


def test_cross_carrier_storage_is_rejected(baseline):
    rows, _ = baseline
    with pytest.raises(ValueError, match="physical owner"):
        storage.replay(rows, 1280, 27, slot_override=lambda eid,op,c,l,p: (p+23) % (1280*23) if op == storage.CAPTURE else p)


def test_stale_writer_is_rejected(baseline):
    rows, _ = baseline
    bad = rows.copy()
    eid = next(i for i,r in enumerate(rows) if r[0] == storage.EXPORT)
    row = list(bad[eid]); row[5] = eid; bad[eid] = tuple(row)
    with pytest.raises(ValueError, match="stale writer"):
        storage.replay(bad, 1280, 27)


def test_retired_version_cannot_be_read_again(baseline):
    rows, _ = baseline
    captures = [r for r in rows if r[0] == storage.CAPTURE]
    first = captures[0]
    # Retiring a captured version is safe only while no later read asks for it.
    # Append a request for that exact writer: the independent last-use pass
    # must reject slot reuse before that read.
    op, a, b, old, owner, wa, wb, val = first
    writer = rows.index(first)
    request = (storage.EXPORT, old, storage.ABSENT, 12*owner, owner,
               writer, storage.ABSENT, val)
    with pytest.raises(ValueError, match="live|retired"):
        storage.replay(rows+[request], 1280, 27)


@pytest.mark.parametrize("field", [1, 5, 7])
def test_independent_physical_replay_rejects_forged_reads(baseline, field):
    rows, _ = baseline
    data, _ = storage.replay(rows, 1280, 27)
    eid = next(i for i,r in enumerate(rows) if r[0] == storage.EXPORT)
    bad = bytearray(data)
    row = list(struct.unpack_from("<8Q", bad, eid*64))
    row[field] = row[field]+1
    struct.pack_into("<8Q", bad, eid*64, *row)
    with pytest.raises(ValueError):
        storage.replay_physical(bytes(bad), 1280)


@pytest.mark.parametrize("fault", ["unknown_opcode", "protected_destination", "extra_operand", "mean_owner"])
def test_independent_physical_replay_rejects_bad_interfaces(baseline, fault):
    rows, _ = baseline
    data, _ = storage.replay(rows, 1280, 27)
    opcode = storage.MEAN if fault == "mean_owner" else storage.EXPORT
    eid = next(i for i,r in enumerate(rows) if r[0] == opcode)
    bad = bytearray(data)
    row = list(struct.unpack_from("<8Q", bad, eid*64))
    if fault == "unknown_opcode":
        row[0] = 99
    elif fault == "protected_destination":
        row[3] = row[4]*23+12
    elif fault == "extra_operand":
        row[2], row[6] = row[1], row[5]
    else:
        row[4] = row[1]//23
    struct.pack_into("<8Q", bad, eid*64, *row)
    with pytest.raises(ValueError):
        storage.replay_physical(bytes(bad), 1280)


def test_valid_store_execution_of_wrong_commit_is_not_a_refinement(baseline):
    rows, _ = baseline
    data, _ = storage.replay(rows, 1280, 27)
    bad = bytearray(data)
    assert rows[-1][0] == storage.COMMIT
    value, = struct.unpack_from("<Q", bad, len(bad)-8)
    struct.pack_into("<Q", bad, len(bad)-8, value+2)
    # Local commit intervention is permitted by the target's store semantics.
    # It must match the particular independently checked source word.
    assert storage.replay_physical(bad, 1280)["events"] == len(rows)
    with pytest.raises(ValueError, match="semantic field drift"):
        storage.check_refinement(rows, bad, 1280)


def test_port_permutation_cannot_change_the_captured_wiring(baseline):
    rows, _ = baseline
    data, _ = storage.replay(rows, 1280, 27)
    changed = []
    for row in struct.iter_unpack("<8Q", data):
        row = list(row)
        for field in (1,2,3):
            address = row[field]
            if address != storage.ABSENT and address%23 in (0,1):
                row[field] = address-address%23+(1-address%23)
        changed.append(struct.pack("<8Q", *row))
    bad = b"".join(changed)
    assert bad != data
    # A global port rename computes the same values and writer graph.
    assert storage.replay_physical(bad, 1280)["events"] == len(rows)
    with pytest.raises(ValueError, match="port identity drift"):
        storage.check_refinement(rows, bad, 1280)


def test_target_store_enforces_the_declared_carrier_set():
    row = (storage.ZERO, storage.ABSENT, storage.ABSENT, 23*1280+12,
           1280, storage.ABSENT, storage.ABSENT, 0)
    with pytest.raises(ValueError, match="local destination"):
        storage.replay_physical(struct.pack("<8Q", *row), 1280)


def test_target_source_comparison_rejects_event_elision(baseline):
    rows, _ = baseline
    data, _ = storage.replay(rows, 1280, 27)
    with pytest.raises(ValueError, match="event count"):
        storage.check_refinement(rows, data[:-64], 1280)


def test_retained_tape_decoder_is_actually_exercised(monkeypatch):
    monkeypatch.setattr(verify.codec, "decode", lambda source: iter([b""]))
    with pytest.raises(ValueError, match="codec round trip"):
        verify.control_receipt()


@pytest.mark.parametrize("text", ['{"events":1,"events":2}', '{"events":1.0}', '{"events":NaN}'])
def test_receipt_json_rejects_ambiguous_numbers_and_keys(tmp_path, text):
    path = tmp_path/"invalid.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        verify.codec.load_json(path)


def test_inherited_imports_are_isolated(monkeypatch):
    sentinels = {name: ModuleType(name) for name in ("verify", "pack", "tape")}
    for name, sentinel in sentinels.items():
        monkeypatch.setitem(sys.modules, name, sentinel)
    isolated = verify.module("m1_import_isolation_control", HERE/"verify.py")
    assert sum(row["events"] for row in isolated.control_receipt()) == 91272
    assert all(sys.modules[name] is sentinel for name, sentinel in sentinels.items())


def test_dedicated_ci_preserves_the_frozen_mandatory_runner():
    projection = verify.codec.load_json(verify.ROOT/"code/invariant_mining/outputs/source_projection.json")
    pin = next(row for row in projection["control_documents"] if row["path"] == "tools/run_mandatory_suite.py")
    runner = (verify.ROOT/pin["path"]).read_bytes()
    assert "sha256:"+hashlib.sha256(runner).hexdigest() == pin["sha256"]
    workflow = yaml.load((verify.ROOT/".github/workflows/source-routing-refinement.yml").read_text(encoding="utf-8"),
                         Loader=yaml.BaseLoader)
    job = workflow["jobs"]["controls"]
    assert set(job["strategy"]["matrix"]["os"]) == {"ubuntu-latest"}
    assert any(s.get("run") == "python -m pytest -q code/source_routing_refinement" for s in job["steps"])
    for trigger in ("push", "pull_request"):
        assert "code/source_routing_refinement/**" in workflow["on"][trigger]["paths"]
        assert "code/source_read_routing/**" in workflow["on"][trigger]["paths"]


def test_axiom_audit_covers_every_refinement_theorem():
    directory = verify.ROOT/"Lean/Geometry"
    declarations = set()
    for stem in ("SourceRoutingStorage", "SourceRoutingBudget", "SourceRoutingHierarchy"):
        source = (directory/f"{stem}.lean").read_text(encoding="utf-8")
        declarations.update(f"OPH.{stem}.{name}" for name in re.findall(r"^theorem ([\w.]+)", source, re.M))
    audit = (directory/"SourceRoutingRefinementAxiomAudit.lean").read_text(encoding="utf-8")
    checked = re.findall(r"^audit_routing_axioms (OPH\.[\w.]+)$", audit, re.M)
    assert len(checked) == len(set(checked)) == 26
    assert set(checked) == declarations


def test_subdivision_mutation_fails():
    fine = list(hierarchy.refine(hierarchy.BASE))
    fine[0] = (999,998,997)
    with pytest.raises(ValueError, match="subdivision"):
        hierarchy.refinement_certificate(hierarchy.BASE, fine)


def test_base_paths_reach_every_carrier():
    graph = hierarchy.graph(hierarchy.BASE)
    distances = [len(p)-1 for root in range(20) for p in hierarchy.paths(graph,root)]
    assert max(distances) == 3
    assert len(distances) == 400


@pytest.mark.parametrize("n", [21,27,125,2197,9261,1000000])
def test_minimal_support_capacity(n):
    level = hierarchy.minimal_level(n)
    C = 20*4**level
    assert n <= C < 4*n
    assert 14*C+9*n <= 23*C


def test_raw_count_bound_uses_all_operations():
    for row in verify.production_bounds():
        assert row["inherited_events"] >= 3*row["logical_events"]
        assert row["derived_working_slot_upper_bound"] < row["inherited_archived_registers"]
