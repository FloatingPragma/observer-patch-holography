"""Adversarial controls for the full-family compiler and retained event codec."""
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import struct

import numpy as np
import pytest

import check_control
import pack
import prepare
import tape
import verify

HERE = Path(__file__).resolve().parent
CONTROLS = HERE/"controls"


@pytest.fixture(scope="module")
def baseline():
    packet = verify.load(CONTROLS/"q3_baseline.json")
    return packet, list(check_control.rows_for(packet, CONTROLS))


@pytest.mark.parametrize("variant", ["baseline", "source", "branch", "scratch"])
def test_every_control_event_and_all_induced_pairs(variant):
    result = check_control.check(verify.load(CONTROLS/f"q3_{variant}.json"), CONTROLS)
    assert result["events"] == 22818
    assert result["logical_pairs"] == 6561


@pytest.mark.parametrize("kind", ["mean_value", "nonlocal_export", "stale_writer", "non_support_mean",
                                  "mutable_capture", "reset_payload", "missing_reset", "missing_read",
                                  "duplicate_read", "old_version", "bad_initial", "extra_event"])
def test_semantic_mutations_without_relying_on_hashes(baseline, kind):
    packet, original = baseline
    rows = [list(row) for row in original]
    at = lambda op: next(i for i, row in enumerate(rows) if row[0] == op)
    if kind == "mean_value":
        rows[at(5)][7] += 1
    elif kind == "nonlocal_export":
        rows[at(3)][4] += 1
    elif kind == "stale_writer":
        rows[at(5)][5] = 0
    elif kind == "non_support_mean":
        rows[at(5)][2] += 12
    elif kind == "mutable_capture":
        rows[at(6)][0] = 1
    elif kind == "reset_payload":
        i = at(4)
        rows[i][1] = rows[i-1][1]
        rows[i][5] = rows[i-1][5]
    elif kind in ("missing_reset", "missing_read"):
        removed = at(6)+1 if kind == "missing_reset" else at(8)
        rows.pop(removed)
        # Coherently renumber the surviving custody, not just a broken hash.
        for row in rows:
            for col in (5,6):
                if removed < row[col] < check_control.ABSENT:
                    row[col] -= 1
    elif kind == "duplicate_read":
        rows.insert(at(8), rows[at(8)][:])
    elif kind == "old_version":
        first_commit = at(9)
        i = next(j for j, row in enumerate(rows) if j > first_commit and row[0] == 8)
        other = next(j for j, row in enumerate(rows) if row[0] == 2)
        rows[i][2] = rows[other][3]
        rows[i][6] = other
    elif kind == "bad_initial":
        rows[0][7] += 2
    else:
        rows.append(rows[-1][:])
    assert rows != [list(row) for row in original], "negative control must actually change the tape"
    with pytest.raises(ValueError):
        check_control.check(packet, CONTROLS, rows)


def test_capacity_undercount_is_recomputed(baseline):
    packet, rows = baseline
    changed = deepcopy(packet)
    changed["costs"]["registers"] -= 1
    with pytest.raises(ValueError, match="resource count"):
        check_control.check(changed, CONTROLS, rows)


@pytest.mark.parametrize("mutation", ["axiom_promotion", "universal_no_go", "physical_clock", "closed_owner", "integer_flag",
                                      "missing_exit", "derived_exit", "missing_read_law", "blank_derivation"])
def test_scientific_scope_and_live_m1_transfer_fail_closed(mutation):
    spec = verify.load(HERE/"specification.json")
    if mutation == "closed_owner":
        spec["M1"]["derivation_owners"] = [740,778]
    elif mutation == "integer_flag":
        spec["scope"]["full_q13_q21_routed_execution"] = 1
    elif mutation == "missing_exit":
        del spec["exit"]
    elif mutation == "derived_exit":
        spec["exit"]["selected_route"] = "axiomatic_read_law_derived"
    elif mutation == "missing_read_law":
        del spec["M1"]["read_law"]
    elif mutation == "blank_derivation":
        spec["M1"]["derivation_obligation"] = "  "
    else:
        key = {"axiom_promotion":"axiomatic_read_law_derived", "universal_no_go":"universal_routing_impossibility",
               "physical_clock":"physical_clock_identified"}[mutation]
        spec["scope"][key] = True
    with pytest.raises(ValueError):
        verify.check_specification(spec)


@pytest.mark.parametrize("mutation", ["missing_reads", "event_undercount", "boolean_hops",
                                       "allocation_undercount", "missing_rss", "controller_undercount"])
def test_complete_resource_accounting_fails_closed(mutation):
    packet = verify.load(CONTROLS/"q3_baseline.json")
    verify.check_resources(packet)
    if mutation == "missing_reads":
        del packet["costs"]["register_reads"]
    elif mutation == "event_undercount":
        packet["costs"]["events"] -= 1
    elif mutation == "boolean_hops":
        packet["costs"]["hops"] = True
    elif mutation == "allocation_undercount":
        packet["implementation_resources"]["register_heap_peak_requested_bytes"] -= 1
    elif mutation == "missing_rss":
        del packet["implementation_resources"]["native_peak_resident_bytes"]
    else:
        packet["routing_control_work"]["bfs_directed_edge_examinations"] -= 1
    with pytest.raises(ValueError):
        verify.check_resources(packet)


def test_codec_retains_every_bit_and_block_boundary():
    rng = np.random.default_rng(777)
    rows = rng.integers(0, 2**64-1, size=(tape.ROWS+3,8), dtype=np.uint64)
    rows[:,0] %= 11
    raw = rows.astype("<u8").tobytes()
    encoded = io.BytesIO()
    summary = tape.encode(io.BytesIO(raw), encoded)
    assert summary["events"] == len(rows)
    assert summary["decoded_sha256"] == hashlib.sha256(raw).hexdigest()
    assert b"".join(tape.decode(io.BytesIO(encoded.getvalue()))) == raw


@pytest.mark.parametrize("mutation", ["magic", "truncate", "extra", "oversize", "compression_suffix"])
def test_codec_rejects_noncanonical_or_incomplete_tapes(mutation):
    output = io.BytesIO()
    tape.encode(io.BytesIO(bytes(64)), output)
    data = bytearray(output.getvalue())
    if mutation == "magic":
        data[0] ^= 1
    elif mutation == "truncate":
        data.pop()
    elif mutation == "extra":
        data.extend(b"extra")
    elif mutation == "oversize":
        struct.pack_into("<I",data,8,tape.ROWS+1)
    else:
        length = struct.unpack_from("<I",data,12)[0]
        struct.pack_into("<I",data,12,length+1)
        data.append(0)
    with pytest.raises(ValueError):
        list(tape.decode(io.BytesIO(data)))


def test_delta_storage_is_addition_of_literal_rows():
    packet = verify.load(CONTROLS/"q3_source.json")
    blocks = list(pack.segment_rows(CONTROLS, packet["tape"]["segments"][1]["path"]))
    assert hashlib.sha256(b"".join(blocks)).hexdigest() == verify.load(CONTROLS/"q3_source_phase1.segment.json")["decoded_sha256"]
    base = np.frombuffer(b"".join(pack.segment_rows(CONTROLS,"q3_baseline_phase1.segment.json")),dtype="<u8").reshape(-1,8)
    changed = np.frombuffer(b"".join(blocks),dtype="<u8").reshape(-1,8)
    assert np.array_equal(base[:,:7], changed[:,:7])
    assert np.any(base[:,7] != changed[:,7])


@pytest.mark.parametrize("value", ['{"x":1,"x":2}', '{"x":NaN}', '{"x":1.0}'])
def test_strict_manifest_json(tmp_path, value):
    path = tmp_path/"bad.json"
    path.write_text(value,encoding="utf-8")
    with pytest.raises(ValueError):
        tape.load_json(path)


@pytest.mark.parametrize("name", ["../escape.json", "..\\escape.json", "C:escape.json"])
def test_segment_paths_are_canonical_on_every_platform(name):
    with pytest.raises(ValueError, match="nonlocal"):
        list(pack.segment_rows(CONTROLS,name))


def test_cyclic_segment_dependency_fails_closed(tmp_path):
    (tmp_path/"cycle.json").write_text(json.dumps({"schema":"oph.source_read_routing.segment.v1", "events":1,
        "parts":[], "decoded_sha256":"0"*64, "delta_base":"cycle.json"}),encoding="utf-8")
    with pytest.raises(ValueError,match="cyclic"):
        list(pack.segment_rows(tmp_path,"cycle.json"))


def test_frozen_input_rebuild_is_byte_identical(tmp_path):
    path = tmp_path/"q3.input"
    result = prepare.prepare(3,path)
    assert path.read_bytes() == (CONTROLS/"q3.input").read_bytes()
    assert "\\" not in result["support_path"]


def test_metric_input_mutation_is_not_accepted_after_resealing(tmp_path):
    packet = verify.load(CONTROLS/"q3_baseline.json")
    data = bytearray((CONTROLS/"q3.input").read_bytes())
    # Change a golden address while preserving the entire container schema.
    data[-8] ^= 1
    path = tmp_path/"q3.input"
    path.write_bytes(data)
    packet["inputs"]["input_sha256"] = hashlib.sha256(data).hexdigest()
    with pytest.raises(ValueError,match="golden address"):
        verify.inspect_input(path,packet["inputs"])
