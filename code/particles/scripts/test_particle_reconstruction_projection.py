"""Cheap custody/scope tests; never run the structural or scientific builders."""
from copy import deepcopy
import hashlib
import json

import pytest

import build_postdiction_ledger as ledger
from tools.renew_native_defect_live import RETAINED as NATIVE_ARCHIVE


@pytest.fixture
def pinned_receipts(tmp_path):
    source = tmp_path/"source.py"
    source.write_text("# fixture source\n")
    pins = {"source.py": hashlib.sha256(source.read_bytes()).hexdigest()}
    inputs = {
        "native_defect_decision_live_results.json": (),
        "particle_identification/results.json": ("verdict", "finite_inverse_response_kinetics"),
        "chiral_candidate/record_two_point_results.json": ("scope", "joint_response", "held_out_validation"),
        "chiral_candidate/results.json": ("verdict", "kinematics"),
    }
    for name, fields in inputs.items():
        original = json.loads((ledger.REPO/"model/common_source"/name).read_text())
        packet = {key: deepcopy(original[key]) for key in ("schema", *fields)}
        packet["source_bindings"] = pins
        if name == 'native_defect_decision_live_results.json':
            packet = deepcopy(original)
            for relative in packet['source_sha256']:
                target = tmp_path / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes((ledger.REPO / relative).read_bytes())
            retained = tmp_path / NATIVE_ARCHIVE
            retained.parent.mkdir(parents=True, exist_ok=True)
            retained.write_bytes((ledger.REPO / retained.relative_to(tmp_path)).read_bytes())
        if "finite_inverse_response_kinetics" in packet:
            packet["finite_inverse_response_kinetics"]["source_sha256"] = pins
        path = tmp_path/"model/common_source"/name
        path.parent.mkdir(parents=True, exist_ok=True)
        if "artifact" in original:
            artifact = path.with_suffix(".npz")
            artifact.write_bytes(b"independent custody fixture, no scientific arrays")
            packet["artifact"] = {"filename": artifact.name,
                                  "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
        path.write_text(json.dumps(packet))
    return tmp_path


def test_projection_preserves_separate_objects_and_replay_boundary(pinned_receipts):
    result = ledger._particle_reconstruction_decisions(repo=pinned_receipts)
    rows = {row["object_id"]: row for row in result["rows"]}
    assert len(rows) == 5
    collision = rows["native_total_two_load_collision"]
    assert collision["evidence"]["mean_attempts_to_split"] == "6"
    assert collision["evidence"]["closed_component_sizes"] == [66]
    assert rows["native_active_matching_memory"]["evidence"]["passive_code_escape_probability"] == {"C": "31/108", "E": "55/108"}
    assert rows["supplied_four_mode_CAR_benchmark"]["parent"] != collision["parent"]
    assert "not a fresh scientific replay" in result["verification"]
    assert result['native_retained_evidence']['historical_documentation_pin_claimed_current'] is False
    assert result['native_retained_evidence']['scientific_payload_exactly_preserved'] is True


@pytest.mark.parametrize("damage", ["source", "artifact", "missing_pins", "mass_claim", "record_claim"])
def test_projection_rejects_custody_or_scope_damage(pinned_receipts, damage):
    path = pinned_receipts/"model/common_source/particle_identification/results.json"
    packet = json.loads(path.read_text())
    if damage == "source":
        (pinned_receipts/"source.py").write_text("changed")
    elif damage == "artifact":
        path.with_suffix(".npz").write_bytes(b"changed")
    elif damage == "missing_pins":
        packet["source_bindings"] = {}
    elif damage == "mass_claim":
        packet["verdict"]["physical_particle_identification_complete"] = True
    else:
        packet["finite_inverse_response_kinetics"]["response_inferred_from_observer_records"] = True
    path.write_text(json.dumps(packet))
    with pytest.raises(ValueError):
        ledger._particle_reconstruction_decisions(repo=pinned_receipts)


@pytest.mark.parametrize('damage', ('science', 'type', 'coverage', 'retained'))
def test_live_native_projection_cannot_rewrite_retained_science(pinned_receipts, damage):
    path = pinned_receipts / 'model/common_source/native_defect_decision_live_results.json'
    packet = json.loads(path.read_text())
    if damage == 'science':
        packet['carrier']['states'] = 79
    elif damage == 'type':
        packet['carrier']['states'] = 78.0
    elif damage == 'coverage':
        del packet['source_sha256'][next(iter(packet['source_sha256']))]
    else:
        retained = pinned_receipts / NATIVE_ARCHIVE
        retained.write_text(retained.read_text() + '\n')
    path.write_text(json.dumps(packet))
    with pytest.raises(ValueError):
        ledger._particle_reconstruction_decisions(repo=pinned_receipts)
