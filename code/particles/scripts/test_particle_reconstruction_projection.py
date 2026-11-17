"""Cheap custody/scope tests; never run the structural or scientific builders."""
from copy import deepcopy
import hashlib
import json

import pytest

import build_postdiction_ledger as ledger


@pytest.fixture
def pinned_receipts(tmp_path):
    source = tmp_path/"source.py"
    source.write_text("# fixture source\n")
    pins = {"source.py": hashlib.sha256(source.read_bytes()).hexdigest()}
    inputs = {
        "native_defect_decision_results.json": ("scope", "decision", "sectors", "collision_candidate",
                                                "quantum_collision_survival", "transport_cost",
                                                "corrected_memory_control"),
        "particle_identification/results.json": ("verdict", "finite_inverse_response_kinetics"),
        "chiral_candidate/record_two_point_results.json": ("scope", "joint_response", "held_out_validation"),
        "chiral_candidate/results.json": ("verdict", "kinematics"),
    }
    for name, fields in inputs.items():
        original = json.loads((ledger.REPO/"model/common_source"/name).read_text())
        packet = {key: deepcopy(original[key]) for key in ("schema", *fields)}
        packet["source_bindings"] = pins
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
