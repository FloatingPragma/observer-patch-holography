from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


pytest.importorskip("numpy")

PROGRAMS = Path(__file__).resolve().parents[1] / "programs"
sys.path.insert(0, str(PROGRAMS))

import cayley_blind_likelihood_analysis as analysis  # noqa: E402
import record_gated_cayley_runtime as runtime  # noqa: E402
import runtime_analysis_packet as adapter  # noqa: E402


ANALYSIS_TEST_PATH = Path(__file__).with_name("test_cayley_blind_likelihood_analysis.py")
ANALYSIS_TEST_SPEC = importlib.util.spec_from_file_location(
    "runtime_packet_analysis_fixture", ANALYSIS_TEST_PATH
)
assert ANALYSIS_TEST_SPEC is not None and ANALYSIS_TEST_SPEC.loader is not None
analysis_fixture = importlib.util.module_from_spec(ANALYSIS_TEST_SPEC)
sys.modules[ANALYSIS_TEST_SPEC.name] = analysis_fixture
ANALYSIS_TEST_SPEC.loader.exec_module(analysis_fixture)


def _runtime_fixture(tmp_path: Path):
    lock = analysis_fixture._build_lock(shots=192)
    lock = json.loads(json.dumps(lock))
    for row in lock["expected_rows"]:
        row["backend_role"] = (
            "development" if row["backend_name"] == "heldout_backend_a" else "heldout"
        )
    unhashed = dict(lock)
    unhashed.pop("analysis_lock_sha256")
    lock["analysis_lock_sha256"] = analysis.sha256_json(unhashed)
    analysis.validate_analysis_lock(lock, verify_code_hash=True)

    manifest_sha = "e" * 64
    catalog_sha = lock["catalog_precommitment_sha256"]
    expected_by_id = {row["opaque_id"]: row for row in lock["expected_rows"]}
    diagnostic_to_calibration = {
        opaque_id: (calibration_id, calibration)
        for calibration_id, calibration in lock["calibrations"].items()
        for opaque_id in calibration["control_rule"]["diagnostic_opaque_ids"]
    }
    public_rows = []
    descriptors = {}
    families = {}
    for opaque_id, row in expected_by_id.items():
        public_rows.append(
            {
                "opaque_id": opaque_id,
                "logical_circuit_sha256": row["logical_circuit_sha256"],
                "shots": row["shots"],
                "backend_role": f"opaque-{row['backend_role']}",
            }
        )
        families[opaque_id] = row["family"]
        descriptors[opaque_id] = {"family": row["family"], "shots": row["shots"]}
    for opaque_id, (_, calibration) in diagnostic_to_calibration.items():
        logical_sha = runtime.sha256_json({"diagnostic": opaque_id})
        code = calibration["control_rule"]["expected_basis_code_by_opaque_id"][opaque_id]
        public_rows.append(
            {
                "opaque_id": opaque_id,
                "logical_circuit_sha256": logical_sha,
                "shots": 192,
                "backend_role": (
                    "opaque-development" if opaque_id.endswith("-a") else "opaque-heldout"
                ),
            }
        )
        families[opaque_id] = "readout_calibration"
        descriptors[opaque_id] = {
            "family": "readout_calibration",
            "shots": 192,
            "parameters": {
                "basis_code": code,
                "num_qubits": 4,
                "evidentiary_role": "diagnostic_only",
            },
        }
    manifest = {
        "catalog_precommitment_sha256": catalog_sha,
        "circuits": public_rows,
    }
    bundle = SimpleNamespace(
        analysis_lock_document=lock,
        analysis_lock_sha256=lock["analysis_lock_sha256"],
        manifest=manifest,
        manifest_sha256=manifest_sha,
        circuit_families=families,
        circuit_descriptors=descriptors,
    )

    public_by_id = {row["opaque_id"]: row for row in public_rows}
    submission_paths = {
        role: tmp_path / f"runtime_{role}_submission_events.ndjson"
        for role in adapter.SUBMISSION_ROLES
    }
    registrations = {}
    for opaque_id, public in public_by_id.items():
        expected = expected_by_id.get(opaque_id)
        if expected is not None:
            role = expected["backend_role"]
            backend_name = expected["backend_name"]
            layout_id = expected["layout_id"]
            layout = expected["physical_layout"]
            family = expected["family"]
            shots = expected["shots"]
        else:
            role = "development" if opaque_id.endswith("-a") else "heldout"
            reference = next(
                row for row in lock["expected_rows"] if row["backend_role"] == role
            )
            backend_name = reference["backend_name"]
            layout_id = reference["layout_id"]
            layout = reference["physical_layout"]
            family = "readout_calibration"
            shots = 192
        group_id = runtime.sha256_json({"group": opaque_id})
        job_id = f"job-{opaque_id}"
        binding = {
            "opaque_id": opaque_id,
            "logical_circuit_sha256": public["logical_circuit_sha256"],
            "compiled_qpy_sha256": runtime.sha256_json({"compiled": opaque_id}),
        }
        common = {
            "schema_version": runtime.RUNTIME_SCHEMA,
            "timestamp_utc": "2026-07-11T00:00:00+00:00",
            "manifest_sha256": manifest_sha,
            "analysis_lock_sha256": lock["analysis_lock_sha256"],
            "operator_source_sha256": "b" * 64,
            "plan_sha256": "a" * 64,
            "group_id": group_id,
            "backend_name": backend_name,
            "backend_role": role,
            "backend_role_opaque_id": f"opaque-{role}",
            "family": family,
            "layout_id": layout_id,
            "physical_layout": layout,
            "properties_last_update": "2026-07-11T00:00:00+00:00",
            "shots": shots,
            "opaque_ids": [opaque_id],
            "circuit_bindings": [binding],
        }
        runtime.append_event(
            submission_paths[role],
            {**common, "event_type": "submission_started"},
        )
        registration = runtime.append_event(
            submission_paths[role],
            {
                **common,
                "event_type": "job_registered",
                "job_id": job_id,
                "status": "DONE",
            },
        )
        registrations[opaque_id] = registration

    harvest_path = tmp_path / "runtime_harvest_events.ndjson"
    for opaque_id, registration in registrations.items():
        expected = expected_by_id.get(opaque_id)
        if expected is not None:
            heated = expected["valid_codes"][0]
            final = expected["valid_codes"][0]
            joined = {f"{final:03b}0{heated:03b}": expected["shots"]}
            width = 7
            layout_text = "final[3]|decision[1]|heated[3]"
        else:
            code = descriptors[opaque_id]["parameters"]["basis_code"]
            joined = {f"{code:04b}": 192}
            width = 4
            layout_text = "calibration_c[4]"
        binding = registration["circuit_bindings"][0]
        result = {
            "opaque_id": opaque_id,
            "joint_counts": joined,
            "raw_joined_counts_sha256": runtime.sha256_json(joined),
            "logical_circuit_sha256": binding["logical_circuit_sha256"],
            "compiled_qpy_sha256": binding["compiled_qpy_sha256"],
            "joined_count_width": width,
            "joined_count_layout": layout_text,
        }
        event = {
            "schema_version": runtime.RUNTIME_SCHEMA,
            "event_type": "job_harvest",
            "timestamp_utc": "2026-07-11T00:11:00+00:00",
            "manifest_sha256": manifest_sha,
            "analysis_lock_sha256": lock["analysis_lock_sha256"],
            "operator_source_sha256": registration["operator_source_sha256"],
            "harvester_source_sha256": "c" * 64,
            "plan_sha256": registration["plan_sha256"],
            "submission_event_sha256": registration["event_sha256"],
            "job_id": registration["job_id"],
            "group_id": registration["group_id"],
            "backend_name": registration["backend_name"],
            "backend_role": registration["backend_role"],
            "backend_role_opaque_id": registration["backend_role_opaque_id"],
            "family": registration["family"],
            "layout_id": registration["layout_id"],
            "physical_layout": registration["physical_layout"],
            "properties_last_update": registration["properties_last_update"],
            "shots": registration["shots"],
            "opaque_ids": registration["opaque_ids"],
            "circuit_bindings": registration["circuit_bindings"],
            "status": "DONE",
            "metrics": {"timestamps": {"running": "2026-07-11T00:10:00+00:00"}},
            "calibration_at_execution": {
                "properties": {"last_update_date": "2026-07-11T00:00:00+00:00"},
                "properties_sha256": "f" * 64,
            },
            "results": [result],
        }
        runtime.append_event(harvest_path, event)
    return bundle, submission_paths, harvest_path


def test_runtime_journals_seal_into_fully_validated_hardware_packet(tmp_path: Path) -> None:
    bundle, submissions, harvest = _runtime_fixture(tmp_path)
    packet = adapter.seal_runtime_analysis_packet(
        bundle=bundle,
        submission_journals=submissions,
        harvest_journal=harvest,
        created_utc="2026-07-11T01:00:00+00:00",
    )
    assert packet["source_kind"] == "ibm_qpu_hardware"
    assert len(packet["rows"]) == 3
    assert set(packet["calibration_results"]) == {"cal_backend_a", "cal_backend_b"}
    assert all(result["gof_p_value"] == 1.0 for result in packet["calibration_results"].values())
    analysis.validate_data_packet(bundle.analysis_lock_document, packet)


def test_packet_adapter_rejects_omitted_harvested_job(tmp_path: Path) -> None:
    bundle, submissions, harvest = _runtime_fixture(tmp_path)
    lines = harvest.read_text(encoding="utf-8").splitlines()
    harvest.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
    with pytest.raises(adapter.PacketAdapterError, match="omits"):
        adapter.seal_runtime_analysis_packet(
            bundle=bundle,
            submission_journals=submissions,
            harvest_journal=harvest,
        )


def test_packet_adapter_rejects_tampered_raw_counts(tmp_path: Path) -> None:
    bundle, submissions, harvest = _runtime_fixture(tmp_path)
    lines = harvest.read_text(encoding="utf-8").splitlines()
    row = json.loads(lines[0])
    row["results"][0]["joint_counts"] = {"0000000": 191}
    lines[0] = json.dumps(row, sort_keys=True, separators=(",", ":"))
    harvest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(runtime.RuntimeSafetyError, match="hash chain"):
        adapter.seal_runtime_analysis_packet(
            bundle=bundle,
            submission_journals=submissions,
            harvest_journal=harvest,
        )


def _dynamic_row_execution(*, family: str, layout: str):
    expected = {
        "row_id": "opaque-dynamic",
        "opaque_id": "opaque-dynamic",
        "logical_circuit_sha256": "d" * 64,
        "family": family,
        "backend_role": "development",
        "backend_name": "ibm-test",
        "layout_id": "layout-test",
        "physical_layout": [1, 2, 3, 4],
        "valid_codes": [0, 1],
        "shots": 3,
    }
    harvest = {
        "family": family,
        "backend_role": expected["backend_role"],
        "backend_name": expected["backend_name"],
        "layout_id": expected["layout_id"],
        "physical_layout": expected["physical_layout"],
        "shots": expected["shots"],
        "job_id": "job-dynamic",
        "group_id": "group-dynamic",
        "submission_event_sha256": "a" * 64,
        "event_sha256": "b" * 64,
        "metrics": {"timestamps": {"running": "2026-07-11T00:10:00+00:00"}},
        "calibration_at_execution": {
            "properties": {"last_update_date": "2026-07-11T00:00:00+00:00"}
        },
    }
    joint_counts = {"0011000": 3}
    result = {
        "opaque_id": expected["opaque_id"],
        "logical_circuit_sha256": expected["logical_circuit_sha256"],
        "compiled_qpy_sha256": "c" * 64,
        "joint_counts": joint_counts,
        "raw_joined_counts_sha256": runtime.sha256_json(joint_counts),
        "joined_count_width": 7,
        "joined_count_layout": layout,
    }
    return expected, harvest, result


def test_dynamic_packet_rows_accept_the_frozen_mh_joined_layout() -> None:
    expected, harvest, result = _dynamic_row_execution(
        family="mh", layout="after[3]|accepted[1]|before[3]"
    )
    rows = adapter._dynamic_rows(
        {expected["opaque_id"]: expected},
        {expected["opaque_id"]: (harvest, result)},
    )
    assert len(rows) == 1
    assert rows[0]["retrieved_shots"] == 3


def test_dynamic_packet_rows_reject_a_layout_from_the_wrong_family() -> None:
    expected, harvest, result = _dynamic_row_execution(
        family="mh", layout="final[3]|decision[1]|heated[3]"
    )
    with pytest.raises(adapter.PacketAdapterError, match="frozen seven-bit"):
        adapter._dynamic_rows(
            {expected["opaque_id"]: expected},
            {expected["opaque_id"]: (harvest, result)},
        )


def test_packet_cli_accepts_frozen_creation_time_for_exact_replay() -> None:
    args = adapter.parse_args(
        [
            "--manifest",
            "manifest.json",
            "--reveal",
            "reveal.json",
            "--analysis-lock",
            "analysis-lock.json",
            "--development-submission-journal",
            "development.ndjson",
            "--heldout-submission-journal",
            "heldout.ndjson",
            "--harvest-journal",
            "harvest.ndjson",
            "--created-utc",
            "2026-07-11T06:14:46.204442+00:00",
            "--json-out",
            "packet.json",
        ]
    )
    assert args.created_utc == "2026-07-11T06:14:46.204442+00:00"
