#!/usr/bin/env python3
"""Tests for the compact-transient receipt scaffold and evidence boundary."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import pathlib
import subprocess
import sys
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "compact_transients" / "build_compact_transient_receipts.py"
SPEC = importlib.util.spec_from_file_location("compact_transient_receipts", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
receipts = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(receipts)


def sha256(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


COMMON_CONTENT: dict[str, Any] = {
    "assumptions": ["declared model domain"],
    "units": {"time": "s"},
    "parameters": {"model": "frozen-v1"},
    "outputs": {"result": sha256(b"checked result")},
}

RECEIPT_CONTENT: dict[str, dict[str, Any]] = {
    "COMPACT_QUOTIENT_RECEIPT": {
        "quotient_space": "Q_r^CR",
        "equivalence_relation": "Gamma_r^CR",
        "canonicalizer_deterministic": True,
        "readouts_factor_through_quotient": True,
    },
    "COMPACT_SOURCE_LAW_RECEIPT": {
        "channel": "test-channel",
        "source_law": "normalized mu_test",
        "normalization_checked": True,
        "normalization_error": 1e-12,
        "normalization_tolerance": 1e-9,
    },
    "PACKETIZED_KERNEL_RECEIPT": {
        "kernel": "K_test(dq,dPi,dell,dtau|q)",
        "normalization_checked": True,
        "coupled_sampler_checked": True,
        "joint_variables": ["state", "packet", "receipt", "clock"],
    },
    "PHYSICAL_CLOCK_RECEIPT": {
        "clock_map": "T_test",
        "calibration_checked": True,
        "input_unit": "step",
        "output_unit": "s",
    },
    "FINITE_PACKET_PARENT_RECEIPT": {
        "finite_schema_checked": True,
        "sectors": ["radio", "gamma", "GW", "optical", "neutrino", "environmental", "recipient"],
    },
    "PROPAGATION_RECEIPT": {
        "channel": "test-channel",
        "propagation_kernel": "P_test",
        "normalization_checked": True,
        "domain_codomain_checked": True,
    },
    "DETECTION_THINNING_RECEIPT": {
        "detector_model": "R_det_test",
        "exposure": "ObsWin_test",
        "fit_checked": True,
        "probability_bounds_checked": True,
    },
    "CENSORING_AND_UPPER_LIMIT_RECEIPT": {
        "censoring_model": "Cens_test",
        "zero_count_exposure_checked": True,
        "known_upper_limits_checked": True,
        "no_double_counting_checked": True,
    },
    "POINT_PROCESS_LIKELIHOOD_RECEIPT": {
        "likelihood": "marked Poisson test likelihood",
        "compensator_included": True,
        "event_count": 3,
        "log_likelihood": -4.5,
    },
    "HELDOUT_LIKELIHOOD_RECEIPT": {
        "proper_score": "held-out log score",
        "preregistration_id": "test-prereg-v1",
        "split_checked": True,
        "heldout_event_count": 2,
        "heldout_score": -2.5,
    },
    "CONTROL_MODEL_RECEIPT": {
        "controls": ["M0", "M1"],
        "acceptance_threshold": 1.0,
        "preregistered": True,
        "absolute_calibration_checked": True,
        "multiple_testing_checked": True,
    },
    "REFINEMENT_STABILITY_RECEIPT": {
        "regulators": ["r1", "r2"],
        "maximum_defect": 1e-5,
        "defect_tolerance": 1e-4,
        "stability_checked": True,
    },
    "FROZEN_HASHES_RECEIPT": {
        "freeze_id": "test-freeze-v1",
        "immutable_checked": True,
    },
}


def write_json(path: pathlib.Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_evidence_bundle(
    root: pathlib.Path, receipt_names: tuple[str, ...]
) -> tuple[pathlib.Path, dict[str, Any], dict[str, pathlib.Path]]:
    evidence_dir = root / "external_evidence"
    input_dir = evidence_dir / "inputs"
    input_dir.mkdir(parents=True)
    offers: dict[str, Any] = {}
    artifact_paths: dict[str, pathlib.Path] = {}
    for receipt in receipt_names:
        input_count = 2 if receipt == "HELDOUT_LIKELIHOOD_RECEIPT" else 1
        input_records = []
        for index in range(input_count):
            input_path = input_dir / f"{receipt.lower()}-{index}.json"
            input_raw = json.dumps({"receipt": receipt, "input": index}).encode()
            input_path.write_bytes(input_raw)
            input_records.append(
                {"path": str(input_path.relative_to(evidence_dir)), "sha256": sha256(input_raw)}
            )
        content = {**copy.deepcopy(COMMON_CONTENT), **copy.deepcopy(RECEIPT_CONTENT[receipt])}
        if receipt == "FROZEN_HASHES_RECEIPT":
            content["frozen_artifacts"] = {
                record["path"]: record["sha256"] for record in input_records
            }
        artifact = {
            "schema_version": receipts.EVIDENCE_SCHEMA,
            "receipt": receipt,
            "kind": receipts.EVIDENCE_KINDS[receipt],
            "producer": "independent-test-producer/v1",
            "inputs": input_records,
            "content": content,
        }
        artifact_path = evidence_dir / f"artifact-{len(artifact_paths)}.json"
        write_json(artifact_path, artifact)
        offers[receipt] = {
            "path": str(artifact_path.relative_to(evidence_dir)),
            "sha256": sha256(artifact_path.read_bytes()),
        }
        artifact_paths[receipt] = artifact_path
    manifest = {"schema_version": receipts.EVIDENCE_MANIFEST_SCHEMA, "receipts": offers}
    manifest_path = evidence_dir / "manifest.json"
    write_json(manifest_path, manifest)
    return manifest_path, manifest, artifact_paths


def run_builder(
    tmp_path: pathlib.Path,
    *,
    evidence_manifest: pathlib.Path | None = None,
    config: pathlib.Path | None = None,
    env: dict[str, str] | None = None,
) -> tuple[dict[str, Any], dict[str, Any], pathlib.Path]:
    out_dir = tmp_path / "output"
    command = [sys.executable, str(SCRIPT), "--output", str(out_dir)]
    if evidence_manifest is not None:
        command.extend(["--evidence-manifest", str(evidence_manifest)])
    if config is not None:
        command.extend(["--config", str(config)])
    subprocess.run(command, check=True, cwd=ROOT, env=env, capture_output=True, text=True)
    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    audit = json.loads((out_dir / "promotion_audit.json").read_text(encoding="utf-8"))
    return manifest, audit, out_dir


def rewrite_artifact(
    receipt: str,
    manifest_path: pathlib.Path,
    manifest: dict[str, Any],
    artifact_paths: dict[str, pathlib.Path],
    mutate: Any,
) -> None:
    path = artifact_paths[receipt]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    mutate(artifact)
    write_json(path, artifact)
    manifest["receipts"][receipt]["sha256"] = sha256(path.read_bytes())
    write_json(manifest_path, manifest)


def test_build_compact_transient_receipts_defaults_to_stated_cr1_boundary(
    tmp_path: pathlib.Path,
) -> None:
    manifest, audit, out_dir = run_builder(tmp_path)

    assert manifest["milestone"] == "COMPACT_TRANSIENT_RECEIPT_SCAFFOLD"
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert manifest["first_blocked_gate"] == "COMPACT_QUOTIENT_RECEIPT"
    assert manifest["promotion_allowed"] is False
    assert manifest["physical_claim"] is False
    assert manifest["missing_files"] == []
    assert [
        name for name, passed in audit["readiness_gates"].items() if passed
    ] == ["NO_GENERATION_LEAKAGE_RECEIPT"]
    assert len(audit["readiness_gates"]) == 27
    assert audit["verdict_checking_tier"] == "SCHEMA_AND_PROVENANCE_ONLY"
    assert audit["scientific_validation_performed"] is False

    for rel_path in manifest["required_files"]:
        assert (out_dir / rel_path).is_file(), rel_path
    assert "manifest.json" not in manifest["file_hashes"]
    for rel_path, expected in manifest["file_hashes"].items():
        assert sha256((out_dir / rel_path).read_bytes()) == expected, rel_path

    frb = json.loads((out_dir / "frb_controls.json").read_text(encoding="utf-8"))
    assert frb["controls"]["M2"] == "young_plus_old_gc_repair_reload_timing"

    bh = json.loads((out_dir / "bh_recycling.json").read_text(encoding="utf-8"))
    assert bh["genealogy_dag_required"] is True
    assert "ringdown_residual" in bh["forbidden_path"]


def test_valid_external_artifacts_close_only_the_cr2_receipts(tmp_path: pathlib.Path) -> None:
    manifest_path, _, _ = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)

    assert manifest["strongest_allowed_claim"] == "CR2_CONDITIONAL_PHENOMENOLOGY"
    assert manifest["first_blocked_gate"] == "CONTROLS"
    assert all(audit["readiness_gates"][name] for name in receipts.CR2_RECEIPTS)


def test_cr3_receipts_use_the_same_external_evidence_mechanism(tmp_path: pathlib.Path) -> None:
    names = receipts.CR2_RECEIPTS + (
        "CONTROL_MODEL_RECEIPT",
        "REFINEMENT_STABILITY_RECEIPT",
        "FROZEN_HASHES_RECEIPT",
    )
    manifest_path, _, _ = write_evidence_bundle(tmp_path, names)
    manifest, _, _ = run_builder(tmp_path, evidence_manifest=manifest_path)

    assert manifest["strongest_allowed_claim"] == "CR3_FROZEN_PHYSICAL_PREDICTION"
    assert manifest["first_blocked_gate"] == "COMPACT_SOURCE_ACTION_DERIVED_RECEIPT"


def test_non_object_offer_is_reported_invalid_without_crashing(tmp_path: pathlib.Path) -> None:
    """The offer-shape guard provides a stable diagnostic and prevents a crash."""
    manifest_path, evidence, _ = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    evidence["receipts"]["COMPACT_SOURCE_LAW_RECEIPT"] = "present"
    write_json(manifest_path, evidence)

    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    check = audit["evidence_checks"]["COMPACT_SOURCE_LAW_RECEIPT"]
    assert check == {
        "status": "INVALID",
        "reason": "offer must contain only path and sha256",
    }


def test_empty_artifact_reports_a_precise_reason(tmp_path: pathlib.Path) -> None:
    """The empty-file guard distinguishes emptiness from malformed JSON."""
    manifest_path, evidence, paths = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    target = "COMPACT_SOURCE_LAW_RECEIPT"
    paths[target].write_bytes(b"")
    evidence["receipts"][target]["sha256"] = sha256(b"")
    write_json(manifest_path, evidence)

    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert audit["evidence_checks"][target]["reason"] == "evidence artifact is empty"


def test_wrong_kind_artifact_cannot_close_a_receipt(tmp_path: pathlib.Path) -> None:
    manifest_path, evidence, paths = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    target = "COMPACT_SOURCE_LAW_RECEIPT"
    rewrite_artifact(target, manifest_path, evidence, paths, lambda artifact: artifact.update(kind="physical_clock_calibration"))

    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert audit["evidence_checks"][target]["reason"] == "evidence receipt or kind mismatch"


def test_stale_input_hash_cannot_close_a_receipt(tmp_path: pathlib.Path) -> None:
    manifest_path, _, paths = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    target = "COMPACT_SOURCE_LAW_RECEIPT"
    artifact = json.loads(paths[target].read_text(encoding="utf-8"))
    stale_input = paths[target].parent / artifact["inputs"][0]["path"]
    stale_input.write_text("changed after evidence generation\n", encoding="utf-8")

    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert audit["evidence_checks"][target]["reason"] == "recorded input hash is stale"


def test_circular_input_cannot_close_a_receipt(tmp_path: pathlib.Path) -> None:
    manifest_path, evidence, paths = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    target = "COMPACT_SOURCE_LAW_RECEIPT"

    def point_input_at_artifact(artifact: dict[str, Any]) -> None:
        artifact["inputs"] = [
            {
                "path": paths[target].name,
                "sha256": sha256(paths[target].read_bytes()),
            }
        ]

    rewrite_artifact(target, manifest_path, evidence, paths, point_input_at_artifact)
    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)

    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert audit["evidence_checks"][target]["reason"] == "input is missing, duplicate, or circular"


def test_builder_authored_artifact_cannot_close_a_receipt(tmp_path: pathlib.Path) -> None:
    manifest_path, evidence, paths = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    target = "COMPACT_SOURCE_LAW_RECEIPT"

    builder_output = tmp_path / "builder_output"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(builder_output)],
        check=True,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    generated_source_law = builder_output / "source_law.json"
    evidence["receipts"][target] = {
        "path": os.path.relpath(generated_source_law, manifest_path.parent),
        "sha256": sha256(generated_source_law.read_bytes()),
    }
    write_json(manifest_path, evidence)
    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert audit["evidence_checks"][target]["reason"] == "evidence schema version mismatch"

    evidence["receipts"][target] = {
        "path": str(paths[target].relative_to(manifest_path.parent)),
        "sha256": sha256(paths[target].read_bytes()),
    }
    rewrite_artifact(target, manifest_path, evidence, paths, lambda artifact: artifact.update(producer=receipts.BUILDER_ID))

    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert audit["evidence_checks"][target]["reason"] == "evidence producer is missing or self-authored"


def test_unsupported_legacy_enable_flags_are_reported_invalid(tmp_path: pathlib.Path) -> None:
    """Strict manifest shape rejects rather than silently ignores legacy flags."""
    evidence_dir = tmp_path / "external_evidence"
    manifest_path = evidence_dir / "manifest.json"
    write_json(
        manifest_path,
        {
            "schema_version": receipts.EVIDENCE_MANIFEST_SCHEMA,
            "receipts": {},
            "enable_cr2": True,
        },
    )
    config = tmp_path / "config.json"
    write_json(config, {"enable_cr2": True})
    env = {**os.environ, "COMPACT_TRANSIENT_ENABLE_CR2": "1"}

    manifest, audit, _ = run_builder(
        tmp_path, evidence_manifest=manifest_path, config=config, env=env
    )
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert set(
        (check["status"], check["reason"])
        for check in audit["evidence_checks"].values()
    ) == {("INVALID", "manifest must contain only schema_version and receipts")}


def test_schema_conforming_but_meaningless_evidence_currently_reaches_cr2(
    tmp_path: pathlib.Path,
) -> None:
    manifest_path, evidence, paths = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    semantic_string_fields = {
        "COMPACT_QUOTIENT_RECEIPT": ("quotient_space", "equivalence_relation"),
        "COMPACT_SOURCE_LAW_RECEIPT": ("channel", "source_law"),
        "PACKETIZED_KERNEL_RECEIPT": ("kernel",),
        "PHYSICAL_CLOCK_RECEIPT": ("clock_map", "input_unit", "output_unit"),
        "PROPAGATION_RECEIPT": ("channel", "propagation_kernel"),
        "DETECTION_THINNING_RECEIPT": ("detector_model", "exposure"),
        "CENSORING_AND_UPPER_LIMIT_RECEIPT": ("censoring_model",),
        "POINT_PROCESS_LIKELIHOOD_RECEIPT": ("likelihood",),
        "HELDOUT_LIKELIHOOD_RECEIPT": ("proper_score", "preregistration_id"),
    }

    for receipt in receipts.CR2_RECEIPTS:
        def replace_meaning(artifact: dict[str, Any], receipt: str = receipt) -> None:
            artifact["producer"] = "arbitrary external string"
            content = artifact["content"]
            content["assumptions"] = ["plausible but meaningless"]
            content["units"] = {"plausible but meaningless": "plausible but meaningless"}
            content["parameters"] = {"plausible but meaningless": "plausible but meaningless"}
            for field in semantic_string_fields.get(receipt, ()):
                content[field] = "plausible but meaningless"

        rewrite_artifact(receipt, manifest_path, evidence, paths, replace_meaning)

    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)
    assert manifest["strongest_allowed_claim"] == "CR2_CONDITIONAL_PHENOMENOLOGY"
    assert audit["verdict_checking_tier"] == "SCHEMA_AND_PROVENANCE_ONLY"
    assert audit["scientific_validation_performed"] is False


def test_schema_conforming_but_physically_meaningless_numerics_currently_reach_cr2(
    tmp_path: pathlib.Path,
) -> None:
    manifest_path, evidence, paths = write_evidence_bundle(tmp_path, receipts.CR2_RECEIPTS)
    numeric_changes = {
        "COMPACT_SOURCE_LAW_RECEIPT": {
            "normalization_error": 1e9,
            "normalization_tolerance": 1e9,
        },
        "POINT_PROCESS_LIKELIHOOD_RECEIPT": {
            "event_count": 0,
            "log_likelihood": 1e9,
        },
        "HELDOUT_LIKELIHOOD_RECEIPT": {
            "heldout_event_count": 1,
            "heldout_score": 1e9,
        },
    }
    for receipt, changes in numeric_changes.items():
        rewrite_artifact(
            receipt,
            manifest_path,
            evidence,
            paths,
            lambda artifact, changes=changes: artifact["content"].update(changes),
        )

    manifest, audit, _ = run_builder(tmp_path, evidence_manifest=manifest_path)
    assert manifest["strongest_allowed_claim"] == "CR2_CONDITIONAL_PHENOMENOLOGY"
    assert audit["verdict_checking_tier"] == "SCHEMA_AND_PROVENANCE_ONLY"
    assert audit["scientific_validation_performed"] is False


def test_rejects_ringdown_residual_as_generation_prior_input(tmp_path: pathlib.Path) -> None:
    config = tmp_path / "source_config.json"
    write_json(config, {"inputs": ["ringdown_residual"]})

    manifest, audit, _ = run_builder(tmp_path, config=config)
    assert manifest["strongest_allowed_claim"] == "CR1_QUOTIENT_DIAGNOSTIC"
    assert manifest["first_blocked_gate"] == "NO_GENERATION_LEAKAGE_RECEIPT"
    assert "ringdown_residual" in manifest["target_leak_hits"]
    assert audit["readiness_gates"]["NO_GENERATION_LEAKAGE_RECEIPT"] is False
