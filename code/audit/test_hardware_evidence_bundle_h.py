"""Adversarial tests for the non-promoting class-H evidence scaffold (#325)."""

from __future__ import annotations

import copy
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
VERIFIER_PATH = REPO_ROOT / "tools" / "verify_hardware_evidence_bundle_h.py"
FIXTURE = (
    Path(__file__).resolve().parent
    / "fixtures"
    / "hardware_evidence_bundle_h"
    / "reference_nonphysical"
)

spec = importlib.util.spec_from_file_location("verify_hardware_evidence_bundle_h", VERIFIER_PATH)
verifier = importlib.util.module_from_spec(spec)
sys.modules["verify_hardware_evidence_bundle_h"] = verifier
spec.loader.exec_module(verifier)


def fixture_copy(tmp_path: Path) -> tuple[Path, dict]:
    root = tmp_path / "bundle"
    shutil.copytree(FIXTURE, root)
    bundle_path = root / "bundle.json"
    return bundle_path, json.loads(bundle_path.read_text(encoding="utf-8"))


def artifact_row(bundle: dict, artifact_id: str) -> dict:
    return next(row for row in bundle["artifacts"] if row["artifact_id"] == artifact_id)


def write_bundle(bundle_path: Path, bundle: dict) -> None:
    bundle["bundle_binding"]["binding_sha256"] = verifier.recompute_bundle_binding(bundle)
    bundle_path.write_text(json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_artifact(
    bundle_path: Path,
    bundle: dict,
    artifact_id: str,
    value,
) -> None:
    row = artifact_row(bundle, artifact_id)
    path = bundle_path.parent / row["path"]
    if isinstance(value, str):
        path.write_text(value, encoding="utf-8")
    else:
        path.write_text(
            json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )
    row["sha256"] = verifier.sha256_path(path)


def verify(
    bundle_path: Path,
    *,
    replay_registry: Path | None | bool = True,
) -> dict:
    if replay_registry is True:
        replay_registry = bundle_path.parent / "replay_registry.json"
    return verifier.verify_bundle(
        bundle_path,
        replay_registry_path=replay_registry if isinstance(replay_registry, Path) else None,
    )


def assert_rejected(report: dict, code: str, *, verdict: str = "INSUFFICIENT") -> None:
    assert report["schema_valid"] is True
    assert report["verdict"] == verdict
    assert code in report["rejection_codes"]


def make_declared_physical_attestation(
    bundle_path: Path,
    bundle: dict,
    *,
    mode: str = "independent_end_to_end_witness",
    compromised: bool = False,
    signer_count: int = 1,
) -> None:
    """Build a self-authored *declaration*, never real evidence, for attack tests."""
    bundle["claim_boundary"] = {
        "fixture_kind": "physical_bundle",
        "physical_claim": True,
        "theory_discriminating": True,
        "allowed_conclusion": "Adversarial fabricated physical claim; must not promote.",
    }
    bundle["claim"]["effect_statement"] = "The synthetic device levitated the Moon."
    bundle["claim"]["magnitude"] = {"value": 1.0, "unit": "moon"}
    bundle["claim"]["uncertainty"] = {
        "kind": "fabricated",
        "value": 0.0,
        "unit": "moon",
    }
    bundle["claim"]["extraordinary_effect"] = True
    write_artifact(
        bundle_path,
        bundle,
        "claim-text-001",
        "The synthetic device levitated the Moon.\n",
    )
    parties = [f"fake-independent-party-{index}" for index in range(signer_count)]
    signer_ids = [f"fake-signer-{index}" for index in range(signer_count)]
    bundle["attestation"] = {
        "mode": mode,
        "artifact_ids": ["attestation-001"],
        "party_ids": parties,
        "independent_of_claimant": True,
        "fresh_device_or_run_series": True,
        "witnessed_end_to_end": True,
        "signer_ids": signer_ids,
    }
    bundle["signers"] = [
        {
            "signer_id": signer_id,
            "party_id": party_id,
            "key_id": f"fake-key-{index}",
            "compromised": compromised,
            "signed_artifact_ids": ["attestation-001"],
        }
        for index, (signer_id, party_id) in enumerate(zip(signer_ids, parties))
    ]
    write_artifact(
        bundle_path,
        bundle,
        "attestation-001",
        {
            "bundle_id": bundle["bundle_id"],
            "mode": mode,
            "note": "Self-authored fixture bytes; no signature or trust root.",
        },
    )


def test_reference_fixture_is_integrity_valid_but_explicitly_nonphysical(tmp_path: Path) -> None:
    bundle_path, _ = fixture_copy(tmp_path)
    report = verify(bundle_path)
    assert report["integrity_valid"] is True
    assert_rejected(report, "NON_PHYSICAL_FIXTURE")
    assert "TRUST_ROOT_SIGNATURE_VERIFICATION_OPEN" in report["open_gates"]
    assert "FRESH_REPRODUCTION_BUNDLE_VERIFICATION_OPEN" in report["open_gates"]
    assert "ANALYSIS_TO_CLAIM_EXECUTION_OPEN" in report["open_gates"]
    assert report["ignored_producer_assertions"] == [
        "analysis_passed",
        "calibration_passed",
        "physical_effect_established",
    ]


def test_replay_registry_is_required_when_replay_is_in_scope(tmp_path: Path) -> None:
    bundle_path, _ = fixture_copy(tmp_path)
    report = verify(bundle_path, replay_registry=None)
    assert_rejected(report, "REPLAY_REGISTRY_REQUIRED")
    assert report["predicates"]["replay_protection"] is False


def test_seen_capture_nonce_is_rejected(tmp_path: Path) -> None:
    bundle_path, _ = fixture_copy(tmp_path)
    registry = tmp_path / "seen.json"
    registry.write_text(
        json.dumps(
            {
                "registry_id": "class-h-fixture-registry",
                "seen_nonces": ["nonce-synthetic-001"],
            }
        ),
        encoding="utf-8",
    )
    report = verify(bundle_path, replay_registry=registry)
    assert_rejected(report, "REPLAY_NONCE_SEEN")


def test_anonymous_replay_nonce_list_is_not_an_authoritative_registry(
    tmp_path: Path,
) -> None:
    bundle_path, _ = fixture_copy(tmp_path)
    registry = tmp_path / "anonymous-list.json"
    registry.write_text("[]\n", encoding="utf-8")
    assert_rejected(verify(bundle_path, replay_registry=registry), "REPLAY_REGISTRY_INVALID")


def test_selective_reporting_of_a_scheduled_run_is_rejected(tmp_path: Path) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    bundle["run_schedule"]["run_ids"].append("run-omitted-002")
    write_artifact(
        bundle_path,
        bundle,
        "schedule-001",
        {"run_ids": bundle["run_schedule"]["run_ids"]},
    )
    write_bundle(bundle_path, bundle)
    assert_rejected(verify(bundle_path), "SELECTIVE_REPORTING")


def test_stale_calibration_is_rejected(tmp_path: Path) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    calibration = bundle["calibrations"][0]
    calibration["valid_until_utc"] = "2026-07-20T11:59:59Z"
    write_artifact(
        bundle_path,
        bundle,
        "calibration-001",
        {
            "calibration_id": calibration["calibration_id"],
            "device_id": calibration["device_id"],
            "reference_id": calibration["reference_id"],
            "valid_from_utc": calibration["valid_from_utc"],
            "valid_until_utc": calibration["valid_until_utc"],
        },
    )
    write_bundle(bundle_path, bundle)
    assert_rejected(verify(bundle_path), "CALIBRATION_OUT_OF_WINDOW")


def test_calibration_metadata_cannot_drift_from_bound_artifact(tmp_path: Path) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    bundle["calibrations"][0]["reference_id"] = "self-authored-different-reference"
    write_bundle(bundle_path, bundle)
    assert_rejected(verify(bundle_path), "CALIBRATION_ARTIFACT_MISMATCH")


def test_device_substitution_is_rejected_even_when_raw_bytes_are_rebound(tmp_path: Path) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    run = bundle["runs"][0]
    run["device_id"] = "substituted-device"
    write_artifact(
        bundle_path,
        bundle,
        "raw-001",
        {
            "capture_nonce": run["capture_nonce"],
            "channel": "synthetic_voltage",
            "device_id": run["device_id"],
            "run_id": run["run_id"],
            "samples": [0.0, 0.0, 0.0],
        },
    )
    write_bundle(bundle_path, bundle)
    assert_rejected(verify(bundle_path), "DEVICE_SUBSTITUTION")


def test_custody_break_is_rejected(tmp_path: Path) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    bundle["custody"]["continuous_declared"] = False
    write_artifact(
        bundle_path,
        bundle,
        "custody-001",
        {
            "continuous_declared": False,
            "segments": bundle["custody"]["segments"],
        },
    )
    write_bundle(bundle_path, bundle)
    assert_rejected(verify(bundle_path), "CUSTODY_BREAK")


@pytest.mark.parametrize("artifact_id", ["analysis-001", "claim-text-001"])
def test_altered_analysis_or_claim_text_fails_hash_check(
    tmp_path: Path,
    artifact_id: str,
) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    row = artifact_row(bundle, artifact_id)
    path = bundle_path.parent / row["path"]
    path.write_bytes(path.read_bytes() + b"\nunauthorized mutation\n")
    report = verify(bundle_path)
    assert_rejected(report, "ARTIFACT_HASH_MISMATCH", verdict="INVALID")
    assert report["integrity_valid"] is False


def test_compromised_single_signer_cannot_attest(tmp_path: Path) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    make_declared_physical_attestation(bundle_path, bundle, compromised=True)
    write_bundle(bundle_path, bundle)
    report = verify(bundle_path)
    assert_rejected(report, "SIGNER_COMPROMISED")
    assert "TRUST_ROOT_SIGNATURE_VERIFICATION_OPEN" in report["rejection_codes"]


def test_multiple_signatures_over_same_capture_are_not_independent_attestation(
    tmp_path: Path,
) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    make_declared_physical_attestation(
        bundle_path,
        bundle,
        mode="same_capture_multi_signature",
        signer_count=2,
    )
    write_bundle(bundle_path, bundle)
    report = verify(bundle_path)
    assert_rejected(report, "SIGNATURES_NOT_INDEPENDENT_ATTESTATION")


def test_fabricated_physical_claim_with_recomputed_hashes_never_promotes(
    tmp_path: Path,
) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    make_declared_physical_attestation(bundle_path, bundle)
    # The attacker controls every declared metadata field, recomputes all file
    # hashes and the canonical root, and supplies an empty nonce registry.
    write_bundle(bundle_path, bundle)
    report = verify(bundle_path)
    assert report["integrity_valid"] is True
    assert report["verdict"] == "INSUFFICIENT"
    assert {
        "TRUST_ROOT_SIGNATURE_VERIFICATION_OPEN",
        "FRESH_REPRODUCTION_BUNDLE_VERIFICATION_OPEN",
        "ANALYSIS_TO_CLAIM_EXECUTION_OPEN",
    }.issubset(report["rejection_codes"])
    assert report["verdict"] != "SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL"


def test_producer_pass_booleans_are_ignored(tmp_path: Path) -> None:
    bundle_path, bundle = fixture_copy(tmp_path)
    before = verify(bundle_path)
    changed = copy.deepcopy(bundle)
    changed["producer_assertions"] = {
        key: not value for key, value in changed["producer_assertions"].items()
    }
    write_bundle(bundle_path, changed)
    after = verify(bundle_path)
    assert before["verdict"] == after["verdict"] == "INSUFFICIENT"
    assert before["predicates"] == after["predicates"]
    assert before["rejection_codes"] == after["rejection_codes"]
    assert before["ignored_producer_assertions"] == after["ignored_producer_assertions"]


def test_cli_exit_codes_fail_closed(tmp_path: Path) -> None:
    bundle_path, _ = fixture_copy(tmp_path)
    insufficient = subprocess.run(
        [
            sys.executable,
            str(VERIFIER_PATH),
            str(bundle_path),
            "--replay-registry",
            str(bundle_path.parent / "replay_registry.json"),
        ],
        capture_output=True,
        text=True,
    )
    assert insufficient.returncode == 1
    assert '"verdict": "INSUFFICIENT"' in insufficient.stdout

    malformed = tmp_path / "malformed.json"
    malformed.write_text("{}\n", encoding="utf-8")
    invalid = subprocess.run(
        [sys.executable, str(VERIFIER_PATH), str(malformed)],
        capture_output=True,
        text=True,
    )
    assert invalid.returncode == 2
    assert '"verdict": "INVALID"' in invalid.stdout
