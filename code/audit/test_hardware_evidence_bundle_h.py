"""Adversarial tests for the non-promoting class-H evidence scaffold (#325)."""

from __future__ import annotations

import copy
import base64
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from tools import hardware_evidence_external as external

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
    trust_policy: Path | None = None,
    external_evidence: Path | None = None,
) -> dict:
    if replay_registry is True:
        replay_registry = bundle_path.parent / "replay_registry.json"
    return verifier.verify_bundle(
        bundle_path,
        replay_registry_path=replay_registry if isinstance(replay_registry, Path) else None,
        trust_policy_path=trust_policy,
        external_evidence_path=external_evidence,
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


def make_authenticated_contract_fixture(tmp_path: Path) -> dict:
    """Build a non-vacuous policy fixture without claiming a real experiment."""
    bundle_path, bundle = fixture_copy(tmp_path)
    bundle["bundle_id"] = "class-h-authenticated-contract-v1"
    effect = "Mean synthetic voltage is zero."
    bundle["claim"] = {
        "class": "H",
        "device_id": "synthetic-device-001",
        "protocol_id": "protocol-synthetic-001",
        "conditions": "Decision-procedure fixture; no real experiment is asserted.",
        "effect_statement": effect,
        "magnitude": {"value": 0.0, "unit": "synthetic_volt"},
        "uncertainty": {
            "kind": "max_absolute_deviation",
            "value": 0.0,
            "unit": "synthetic_volt",
        },
        "extraordinary_effect": False,
    }
    bundle["claim_boundary"] = {
        "fixture_kind": "physical_bundle",
        "physical_claim": True,
        "theory_discriminating": False,
        "allowed_conclusion": (
            "The decision procedure has a satisfiable evidence-policy path; "
            "this synthetic fixture is not evidence of a real device."
        ),
    }
    bundle["attestation"] = {
        "mode": "independent_end_to_end_witness",
        "artifact_ids": ["attestation-001"],
        "party_ids": ["independent-witness"],
        "independent_of_claimant": True,
        "fresh_device_or_run_series": False,
        "witnessed_end_to_end": True,
        "signer_ids": ["independent-witness-signer"],
    }
    bundle["signers"] = [
        {
            "signer_id": "independent-witness-signer",
            "party_id": "independent-witness",
            "key_id": "key-independent-attestor",
            "compromised": False,
            "signed_artifact_ids": ["attestation-001"],
        }
    ]
    write_artifact(
        bundle_path,
        bundle,
        "analysis-001",
        {
            "schema": "oph.hardware_evidence_bundle_h.analysis.mean_max_deviation.v1",
            "raw_artifact_ids": ["raw-001"],
            "sample_field": "samples",
            "effect_statement": effect,
            "unit": "synthetic_volt",
            "uncertainty_kind": "max_absolute_deviation",
        },
    )
    artifact_row(bundle, "analysis-001")["media_type"] = (
        "application/vnd.oph.hardware-analysis+json"
    )
    write_artifact(bundle_path, bundle, "claim-text-001", effect + "\n")
    write_artifact(
        bundle_path,
        bundle,
        "attestation-001",
        {
            "bundle_id": bundle["bundle_id"],
            "mode": "independent_end_to_end_witness",
            "witnessed_scope": [
                "device identity",
                "capture",
                "controls",
                "custody",
            ],
        },
    )

    write_bundle(bundle_path, bundle)
    registry_path = bundle_path.parent / "replay_registry.json"
    assignments = [
        {
            "nonce": bundle["runs"][0]["capture_nonce"],
            "bundle_id": bundle["bundle_id"],
            "binding_sha256": bundle["bundle_binding"]["binding_sha256"],
            "registered_utc": "2026-07-20T10:30:00Z",
            "state": "consumed",
            "consumed_utc": "2026-07-20T12:10:00Z",
        }
    ]
    registry_path.write_text(
        json.dumps(
            {
                "registry_id": bundle["replay_protection"]["registry_id"],
                "seen_nonces": [],
                "assignments": assignments,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    role_parties = {
        "claimant": ("synthetic-claimant", "claimant-organization"),
        "measurement_authority": (
            "measurement-authority",
            "measurement-organization",
        ),
        "calibration_authority": (
            "calibration-authority",
            "calibration-organization",
        ),
        "device_authority": ("device-authority", "device-organization"),
        "custody_authority": ("custody-authority", "custody-organization"),
        "preregistration_authority": (
            "preregistration-authority",
            "registry-organization",
        ),
        "replay_authority": ("replay-authority", "registry-organization"),
        "independent_attestor": (
            "independent-witness",
            "independent-witness-organization",
        ),
    }
    private_keys = {
        role: Ed25519PrivateKey.generate() for role in role_parties
    }
    anchors = []
    for role, (party_id, organization_id) in role_parties.items():
        key_id = f"key-{role.replace('_', '-')}"
        public_bytes = private_keys[role].public_key().public_bytes(
            encoding=serialization.Encoding.Raw,
            format=serialization.PublicFormat.Raw,
        )
        anchors.append(
            {
                "key_id": key_id,
                "party_id": party_id,
                "organization_id": organization_id,
                "role": role,
                "public_key_ed25519": base64.b64encode(public_bytes).decode("ascii"),
                "valid_from_utc": "2026-01-01T00:00:00Z",
                "valid_until_utc": "2027-01-01T00:00:00Z",
                "revoked": False,
            }
        )

    def sign(role: str, payload: dict) -> str:
        return base64.b64encode(
            private_keys[role].sign(external.canonical_bytes(payload))
        ).decode("ascii")

    signatures = []

    def add_signature(
        role: str,
        subject_kind: str,
        subject_id: str,
        digest: str,
    ) -> None:
        key_id = f"key-{role.replace('_', '-')}"
        payload = external.signature_payload(
            subject_kind=subject_kind,
            subject_id=subject_id,
            sha256=digest,
            key_id=key_id,
        )
        signatures.append(
            {
                **payload,
                "signature_base64": sign(role, payload),
            }
        )

    add_signature(
        "claimant",
        "bundle_binding",
        bundle["bundle_id"],
        bundle["bundle_binding"]["binding_sha256"],
    )
    add_signature(
        "independent_attestor",
        "bundle_binding",
        bundle["bundle_id"],
        bundle["bundle_binding"]["binding_sha256"],
    )
    for row in bundle["artifacts"]:
        role = external.ARTIFACT_AUTHORITY.get(row["role"])
        if role is not None:
            add_signature(role, "artifact", row["artifact_id"], row["sha256"])

    commitments = []
    for kind, artifact_id in (
        ("run_schedule", bundle["run_schedule"]["artifact_id"]),
        ("analysis", bundle["analysis_binding"]["analysis_artifact_id"]),
    ):
        key_id = "key-preregistration-authority"
        payload = external.commitment_payload(
            kind=kind,
            artifact_id=artifact_id,
            sha256=artifact_row(bundle, artifact_id)["sha256"],
            committed_utc="2026-07-20T10:00:00Z",
            key_id=key_id,
            policy_id="class-h-test-policy",
            bundle_id=bundle["bundle_id"],
            claimant_id=bundle["claimant_id"],
            device_id=bundle["claim"]["device_id"],
            protocol_id=bundle["claim"]["protocol_id"],
        )
        commitments.append(
            {
                **payload,
                "signature_base64": sign("preregistration_authority", payload),
            }
        )

    registry_payload = external.registry_payload(
        registry_id=bundle["replay_protection"]["registry_id"],
        snapshot_utc="2026-07-20T12:15:00Z",
        seen_nonces=[],
        assignments=assignments,
        key_id="key-replay-authority",
    )
    evidence = {
        "schema": external.EVIDENCE_SCHEMA,
        "policy_id": "class-h-test-policy",
        "signatures": signatures,
        "commitments": commitments,
        "nonce_registry": {
            **registry_payload,
            "signature_base64": sign("replay_authority", registry_payload),
        },
    }
    policy = {
        "schema": external.TRUST_SCHEMA,
        "policy_id": "class-h-test-policy",
        "anchors": anchors,
    }
    policy_path = bundle_path.parent / "operator_trust_policy.json"
    evidence_path = bundle_path.parent / "external_evidence.json"
    policy_path.write_text(
        json.dumps(policy, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    evidence_path.write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "bundle_path": bundle_path,
        "bundle": bundle,
        "registry_path": registry_path,
        "policy_path": policy_path,
        "evidence_path": evidence_path,
        "policy": policy,
        "evidence": evidence,
        "private_keys": private_keys,
    }


def write_external(context: dict) -> None:
    context["policy_path"].write_text(
        json.dumps(context["policy"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    context["evidence_path"].write_text(
        json.dumps(context["evidence"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def refresh_authenticated_signatures(context: dict) -> None:
    """Refresh trusted signatures after an intentional test mutation."""
    write_bundle(context["bundle_path"], context["bundle"])
    role_by_key = {
        f"key-{role.replace('_', '-')}": role
        for role in context["private_keys"]
    }
    for row in context["evidence"]["signatures"]:
        if row["subject_kind"] == "bundle_binding":
            row["sha256"] = context["bundle"]["bundle_binding"]["binding_sha256"]
        else:
            row["sha256"] = artifact_row(
                context["bundle"],
                row["subject_id"],
            )["sha256"]
        payload = {key: value for key, value in row.items() if key != "signature_base64"}
        role = role_by_key[row["key_id"]]
        row["signature_base64"] = base64.b64encode(
            context["private_keys"][role].sign(external.canonical_bytes(payload))
        ).decode("ascii")
    for row in context["evidence"]["commitments"]:
        row["sha256"] = artifact_row(
            context["bundle"],
            row["artifact_id"],
        )["sha256"]
        payload = {key: value for key, value in row.items() if key != "signature_base64"}
        row["signature_base64"] = base64.b64encode(
            context["private_keys"]["preregistration_authority"].sign(
                external.canonical_bytes(payload)
            )
        ).decode("ascii")
    registry = context["evidence"]["nonce_registry"]
    for assignment in registry["assignments"]:
        assignment["binding_sha256"] = context["bundle"]["bundle_binding"][
            "binding_sha256"
        ]
    registry_payload = {
        key: value
        for key, value in registry.items()
        if key != "signature_base64"
    }
    registry["signature_base64"] = base64.b64encode(
        context["private_keys"]["replay_authority"].sign(
            external.canonical_bytes(registry_payload)
        )
    ).decode("ascii")
    stored_registry = json.loads(
        context["registry_path"].read_text(encoding="utf-8")
    )
    stored_registry["assignments"] = registry["assignments"]
    context["registry_path"].write_text(
        json.dumps(stored_registry, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_external(context)


def verify_authenticated(context: dict) -> dict:
    return verify(
        context["bundle_path"],
        replay_registry=context["registry_path"],
        trust_policy=context["policy_path"],
        external_evidence=context["evidence_path"],
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


def test_authenticated_contract_has_a_nonvacuous_sufficient_path(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    report = verify_authenticated(context)
    assert report["verdict"] == "SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL"
    assert report["open_gates"] == []
    assert report["integrity_valid"] is True
    assert all(report["predicates"].values())
    assert "key-independent-attestor" in report["verified_external_anchor_ids"]
    cli = subprocess.run(
        [
            sys.executable,
            str(VERIFIER_PATH),
            str(context["bundle_path"]),
            "--replay-registry",
            str(context["registry_path"]),
            "--trust-policy",
            str(context["policy_path"]),
            "--external-evidence",
            str(context["evidence_path"]),
        ],
        capture_output=True,
        text=True,
    )
    assert cli.returncode == 0
    assert '"verdict": "SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL"' in cli.stdout


def test_bad_ed25519_signature_is_invalid(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    context["evidence"]["signatures"][0]["signature_base64"] = base64.b64encode(
        b"\x00" * 64
    ).decode("ascii")
    write_external(context)
    assert_rejected(
        verify_authenticated(context),
        "SIGNATURE_INVALID",
        verdict="INVALID",
    )


def test_witness_signature_cannot_be_reused_after_bundle_rebinding(
    tmp_path: Path,
) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    context["bundle"]["claim"]["conditions"] = "Rebound claimant-controlled conditions."
    write_bundle(context["bundle_path"], context["bundle"])
    claimant_row = next(
        item
        for item in context["evidence"]["signatures"]
        if item["subject_kind"] == "bundle_binding"
        and item["key_id"] == "key-claimant"
    )
    claimant_row["sha256"] = context["bundle"]["bundle_binding"]["binding_sha256"]
    payload = {
        key: value
        for key, value in claimant_row.items()
        if key != "signature_base64"
    }
    claimant_row["signature_base64"] = base64.b64encode(
        context["private_keys"]["claimant"].sign(external.canonical_bytes(payload))
    ).decode("ascii")
    write_external(context)
    assert_rejected(
        verify_authenticated(context),
        "SIGNATURE_SUBJECT_MISMATCH",
        verdict="INVALID",
    )


def test_self_authored_external_key_is_not_an_operator_trust_root(
    tmp_path: Path,
) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    row = next(
        item
        for item in context["evidence"]["signatures"]
        if item["subject_kind"] == "artifact"
        and artifact_row(context["bundle"], item["subject_id"])["role"]
        == "attestation"
    )
    row["key_id"] = "attacker-self-authored-key"
    row["signature_base64"] = base64.b64encode(b"\x00" * 64).decode("ascii")
    write_external(context)
    assert_rejected(
        verify_authenticated(context),
        "SIGNATURE_KEY_UNTRUSTED",
    )


def test_same_organization_witness_is_not_independent(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    claimant_org = next(
        row["organization_id"]
        for row in context["policy"]["anchors"]
        if row["role"] == "claimant"
    )
    next(
        row
        for row in context["policy"]["anchors"]
        if row["role"] == "independent_attestor"
    )["organization_id"] = claimant_org
    write_external(context)
    assert_rejected(
        verify_authenticated(context),
        "FRESH_REPRODUCTION_BUNDLE_VERIFICATION_OPEN",
    )


def test_post_capture_schedule_commitment_is_rejected(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    row = next(
        item
        for item in context["evidence"]["commitments"]
        if item["kind"] == "run_schedule"
    )
    row["committed_utc"] = "2026-07-20T12:05:00Z"
    payload = {key: value for key, value in row.items() if key != "signature_base64"}
    row["signature_base64"] = base64.b64encode(
        context["private_keys"]["preregistration_authority"].sign(
            external.canonical_bytes(payload)
        )
    ).decode("ascii")
    write_external(context)
    assert_rejected(verify_authenticated(context), "COMMITMENT_NOT_PRE_RUN")


def test_preregistration_commitment_cannot_be_reused_for_another_bundle(
    tmp_path: Path,
) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    context["bundle"]["bundle_id"] = "different-class-h-bundle"
    write_artifact(
        context["bundle_path"],
        context["bundle"],
        "attestation-001",
        {
            "bundle_id": context["bundle"]["bundle_id"],
            "mode": "independent_end_to_end_witness",
            "witnessed_scope": [
                "device identity",
                "capture",
                "controls",
                "custody",
            ],
        },
    )
    refresh_authenticated_signatures(context)
    assert_rejected(
        verify_authenticated(context),
        "COMMITMENT_SUBJECT_MISMATCH",
        verdict="INVALID",
    )


def test_commitment_key_must_be_valid_at_commitment_time(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    anchor = next(
        row
        for row in context["policy"]["anchors"]
        if row["role"] == "preregistration_authority"
    )
    anchor["valid_from_utc"] = "2026-07-20T11:00:00Z"
    write_external(context)
    assert_rejected(
        verify_authenticated(context),
        "COMMITMENT_KEY_NOT_VALID_AT_EVENT",
    )


def test_nonce_reserved_for_another_bundle_is_rejected(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    registry = context["evidence"]["nonce_registry"]
    registry["assignments"][0]["bundle_id"] = "different-bundle"
    payload = {key: value for key, value in registry.items() if key != "signature_base64"}
    registry["signature_base64"] = base64.b64encode(
        context["private_keys"]["replay_authority"].sign(
            external.canonical_bytes(payload)
        )
    ).decode("ascii")
    stored = json.loads(context["registry_path"].read_text(encoding="utf-8"))
    stored["assignments"] = registry["assignments"]
    context["registry_path"].write_text(
        json.dumps(stored, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_external(context)
    assert_rejected(
        verify_authenticated(context),
        "REPLAY_NONCE_NOT_ATOMICALLY_CONSUMED",
    )


def test_consumed_nonce_is_bound_to_the_canonical_bundle_root(
    tmp_path: Path,
) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    registry = context["evidence"]["nonce_registry"]
    registry["assignments"][0]["binding_sha256"] = "0" * 64
    payload = {key: value for key, value in registry.items() if key != "signature_base64"}
    registry["signature_base64"] = base64.b64encode(
        context["private_keys"]["replay_authority"].sign(
            external.canonical_bytes(payload)
        )
    ).decode("ascii")
    stored = json.loads(context["registry_path"].read_text(encoding="utf-8"))
    stored["assignments"] = registry["assignments"]
    context["registry_path"].write_text(
        json.dumps(stored, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_external(context)
    assert_rejected(
        verify_authenticated(context),
        "REPLAY_NONCE_NOT_ATOMICALLY_CONSUMED",
    )


def test_stale_registry_snapshot_is_rejected(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    registry = context["evidence"]["nonce_registry"]
    registry["snapshot_utc"] = "2026-07-20T11:00:00Z"
    payload = {key: value for key, value in registry.items() if key != "signature_base64"}
    registry["signature_base64"] = base64.b64encode(
        context["private_keys"]["replay_authority"].sign(
            external.canonical_bytes(payload)
        )
    ).decode("ascii")
    write_external(context)
    assert_rejected(
        verify_authenticated(context),
        "REPLAY_REGISTRY_SNAPSHOT_STALE",
    )


def test_unsigned_seen_nonce_state_drift_is_invalid(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    stored = json.loads(context["registry_path"].read_text(encoding="utf-8"))
    stored["seen_nonces"] = ["unrelated-prior-nonce"]
    context["registry_path"].write_text(
        json.dumps(stored, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    assert_rejected(
        verify_authenticated(context),
        "REPLAY_REGISTRY_SNAPSHOT_MISMATCH",
        verdict="INVALID",
    )


def test_signed_raw_mutation_cannot_bypass_analysis_replay(tmp_path: Path) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    write_artifact(
        context["bundle_path"],
        context["bundle"],
        "raw-001",
        {
            "capture_nonce": "nonce-synthetic-001",
            "channel": "synthetic_voltage",
            "device_id": "synthetic-device-001",
            "run_id": "run-synthetic-001",
            "samples": [1.0, 1.0, 1.0],
        },
    )
    refresh_authenticated_signatures(context)
    assert_rejected(
        verify_authenticated(context),
        "ANALYSIS_TO_CLAIM_EXECUTION_OPEN",
    )


def test_malformed_nested_external_input_returns_invalid_not_traceback(
    tmp_path: Path,
) -> None:
    context = make_authenticated_contract_fixture(tmp_path)
    context["policy"]["anchors"][0]["role"] = []
    write_external(context)
    report = verify_authenticated(context)
    assert_rejected(
        report,
        "EXTERNAL_EVIDENCE_VERIFICATION_EXCEPTION",
        verdict="INVALID",
    )
    cli = subprocess.run(
        [
            sys.executable,
            str(VERIFIER_PATH),
            str(context["bundle_path"]),
            "--replay-registry",
            str(context["registry_path"]),
            "--trust-policy",
            str(context["policy_path"]),
            "--external-evidence",
            str(context["evidence_path"]),
        ],
        capture_output=True,
        text=True,
    )
    assert cli.returncode == 2
    assert "Traceback" not in cli.stderr
