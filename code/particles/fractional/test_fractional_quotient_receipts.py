#!/usr/bin/env python3
"""Smoke-test the fractional quotient-sector receipt scaffold."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "fractional" / "build_fractional_quotient_receipts.py"


def test_build_fractional_quotient_receipts(tmp_path: pathlib.Path) -> None:
    out_dir = tmp_path / "fractional"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(out_dir)],
        check=True,
        cwd=ROOT,
    )

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["milestone"] == "FRACTIONAL_QUOTIENT_SECTOR_SANDBOX"
    assert manifest["strongest_allowed_claim"] == "FRACTIONAL_QUOTIENT_SANDBOX_DIAGNOSTIC"
    assert manifest["first_blocked_gate"] == "MATERIAL_SPECIFIC_HAMILTONIAN_PROOF_RECEIPT"
    assert manifest["promotion_allowed"] is False
    assert manifest["material_claim"] is False
    assert manifest["missing_files"] == []

    for rel_path in manifest["required_files"]:
        assert (out_dir / rel_path).is_file(), rel_path

    receipts = json.loads((out_dir / "receipts.json").read_text(encoding="utf-8"))
    gates = receipts["readiness_gates"]
    assert gates["SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT"] is True
    assert gates["OPTICAL_MODULE_CERTIFICATE"] is True
    assert gates["MATERIAL_SPECIFIC_HAMILTONIAN_PROOF_RECEIPT"] is False

    provenance = receipts["receipt_provenance"]
    for name in (
        "QUOTIENT_LUMPABILITY",
        "CANONICALIZER_IDEMPOTENCE",
        "REPRESENTATIVE_INVARIANCE",
        "NO_ORBIT_SIZE_BIAS",
    ):
        assert gates[name] is True
        assert provenance[name] == "COMPUTED_FROM_LEAN_CERTIFICATE"
    assert provenance["SOURCE_HAMILTONIAN_FROZEN"] == "DECLARED_SANDBOX_SCAFFOLD"

    binding = json.loads((out_dir / "lean_certificate.json").read_text(encoding="utf-8"))
    assert binding["certificate_status"] == "CERTIFICATE_PINNED"
    assert binding["certificate_verdicts"]["QUOTIENT_LUMPABILITY"] is True

    dag = json.loads((out_dir / "no_target_leak_dag.json").read_text(encoding="utf-8"))
    assert dag["status"] == "PASS_EMPTY_COMPARISON_DAG"


def test_rejects_optical_target_as_source_input(tmp_path: pathlib.Path) -> None:
    config = tmp_path / "source_config.json"
    config.write_text(json.dumps({"source_inputs": ["optical_peak_measurement"]}), encoding="utf-8")
    out_dir = tmp_path / "leaky"

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--output",
            str(out_dir),
            "--config",
            str(config),
        ],
        check=True,
        cwd=ROOT,
    )

    dag = json.loads((out_dir / "no_target_leak_dag.json").read_text(encoding="utf-8"))
    receipts = json.loads((out_dir / "receipts.json").read_text(encoding="utf-8"))
    assert dag["status"] == "FAIL_FORBIDDEN_SOURCE_INPUT"
    assert "optical_peak_measurement" in dag["target_leak_hits"]
    assert receipts["readiness_gates"]["NO_TARGET_LEAK_DAG"] is False
    assert receipts["readiness_gates"]["SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT"] is False


def test_negative_control_certificate_fails_receipts(tmp_path: pathlib.Path) -> None:
    negative = (
        ROOT
        / "particles"
        / "fractional"
        / "fractional_quotient_negative_control_certificate.json"
    )
    out_dir = tmp_path / "negative"

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--output",
            str(out_dir),
            "--certificate",
            str(negative),
        ],
        check=True,
        cwd=ROOT,
    )

    receipts = json.loads((out_dir / "receipts.json").read_text(encoding="utf-8"))
    gates = receipts["readiness_gates"]
    for name in (
        "QUOTIENT_LUMPABILITY",
        "CANONICALIZER_IDEMPOTENCE",
        "REPRESENTATIVE_INVARIANCE",
        "NO_ORBIT_SIZE_BIAS",
        "SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT",
    ):
        assert gates[name] is False, name

    binding = json.loads((out_dir / "lean_certificate.json").read_text(encoding="utf-8"))
    assert binding["certificate_status"] == "CERTIFICATE_PINNED"
    assert binding["certificate_witnesses"]["QUOTIENT_LUMPABILITY"] == {
        "pair": ["a0", "a1"],
        "fibre": "vacuum",
    }

    ladder = json.loads((out_dir / "claim_ladder.json").read_text(encoding="utf-8"))
    assert ladder["claim"] == "DIAGNOSTIC_ONLY"
    assert ladder["first_blocked_gate"] == "CANONICALIZER_IDEMPOTENCE"


def test_certificate_not_matching_lean_pin_fails_closed(tmp_path: pathlib.Path) -> None:
    pinned = ROOT / "particles" / "fractional" / "fractional_quotient_certificate.json"
    tampered = json.loads(pinned.read_text(encoding="utf-8"))
    tampered["kernel"]["a0"]["vac"] = "9/10"
    cert_path = tmp_path / "tampered.json"
    cert_path.write_text(json.dumps(tampered), encoding="utf-8")
    out_dir = tmp_path / "tampered"

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--output",
            str(out_dir),
            "--certificate",
            str(cert_path),
        ],
        check=True,
        cwd=ROOT,
    )

    receipts = json.loads((out_dir / "receipts.json").read_text(encoding="utf-8"))
    gates = receipts["readiness_gates"]
    for name in (
        "QUOTIENT_LUMPABILITY",
        "CANONICALIZER_IDEMPOTENCE",
        "REPRESENTATIVE_INVARIANCE",
        "NO_ORBIT_SIZE_BIAS",
    ):
        assert gates[name] is False, name

    binding = json.loads((out_dir / "lean_certificate.json").read_text(encoding="utf-8"))
    assert binding["certificate_status"] == "CERTIFICATE_NOT_PINNED"


def test_runtime_kernel_harness_in_sync() -> None:
    """The checked-in Lean harness must be the exact render of the capture."""
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "particles" / "fractional" / "generate_runtime_kernel_harness.py"),
            "--check",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_runtime_kernel_harness_check_detects_drift(tmp_path: pathlib.Path) -> None:
    harness = (
        ROOT.parent
        / "Lean"
        / "ObserverPatchHolography"
        / "QuotientLumpabilityRuntimeHarness.lean"
    )
    original = harness.read_text(encoding="utf-8")
    try:
        harness.write_text(original.replace("1 / 2", "1 / 3", 1), encoding="utf-8")
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "particles" / "fractional" / "generate_runtime_kernel_harness.py"),
                "--check",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 1
        assert "STALE" in result.stdout
    finally:
        harness.write_text(original, encoding="utf-8")
