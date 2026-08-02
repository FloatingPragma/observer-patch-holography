from __future__ import annotations

import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

import a5_multipole_fixed_point_certificate as base
import fz12_joint_threshold_certificate as producer
import verify_fz12_joint_threshold_independent as independent


def live_receipt() -> dict:
    return json.loads(producer.RECEIPT_PATH.read_bytes())


def resign(value: dict) -> dict:
    changed = copy.deepcopy(value)
    changed.pop("receipt_sha256", None)
    changed["receipt_sha256"] = base.tagged_sha256(
        base.canonical_json_bytes(changed)
    )
    return changed


def write_receipt(path: Path, value: dict) -> None:
    path.write_bytes(base.canonical_json_bytes(value))


def test_joint_threshold_producer_replays_committed_bytes() -> None:
    assert base.canonical_json_bytes(producer.build_receipt()) == producer.RECEIPT_PATH.read_bytes()
    assert producer.verify_committed_receipt()["status"] == producer.STATUS


def test_joint_threshold_independent_verifier_replays() -> None:
    assert independent.verify(producer.RECEIPT_PATH)["status"] == producer.STATUS


def test_joint_threshold_cli_replays() -> None:
    commands = (
        [sys.executable, str(Path(producer.__file__)), "--verify"],
        [sys.executable, str(Path(independent.__file__))],
    )
    for command in commands:
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        assert completed.returncode == 0, completed.stdout + completed.stderr
        assert producer.STATUS in completed.stdout


@pytest.mark.parametrize(
    ("mutator", "message"),
    (
        (
            lambda value: value["exact_leading_derivation"].__setitem__(
                "equal_share_linear_rank", 3
            ),
            "rank drift",
        ),
        (
            lambda value: value["exact_leading_derivation"].__setitem__(
                "equal_share_coefficient_fiber_dimension", 0
            ),
            "fiber dimension drift",
        ),
        (
            lambda value: value["exact_photon_symbol_control"].__setitem__(
                "P8_remainder", "0 <= lambda_hat-P8 <= q^8"
            ),
            "photon-symbol control drift",
        ),
        (
            lambda value: value["exact_li_lepton_consequences"].__setitem__(
                "EFT_truncation_used", True
            ),
            "exact Lorentz-invariant-lepton boundary drift",
        ),
        (
            lambda value: value["exposure_boundary"].__setitem__(
                "existing_Auger_limit_read", True
            ),
            "exposure boundary drift",
        ),
        (
            lambda value: value["fz12_conditional_branches"][
                "universal_principal_symbol_premise"
            ].__setitem__("source_selected", True),
            "conditional-branch contract drift",
        ),
        (
            lambda value: value["parent_pins"][1].__setitem__(
                "sha256", "sha256:" + "0" * 64
            ),
            "parent pin inventory drift",
        ),
        (
            lambda value: value["kinematic_scope"].__setitem__(
                "general_independent_lepton_share_optimization_proved", True
            ),
            "kinematic-scope promotion drift",
        ),
        (
            lambda value: value["exact_leading_derivation"][
                "standard_lepton_share_optimization"
            ].__setitem__(
                "threshold_residual_result",
                "the maximizer depends on the photon coefficient",
            ),
            "standard-lepton share optimization drift",
        ),
        (
            lambda value: value["fz12_conditional_branches"][
                "standard_lepton_premise"
            ].__setitem__("source_selected", True),
            "conditional-branch contract drift",
        ),
        (
            lambda value: value["open_physical_attachments"].__setitem__(
                "electron_action_selected", True
            ),
            "open-attachment promotion drift",
        ),
        (
            lambda value: value.__setitem__("physical_prediction", True),
            "receipt top-level key drift",
        ),
        (
            lambda value: value["parent_pins"][3].update(
                {
                    "path": "Lean/Screen/A2HolonomyBridge.lean",
                    "bytes": 0,
                    "sha256": "sha256:" + "0" * 64,
                    "theorems": [],
                }
            ),
            "parent pin inventory drift",
        ),
    ),
)
def test_resigned_semantic_mutations_fail_closed(
    tmp_path: Path, mutator, message: str
) -> None:
    changed = live_receipt()
    mutator(changed)
    path = tmp_path / "mutated.json"
    write_receipt(path, resign(changed))
    with pytest.raises(independent.VerificationError, match=message):
        independent.verify(path)


def test_noncanonical_receipt_fails_closed(tmp_path: Path) -> None:
    path = tmp_path / "pretty.json"
    path.write_text(json.dumps(live_receipt(), indent=2) + "\n")
    with pytest.raises(independent.VerificationError, match="noncanonical"):
        independent.verify(path)


def test_duplicate_and_nonfinite_json_fail_closed(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"a","schema":"b"}\n', encoding="ascii")
    with pytest.raises(independent.VerificationError, match="duplicate JSON key"):
        independent.load_json(duplicate)

    nonfinite = tmp_path / "nonfinite.json"
    nonfinite.write_text('{"value":NaN}\n', encoding="ascii")
    with pytest.raises(independent.VerificationError, match="non-finite JSON"):
        independent.load_json(nonfinite)


@pytest.mark.parametrize(
    ("module", "error_type"),
    (
        (producer, base.FingerprintError),
        (independent, independent.VerificationError),
    ),
)
def test_lean_kernel_build_gate_rejects_failed_build(
    monkeypatch: pytest.MonkeyPatch, module, error_type
) -> None:
    class FailedBuild:
        returncode = 1
        stdout = ""
        stderr = "Lean elaboration failed"

    module.verify_lean_kernel_build.cache_clear()
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: FailedBuild())
    try:
        with pytest.raises(error_type, match="kernel build failed"):
            module.verify_lean_kernel_build()
    finally:
        module.verify_lean_kernel_build.cache_clear()


@pytest.mark.parametrize(
    ("module", "error_type"),
    (
        (producer, base.FingerprintError),
        (independent, independent.VerificationError),
    ),
)
def test_lean_kernel_build_gate_rejects_sorry_axiom(
    monkeypatch: pytest.MonkeyPatch, module, error_type
) -> None:
    class SorryBuild:
        returncode = 0
        stdout = "sorryAx"
        stderr = ""

    module.verify_lean_kernel_build.cache_clear()
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: SorryBuild())
    try:
        with pytest.raises(error_type, match="contains sorryAx"):
            module.verify_lean_kernel_build()
    finally:
        module.verify_lean_kernel_build.cache_clear()
