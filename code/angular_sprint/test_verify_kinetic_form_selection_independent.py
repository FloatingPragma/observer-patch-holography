#!/usr/bin/env python3
"""Adversarial tests for the independent issue-646 verifier."""

from __future__ import annotations

import copy
import hashlib
import json
import shutil
import sys
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import verify_kinetic_form_selection_independent as independent  # noqa: E402


def _receipt() -> dict:
    return independent.strict_load(independent.DEFAULT_RECEIPT)


def _rehash(receipt: dict) -> None:
    payload = dict(receipt)
    payload.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = independent.tagged_sha256(
        independent.canonical_json_bytes(payload)
    )


def _copy_verifier_inputs(tmp_path: Path) -> None:
    for relative in (
        independent.PORT_RECEIPT_REL,
        independent.CARRIER_MANIFEST_REL,
        independent.MATTER_ATTACHMENT_REL,
        Path("Lean/Screen/RGRepresentationFrontier.lean"),
    ):
        destination = tmp_path / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(independent.REPO_ROOT / relative, destination)


def _repin(receipt: dict, repo_root: Path, relative: Path) -> None:
    payload = (repo_root / relative).read_bytes()
    row = next(
        pin
        for pin in receipt["parent_pins"]
        if pin["path"] == relative.as_posix()
    )
    row["bytes"] = len(payload)
    row["sha256"] = "sha256:" + hashlib.sha256(payload).hexdigest()
    _rehash(receipt)


def test_independent_verifier_reconstructs_the_complete_exact_core() -> None:
    result = independent.verify_receipt(_receipt())
    assert result["receipt"] is True
    assert result["exact_core"] == {
        "jacobi_checks": 1728,
        "ad_invariance_checks": 1728,
        "su2": "1",
        "su3": "1/6",
        "ratio": "6",
    }
    assert result["general_family_higgs_cancellation"]["nG_cancels"] is True
    assert result["registered_matter_indices_u1_su2_su3"] == [
        "10/3",
        "2",
        "2",
    ]
    assert result["comparison_data_read"] is False


def test_independent_verifier_rejects_semantic_mutations_with_valid_self_hash() -> None:
    receipt = _receipt()
    mutated = copy.deepcopy(receipt)
    mutated["matter_trace_branch"]["frozen_rg_statistic"][
        "general_family_higgs_cancellation"
    ]["determinant_cofactors_constant_nG_nH"][0][1] = "1"
    _rehash(mutated)
    result = independent.verify_receipt(mutated)
    assert result["receipt"] is False
    assert "general cofactor certificate mismatch" in result["reasons"]

    promoted = copy.deepcopy(receipt)
    promoted["physical_selection_boundary"]["physical_sector_selected"] = True
    _rehash(promoted)
    result = independent.verify_receipt(promoted)
    assert result["receipt"] is False
    assert "forbidden physical promotion" in result["reasons"][0]

    overstated = copy.deepcopy(receipt)
    overstated["verification"]["classification"] = (
        "independently derived the structure constants from generators"
    )
    _rehash(overstated)
    result = independent.verify_receipt(overstated)
    assert result["receipt"] is False
    assert "verification scope" in result["reasons"][0]


def test_independent_verifier_rejects_parent_pin_and_payload_mutations() -> None:
    receipt = _receipt()
    mutated = copy.deepcopy(receipt)
    mutated["parent_pins"][0]["sha256"] = "sha256:" + "0" * 64
    _rehash(mutated)
    result = independent.verify_receipt(mutated)
    assert result["receipt"] is False
    assert "parent hash drift" in result["reasons"][0]

    stale_hash = copy.deepcopy(receipt)
    stale_hash["ad_invariance"]["verified_basis_triples"] = 1
    result = independent.verify_receipt(stale_hash)
    assert result["receipt"] is False
    assert "receipt self-hash mismatch" in result["reasons"]


def test_independent_verifier_rejects_unbound_carrier_mutation(
    tmp_path: Path,
) -> None:
    _copy_verifier_inputs(tmp_path)
    carrier_path = tmp_path / independent.CARRIER_MANIFEST_REL
    carrier = independent.strict_load(carrier_path)
    carrier["carrier"]["ports"][0] = "mutated-port"
    carrier_path.write_bytes(independent.canonical_json_bytes(carrier))

    result = independent.verify_receipt(_receipt(), repo_root=tmp_path)
    assert result["receipt"] is False
    assert "carrier manifest is not bound" in result["reasons"][0]


def test_independent_verifier_rejects_rehashed_structure_constant_mutation(
    tmp_path: Path,
) -> None:
    _copy_verifier_inputs(tmp_path)
    parent_path = tmp_path / independent.PORT_RECEIPT_REL
    parent = independent.strict_load(parent_path)
    parent["closure"]["structure_constants"]["[0,1]"][0] = "1"
    parent_path.write_bytes(independent.canonical_json_bytes(parent))

    receipt = _receipt()
    _repin(receipt, tmp_path, independent.PORT_RECEIPT_REL)
    result = independent.verify_receipt(receipt, repo_root=tmp_path)
    assert result["receipt"] is False
    assert any(
        marker in result["reasons"][0]
        for marker in ("Jacobi failure", "ad-invariance failure")
    )


def test_strict_loader_rejects_duplicate_keys(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema":"a","schema":"b"}\n', encoding="ascii")
    try:
        independent.strict_load(path)
    except independent.IndependentVerificationError as exc:
        assert "duplicate JSON key" in str(exc)
    else:
        raise AssertionError("duplicate JSON key was accepted")


def test_strict_loader_rejects_nonfinite_json_constants(tmp_path: Path) -> None:
    for token in ("NaN", "Infinity", "-Infinity"):
        path = tmp_path / f"nonfinite-{token.replace('-', 'minus')}.json"
        path.write_text('{"value":' + token + '}\n', encoding="ascii")
        try:
            independent.strict_load(path)
        except independent.IndependentVerificationError as exc:
            assert "non-finite JSON constant" in str(exc)
        else:
            raise AssertionError(f"non-finite JSON constant {token} was accepted")


def test_independent_cli_accepts_the_committed_receipt() -> None:
    import subprocess

    completed = subprocess.run(
        [
            sys.executable,
            str(MODULE_DIR / "verify_kinetic_form_selection_independent.py"),
            "--receipt",
            str(independent.DEFAULT_RECEIPT),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    result = json.loads(completed.stdout)
    assert result["receipt"] is True
    assert result["exact_core"]["ad_invariance_checks"] == 1728
