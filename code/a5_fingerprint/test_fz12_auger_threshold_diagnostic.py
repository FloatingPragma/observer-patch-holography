from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import fz12_auger_threshold_diagnostic as diagnostic


VERIFIER = diagnostic.HERE / "verify_fz12_auger_threshold_diagnostic_independent.py"


def write_canonical(path: Path, value: dict) -> None:
    path.write_bytes(diagnostic.base.canonical_json_bytes(value))


def resign(value: dict, field: str) -> None:
    body = {key: item for key, item in value.items() if key != field}
    value[field] = diagnostic.base.tagged_sha256(
        diagnostic.base.canonical_json_bytes(body)
    )


def independent_run(
    *,
    source: Path = diagnostic.SOURCE_PATH,
    custody: Path = diagnostic.CUSTODY_PATH,
    observation: Path = diagnostic.OBSERVATION_MAP_PATH,
    receipt: Path = diagnostic.RECEIPT_PATH,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [
            sys.executable,
            str(VERIFIER),
            "--source",
            str(source),
            "--custody",
            str(custody),
            "--observation-map",
            str(observation),
            "--receipt",
            str(receipt),
        ],
        cwd=diagnostic.HERE,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def mutated_receipt(tmp_path: Path, mutate) -> Path:
    receipt = json.loads(diagnostic.RECEIPT_PATH.read_text())
    mutate(receipt)
    resign(receipt, "receipt_sha256")
    path = tmp_path / "fz12_auger_threshold_diagnostic_receipt.json"
    write_canonical(path, receipt)
    return path


def test_exact_basis_translation() -> None:
    receipt = diagnostic.build_receipt()
    assert receipt["basis_translation"]["exact_coefficients"] == {
        "c_over_a2": "-1/20",
        "d_iso_over_a4": "1/840",
        "d_I6_over_a4": "-1/12600",
        "delta_gamma2_over_a2": "-1/20",
        "delta_gamma4_iso_over_a4": "-2/525",
        "delta_gamma4_I6_over_a4": "-1/12600",
    }
    assert receipt["basis_translation"]["paper_consumes_delta_gamma4"] is False
    assert receipt["basis_translation"]["higher_order_remainder_control"] is False


def test_exposed_auger_scenario_is_typed_exactly() -> None:
    exposed = diagnostic.build_receipt()["exposed_public_input"]
    assert exposed["exposure_class"] == "EXPOSED_RETROSPECTIVE"
    assert exposed["doi"] == "10.1088/1475-7516/2022/01/023"
    assert exposed["arxiv"] == "2112.06773"
    assert exposed["selected_scenario_maximum_proton_energy_eV"] == "1e20"
    assert exposed["delta_gamma2_strict_lower_bound_eV_minus2"] == "-1e-58"
    assert exposed["confidence_level_for_direct_electromagnetic_bound"] is None
    assert "no confidence level" in exposed["confidence_statement"]
    assert "no electromagnetic bound" in exposed["reference_scenarios"]
    assert exposed["scenario_dependence"] is True


def test_strict_scale_bound_and_codata_conversion() -> None:
    bound = diagnostic.build_receipt()["conditional_bound"]
    assert bound["a_squared_strict_upper_bound_eV_minus2"] == "20e-58"
    assert bound["a_strict_upper_bound_exact_eV_inverse"] == "sqrt(20)*1e-29"
    assert bound["a_upper_approx_eV_inverse"] == "4.472135955e-29"
    assert bound["a_upper_approx_m"] == "8.824730839e-36"
    assert bound["a_upper_approx_planck_lengths"] == "0.5459986722"
    assert bound["headline_m"] == "8.8247e-36"
    assert bound["headline_planck_lengths"] == "0.5460"
    assert bound["unit_conversion"] == {
        "hbar_c_eV_m": "1.973269804e-7",
        "planck_length_m": "1.616255e-35",
        "constant_source": "CODATA 2022 nominal decimal values",
        "constant_uncertainties_propagated": False,
        "rounding_note": (
            "the published coefficient limit has one-significant-digit precision; "
            "additional conversion digits are reproducibility aids"
        ),
    }


def test_physical_attachments_and_scoring_stay_open() -> None:
    receipt = diagnostic.build_receipt()
    attachments = receipt["open_physical_attachments"]
    for key, value in attachments.items():
        if key != "statement":
            assert value is False
    assert all(
        word in attachments["statement"] for word in ("photon", "electron", "positron")
    )
    scope = receipt["scope_boundary"]
    assert scope["public_measurement_read"] is True
    assert scope["leading_EFT_mapping_only"] is True
    for key in (
        "comparison_permitted",
        "comparison_budget_consumed",
        "score_emitted",
        "evidence_claimed",
        "verdict_emitted",
        "OPH_exclusion",
        "frozen_FZ12_modified",
    ):
        assert scope[key] is False
    assert "cannot confirm, support, falsify, or exclude OPH" in scope["statement"]


def test_three_exact_parent_pins_and_no_parent_writes() -> None:
    before = [contract["path"].read_bytes() for contract in diagnostic.PARENT_CONTRACTS]
    pins = diagnostic.build_receipt()["parent_pins"]
    after = [contract["path"].read_bytes() for contract in diagnostic.PARENT_CONTRACTS]
    assert before == after
    assert len(pins) == 3
    assert [pin["bytes"] for pin in pins] == [9296, 3624, 5276]
    assert [pin["sha256"] for pin in pins] == [
        "sha256:0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915",
        "sha256:dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643",
        "sha256:0b34218b075a3d51fe0badb00e3bb889743ae19b6dcf161f7877683e25121d17",
    ]


def test_committed_receipt_is_canonical_byte_exact() -> None:
    assert diagnostic.verify_committed_receipt() == diagnostic.build_receipt()
    assert diagnostic.RECEIPT_PATH.read_bytes() == diagnostic.base.canonical_json_bytes(
        diagnostic.build_receipt()
    )


def test_independent_verifier_is_separate_and_passes() -> None:
    source = VERIFIER.read_text(encoding="utf-8")
    assert "import fz12_auger_threshold_diagnostic" not in source
    assert "from fz12_auger_threshold_diagnostic" not in source
    result = independent_run()
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FZ12_AUGER_DIAGNOSTIC_VERIFY_PASS"


@pytest.mark.parametrize("parent_index", [0, 1, 2])
def test_resigned_parent_mutation_fails_raw_hash_pins(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, parent_index: int
) -> None:
    contract = diagnostic.PARENT_CONTRACTS[parent_index]
    parent = json.loads(contract["path"].read_text())
    parent["audit_mutation"] = True
    resign(parent, contract["self_field"])
    path = tmp_path / contract["path"].name
    write_canonical(path, parent)
    monkeypatch.setitem(contract, "path", path)
    with pytest.raises(diagnostic.base.FingerprintError, match="raw hash drift"):
        diagnostic.build_receipt()

    kwargs = {}
    if parent_index == 0:
        kwargs["source"] = path
    elif parent_index == 1:
        kwargs["custody"] = path
    else:
        kwargs["observation"] = path
    result = independent_run(**kwargs)
    assert result.returncode == 1
    assert "parent raw hash drift" in result.stderr


@pytest.mark.parametrize(
    "mutation",
    [
        "basis_coefficient",
        "energy_basis_map",
        "published_bound",
        "confidence_level",
        "reference_scenario",
        "photon_attachment",
        "interaction_attachment",
        "comparison_permitted",
        "comparison_budget",
        "score",
        "evidence",
        "verdict",
        "oph_exclusion",
        "leading_eft",
        "extra_key",
    ],
)
def test_resigned_semantic_mutations_fail_closed(tmp_path: Path, mutation: str) -> None:
    def mutate(receipt: dict) -> None:
        if mutation == "basis_coefficient":
            receipt["theory_input"]["coefficients"]["c"] = "-a^2/21"
        elif mutation == "energy_basis_map":
            receipt["basis_translation"]["exact_coefficients"][
                "delta_gamma4_iso_over_a4"
            ] = "1/840"
        elif mutation == "published_bound":
            receipt["exposed_public_input"][
                "delta_gamma2_strict_lower_bound_eV_minus2"
            ] = "-1e-59"
        elif mutation == "confidence_level":
            receipt["exposed_public_input"][
                "confidence_level_for_direct_electromagnetic_bound"
            ] = 0.95
        elif mutation == "reference_scenario":
            receipt["exposed_public_input"]["reference_scenarios"] = "bound"
        elif mutation == "photon_attachment":
            receipt["open_physical_attachments"]["photon_sector_attachment_proved"] = (
                True
            )
        elif mutation == "interaction_attachment":
            receipt["open_physical_attachments"][
                "shared_electromagnetic_interaction_kinematics_proved"
            ] = True
        elif mutation == "comparison_permitted":
            receipt["scope_boundary"]["comparison_permitted"] = True
        elif mutation == "comparison_budget":
            receipt["scope_boundary"]["comparison_budget_consumed"] = True
        elif mutation == "score":
            receipt["scope_boundary"]["score_emitted"] = True
        elif mutation == "evidence":
            receipt["scope_boundary"]["evidence_claimed"] = True
        elif mutation == "verdict":
            receipt["scope_boundary"]["verdict_emitted"] = True
        elif mutation == "oph_exclusion":
            receipt["scope_boundary"]["OPH_exclusion"] = True
        elif mutation == "leading_eft":
            receipt["scope_boundary"]["leading_EFT_mapping_only"] = False
        else:
            receipt["scope_boundary"]["unknown_promotion"] = True

    path = mutated_receipt(tmp_path, mutate)
    result = independent_run(receipt=path)
    assert result.returncode == 1
    assert "diagnostic receipt raw hash drift" in result.stderr


def test_producer_verifier_rejects_resigned_receipt_mutation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    path = mutated_receipt(
        tmp_path,
        lambda receipt: receipt["scope_boundary"].update({"score_emitted": True}),
    )
    monkeypatch.setattr(diagnostic, "RECEIPT_PATH", path)
    with pytest.raises(
        diagnostic.base.FingerprintError,
        match="ancestry, arithmetic, schema, or boundary drift",
    ):
        diagnostic.verify_committed_receipt()
