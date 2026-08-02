#!/usr/bin/env python3
"""Regression tests for the issue-646 electroweak pole quotient."""

from __future__ import annotations

import copy
import importlib.util
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import electroweak_pole_quotient_certificate as cert  # noqa: E402
import verify_electroweak_pole_quotient as independent  # noqa: E402


def _receipt() -> dict:
    return cert.build_receipt()


def _rehash(receipt: dict) -> None:
    payload = dict(receipt)
    payload.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = cert.tagged_sha256(
        cert.canonical_json_bytes(payload)
    )


def test_runtime_receipt_is_byte_exact() -> None:
    assert cert.RECEIPT_PATH.read_bytes() == cert.canonical_json_bytes(_receipt())
    cert.verify_runtime()


def test_global_scale_cancellation_is_exact() -> None:
    base = cert.exact_readout(
        common_scale_squared=Fraction(4),
        t_squared=Fraction(3, 5),
        w_mass_factor=Fraction(9, 8),
        w_width_over_mass=Fraction(1, 30),
        z_mass_factor=Fraction(11, 10),
        z_width_over_mass=Fraction(1, 40),
    )
    rescaled = cert.exact_readout(
        common_scale_squared=Fraction(49),
        t_squared=Fraction(3, 5),
        w_mass_factor=Fraction(9, 8),
        w_width_over_mass=Fraction(1, 30),
        z_mass_factor=Fraction(11, 10),
        z_width_over_mass=Fraction(1, 40),
    )
    assert base["absolute_poles"] != rescaled["absolute_poles"]
    assert base["dimensionless_outputs"] == rescaled["dimensionless_outputs"]


def test_factored_coordinates_match_the_pinned_strict_consumer() -> None:
    spec = importlib.util.spec_from_file_location(
        "strict_wz_consumer_for_quotient_test", cert.CONSUMER_SOURCE
    )
    assert spec is not None and spec.loader is not None
    consumer = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = consumer
    spec.loader.exec_module(consumer)

    g = 0.6
    gp = 0.3
    v_f = 2.0
    masses = consumer.tree_masses(g, gp, v_f)
    scale_squared = (g * v_f / 2.0) ** 2
    t_squared = (gp / g) ** 2
    assert math.isclose(masses.w, scale_squared, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(
        masses.z,
        scale_squared * (1.0 + t_squared),
        rel_tol=0.0,
        abs_tol=1e-15,
    )

    w_mass_factor = Fraction(11, 10)
    w_width_ratio = Fraction(1, 20)
    z_mass_factor = Fraction(9, 8)
    z_width_ratio = Fraction(1, 25)
    p_w = cert.pole_factor(w_mass_factor, w_width_ratio)
    p_z = cert.pole_factor(z_mass_factor, z_width_ratio)
    s_w = complex(float(p_w.real), float(p_w.imaginary)) * masses.w
    s_z = complex(float(p_z.real), float(p_z.imaginary)) * masses.z
    w_readout = consumer.energy_pole_coordinates(s_w)
    z_readout = consumer.energy_pole_coordinates(s_z)
    assert math.isclose(
        w_readout["M_GeV"] / z_readout["M_GeV"],
        float(w_mass_factor)
        / (math.sqrt(1.0 + t_squared) * float(z_mass_factor)),
        rel_tol=2e-15,
    )
    assert math.isclose(
        w_readout["Gamma_GeV"] / w_readout["M_GeV"],
        float(w_width_ratio),
        rel_tol=2e-15,
    )
    assert math.isclose(
        z_readout["Gamma_GeV"] / z_readout["M_GeV"],
        float(z_width_ratio),
        rel_tol=2e-15,
    )


def test_every_surviving_consumer_direction_has_an_exact_mutation() -> None:
    rows = {row["name"]: row for row in cert.counterfamilies()["mutations"]}
    assert rows["common_scale"]["changed_dimensionless_outputs"] == []
    assert rows["weak_coupling_ratio"]["changed_dimensionless_outputs"] == [
        "complex_pole_ratio_s_w_over_s_z",
        "mass_ratio_squared",
    ]
    assert "gamma_w_over_m_w" in rows["w_absorptive_pole_factor"][
        "changed_dimensionless_outputs"
    ]
    assert "gamma_z_over_m_z" in rows["z_absorptive_pole_factor"][
        "changed_dimensionless_outputs"
    ]


def test_output_vocabulary_keeps_absent_consumers_open() -> None:
    rows = {row["output"]: row for row in _receipt()["frozen_output_vocabulary"]}
    assert rows["pole-residue ratios"]["status"] == "not_emitted_by_strict_pole_consumer"
    assert rows["asymmetry combinations"]["status"] == "not_emitted_by_strict_pole_consumer"
    residue_inputs = rows["pole-residue ratios"]["consumer_inputs"]
    assert "charged inverse-propagator derivative at the W pole" in residue_inputs
    assert any("full neutral inverse/self-energy matrix" in item for item in residue_inputs)
    assert any("left/right null-vector" in item for item in residue_inputs)


def test_quotient_truncation_and_kernel_fields_are_explicit() -> None:
    receipt = _receipt()
    quotient = receipt["global_quotient"]
    assert "one-loop-truncated" in quotient["truncation_typing"]
    assert quotient["strict_one_loop_expansion"] == (
        "sW/sZ=(1/(1+t^2))*(1+dW-dZ)+O(loop^2)"
    )
    jacobian = receipt["symbolic_log_jacobian"]
    assert "dimension three over C" in jacobian["coordinate_field"]
    assert "four-real-dimensional kernel" in jacobian["physical_real_slice"]


def test_minimal_mass_width_contract_requires_no_absolute_scale() -> None:
    row = _receipt()["minimal_surviving_source_contract"]
    assert row["absolute_scale_required"] is False
    assert len(row["mass_width_triple"]) == 4
    assert any("dW" in item for item in row["mass_width_triple"])
    assert any("dZ" in item for item in row["mass_width_triple"])


def test_independent_verifier_accepts_exact_packet() -> None:
    result = independent.verify_receipt(_receipt())
    assert result == {
        "schema": "oph.electroweak_pole_scale_quotient.verification.v1",
        "receipt": True,
        "reasons": [],
        "exact_global_scale_quotient": True,
        "counterfamily_mutations": 6,
        "comparison_data_read": False,
    }


def test_independent_verifier_rejects_semantic_mutations_with_valid_hash() -> None:
    mutated = copy.deepcopy(_receipt())
    mutated["symbolic_log_jacobian"]["rows"]["d_log_sW_over_sZ"][0] = "1"
    _rehash(mutated)
    result = independent.verify_receipt(mutated)
    assert result["receipt"] is False
    assert "symbolic Jacobian drift" in result["reasons"]

    promoted = copy.deepcopy(_receipt())
    promoted["comparison_boundary"]["comparison_permitted"] = True
    _rehash(promoted)
    result = independent.verify_receipt(promoted)
    assert result["receipt"] is False
    assert "comparison boundary drift" in result["reasons"]


def test_independent_verifier_rejects_re_signed_promotion_and_parent_redirect() -> None:
    promoted = copy.deepcopy(_receipt())
    promoted["physical_prediction"] = True
    _rehash(promoted)
    result = independent.verify_receipt(promoted)
    assert result["receipt"] is False
    assert "top-level key inventory drift" in result["reasons"]

    redirected = copy.deepcopy(_receipt())
    arbitrary = cert.REPO_ROOT / "README.md"
    redirected["parent_pins"] = [
        {
            "path": "README.md",
            "bytes": len(arbitrary.read_bytes()),
            "sha256": cert.tagged_sha256(arbitrary.read_bytes()),
        }
    ] * 4
    _rehash(redirected)
    result = independent.verify_receipt(redirected)
    assert result["receipt"] is False
    assert "parent pin inventory drift" in result["reasons"]


def test_independent_verifier_rejects_claim_surface_and_lean_mutations() -> None:
    cases = (
        (
            lambda value: value["frozen_output_vocabulary"][3].update(
                {"status": "physically_emitted"}
            ),
            "output vocabulary drift",
        ),
        (
            lambda value: value["upstream_dependency_factorization"][
                "active_vacuum_scale"
            ].update({"dimensionless_output_status": "proved_canceled"}),
            "dependency factorization drift",
        ),
        (
            lambda value: value["parent_pins"][3].update({"theorems": []}),
            "parent pin inventory drift",
        ),
        (
            lambda value: value["counterfamilies"].update(
                {"typing": "physical source-admissible family"}
            ),
            "counterfamily typing drift",
        ),
    )
    for mutation, message in cases:
        changed = copy.deepcopy(_receipt())
        mutation(changed)
        _rehash(changed)
        result = independent.verify_receipt(changed)
        assert result["receipt"] is False
        assert message in result["reasons"]


def test_independent_verifier_rejects_counterfamily_overclaim() -> None:
    mutated = copy.deepcopy(_receipt())
    mutated["counterfamilies"]["mutations"][0][
        "changed_dimensionless_outputs"
    ] = ["mass_ratio_squared"]
    _rehash(mutated)
    result = independent.verify_receipt(mutated)
    assert result["receipt"] is False
    assert "mutation claim drift: common_scale" in result["reasons"]


def test_strict_json_loader_rejects_duplicate_and_nonfinite_values(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"a","schema":"b"}\n', encoding="ascii")
    try:
        independent.strict_load(duplicate)
    except independent.VerificationError as exc:
        assert "duplicate JSON key" in str(exc)
    else:
        raise AssertionError("duplicate key was accepted")

    nonfinite = tmp_path / "nonfinite.json"
    nonfinite.write_text('{"value":NaN}\n', encoding="ascii")
    try:
        independent.strict_load(nonfinite)
    except independent.VerificationError as exc:
        assert "non-finite JSON constant" in str(exc)
    else:
        raise AssertionError("non-finite value was accepted")

    pretty = tmp_path / "pretty.json"
    pretty.write_text(json.dumps(_receipt(), indent=2) + "\n", encoding="ascii")
    try:
        independent.strict_load(pretty)
    except independent.VerificationError as exc:
        assert "noncanonical JSON" in str(exc)
    else:
        raise AssertionError("noncanonical JSON was accepted")


def test_receipt_contains_no_float_or_measurement_payload() -> None:
    receipt = _receipt()

    def walk(value: object) -> None:
        if isinstance(value, float):
            raise AssertionError("receipt contains a float")
        if isinstance(value, dict):
            for child in value.values():
                walk(child)
        if isinstance(value, list):
            for child in value:
                walk(child)

    walk(receipt)
    serialized = json.dumps(receipt)
    assert "wz_pdg_2026_target_fixture" not in serialized
    assert receipt["comparison_boundary"]["public_measurement_read"] is False
