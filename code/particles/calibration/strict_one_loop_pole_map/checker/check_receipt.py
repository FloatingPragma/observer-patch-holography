#!/usr/bin/env python3
"""Small independent checker for a strict-one-loop W/Z pole-map receipt.

The checker intentionally does not import src/wz_pole_map.py.  It recomputes
all elementary equations from the receipt and fixture, validates the schema,
and applies fail-closed promotion logic.
"""

from __future__ import annotations

import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import jsonschema


REQUIRED_PROMOTION_GATES = [
    "full_neutral_matrix_supplied",
    "renormalized_FJ_vev_selected",
    "complete_FJ_tadpole_conversion_checked",
    "target_clean_source_matching",
    "source_law_or_deterministic_exactness_receipt",
    "independent_self_energy_engine",
    "finite_order_BRST_ST_Ward_Nielsen_receipt",
    "physical_clock_closed",
    "no_target_ancestry",
]

CONVENTION_ID = "Gamma=s-m2-Delta; s_pole=m2+Delta_at_tree"
SELF_ENERGY_TRANSLATION = "Delta=-Pi; Gamma=s-m2+Pi=s-m2-Delta"
CONDITIONAL_STATUS = "CONDITIONAL_STRICT_1L_POLE_MAP_NOT_OPH_NATIVE_PHYSICAL"
PROMOTION_POLICY = (
    "v1 is an algebraic diagnostic and cannot self-promote; a separate "
    "aggregate checker must verify hash-bound evidence artifacts and subject binding"
)
POLE_CONVENTION = "s_V=(M_V-i*Gamma_V/2)^2; Re(sqrt(s_V))>0; Im(sqrt(s_V))<=0"
COORDINATE_NOTE = (
    "Exact square root of the one-loop-truncated s pole; the nonlinear "
    "coordinate readout contains kinematic terms beyond strict one loop."
)
STRICT_NEUTRAL_RULE = "s_Z=z+Delta_ZZ^(1)(z); off-diagonal product excluded as O(h^2)"
NEUTRAL_EFFECTIVE_TERM = "-Delta_ZA^(1)(s)*Delta_AZ^(1)(s)/s"
MAX_ABS_GEV2_TOL = 1.0e-6
MAX_ABS_GEV_TOL = 1.0e-8
MAX_REL_TOL = 1.0e-10


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def c(pair: list[float]) -> complex:
    return complex(float(pair[0]), float(pair[1]))


def close(a: float, b: float, abs_tol: float, rel_tol: float) -> bool:
    return abs(a - b) <= max(abs_tol, rel_tol * max(1.0, abs(a), abs(b)))


def close_complex(a: complex, b: complex, abs_tol: float, rel_tol: float) -> bool:
    return close(a.real, b.real, abs_tol, rel_tol) and close(a.imag, b.imag, abs_tol, rel_tol)


def physical_root(s: complex) -> complex:
    root = cmath.sqrt(s)
    if root.real < 0.0 or (root.real == 0.0 and root.imag > 0.0):
        root = -root
    return root


def check(receipt_path: Path, fixture_path: Path, schema_path: Path) -> dict[str, Any]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    fixture_schema_path = schema_path.with_name("physical_wz_strict_1l_input_fixture.schema.json")
    fixture_schema = json.loads(fixture_schema_path.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(fixture_schema).validate(fixture)
    jsonschema.Draft202012Validator(schema).validate(receipt)

    assertions: dict[str, bool] = {}

    # Bind the recomputed subject to the exact fixture, not merely to a hash
    # string copied into an otherwise unrelated receipt.
    expected_source = dict(fixture.get("source", {}))
    expected_source["fixture_path"] = fixture_path.name
    expected_source["fixture_sha256"] = sha256(fixture_path)
    assertions["receipt_id_binding"] = receipt["receipt_id"] == fixture["receipt_id"]
    assertions["claim_class_binding"] = receipt["claim_class"] == fixture["claim_class"]
    assertions["source_binding"] = receipt["source"] == expected_source
    assertions["parameter_binding"] = (
        receipt["renormalized_parameters"] == fixture["renormalized_parameters"]
    )
    assertions["contribution_mask_binding"] = (
        receipt["contribution_mask"] == fixture["contribution_mask"]
    )
    assertions["source_covariance_binding"] = (
        receipt["source_covariance"] == fixture.get("source_covariance", {"status": "not_supplied"})
    )
    assertions["tolerance_binding"] = receipt["tolerances"] == fixture["tolerances"]
    assertions["inverse_convention_binding"] = (
        fixture.get("inverse_propagator_convention") == CONVENTION_ID
        and receipt["inverse_propagator_convention"] == CONVENTION_ID
    )
    assertions["self_energy_sign_translation"] = (
        receipt["self_energy_sign_translation"] == SELF_ENERGY_TRANSLATION
    )
    assertions["pole_convention_binding"] = receipt["pole_convention"] == POLE_CONVENTION
    assertions["fixture_order_contract"] = (
        fixture.get("perturbative_order") == "strict_one_loop"
        and fixture.get("loop_factor_included_in_delta") is True
        and fixture.get("contribution_mask", {}).get("maximum_loop_power") == 1
        and fixture.get("contribution_mask", {}).get("one_loop_complete_for_declared_backend") is True
    )

    tol2 = float(fixture["tolerances"]["absolute_GeV2"])
    tol = float(fixture["tolerances"]["absolute_GeV"])
    rel = float(fixture["tolerances"]["relative"])
    assertions["tolerance_safety_cap"] = (
        0.0 < tol2 <= MAX_ABS_GEV2_TOL
        and 0.0 < tol <= MAX_ABS_GEV_TOL
        and 0.0 < rel <= MAX_REL_TOL
    )
    p = receipt["renormalized_parameters"]
    g, gp, v = float(p["g"]), float(p["gp"]), float(p["v_F_GeV"])
    w = g * g * v * v / 4.0
    z = (g * g + gp * gp) * v * v / 4.0
    tree = receipt["tree"]
    assertions["tree_w"] = close(w, float(tree["w_GeV2"]), tol2, rel)
    assertions["tree_z"] = close(z, float(tree["z_GeV2"]), tol2, rel)
    assertions["tree_mW0"] = close(math.sqrt(w), float(tree["mW0_GeV"]), tol, rel)
    assertions["tree_mZ0"] = close(math.sqrt(z), float(tree["mZ0_GeV"]), tol, rel)

    dW = c(receipt["one_loop_coefficients"]["Delta_WW_at_w_GeV2"])
    dZ = c(receipt["one_loop_coefficients"]["Delta_ZZ_at_z_GeV2"])
    fixture_self_energy = fixture["renormalized_one_loop_self_energy"]
    assertions["Delta_WW_binding"] = close_complex(
        dW, c(fixture_self_energy["Delta_WW_at_w_GeV2"]), 0.0, 0.0
    )
    assertions["Delta_ZZ_binding"] = close_complex(
        dZ, c(fixture_self_energy["Delta_ZZ_at_z_GeV2"]), 0.0, 0.0
    )
    sW = complex(w, 0.0) + dW
    sZ = complex(z, 0.0) + dZ
    got_sW = c(receipt["complex_poles"]["W"]["s_pole_truncated_1L_GeV2"])
    got_sZ = c(receipt["complex_poles"]["Z"]["s_pole_truncated_1L_GeV2"])
    assertions["sW"] = close_complex(sW, got_sW, tol2, rel)
    assertions["sZ"] = close_complex(sZ, got_sZ, tol2, rel)

    for name, s, d, m0 in (
        ("W", sW, dW, math.sqrt(w)),
        ("Z", sZ, dZ, math.sqrt(z)),
    ):
        root = physical_root(s)
        readout = receipt["complex_poles"][name]["energy_pole_readout"]
        assertions[f"{name}_sqrt_M"] = close(root.real, float(readout["M_GeV"]), tol, rel)
        assertions[f"{name}_sqrt_Gamma"] = close(-2.0 * root.imag, float(readout["Gamma_GeV"]), tol, rel)
        assertions[f"{name}_sqrt_field"] = close_complex(
            root, c(readout["sqrt_s_GeV"]), tol, rel
        )
        reconstructed = complex(float(readout["M_GeV"]), -0.5 * float(readout["Gamma_GeV"])) ** 2
        assertions[f"{name}_root_reconstruction"] = close_complex(s, reconstructed, tol2, rel)
        assertions[f"{name}_reconstructed_field"] = close_complex(
            s, c(readout["s_reconstructed_GeV2"]), tol2, rel
        )
        assertions[f"{name}_coordinate_note"] = (
            readout["coordinate_transform_note"] == COORDINATE_NOTE
        )
        strict = receipt["one_loop_coefficients"][f"strict_{name}"]
        expected_delta_m = d.real / (2.0 * m0)
        expected_width = -d.imag / m0
        assertions[f"{name}_strict_M0"] = close(m0, float(strict["M0_GeV"]), tol, rel)
        assertions[f"{name}_strict_deltaM"] = close(d.real / (2.0 * m0), float(strict["delta_M_1L_GeV"]), tol, rel)
        assertions[f"{name}_strict_M"] = close(
            m0 + expected_delta_m, float(strict["M_strict_1L_GeV"]), tol, rel
        )
        assertions[f"{name}_strict_width"] = close(-d.imag / m0, float(strict["Gamma_strict_1L_GeV"]), tol, rel)

    neutral_keys = ("Delta_AA_at_z_GeV2", "Delta_AZ_at_z_GeV2", "Delta_ZA_at_z_GeV2")
    point_neutral = all(fixture_self_energy.get(name) is not None for name in neutral_keys)
    expected_claims = {
        name: bool(fixture.get("evidence_gates", {}).get(name, False))
        for name in REQUIRED_PROMOTION_GATES
        if name != "full_neutral_matrix_supplied"
    }
    expected_gates = {
        "strict_one_loop_pole_equation_implemented": True,
        "strict_coefficients_and_sqrt_readout_separated": True,
        "neutral_matrix_point_diagnostic_supplied": point_neutral,
        **{name: False for name in REQUIRED_PROMOTION_GATES},
    }
    assertions["unverified_claim_binding"] = (
        receipt["unverified_evidence_claims"] == expected_claims
    )
    assertions["evidence_gate_recomputation"] = receipt["evidence_gates"] == expected_gates
    assertions["promotion_required_gate_set"] = (
        receipt["promotion_required_gates"] == REQUIRED_PROMOTION_GATES
    )
    assertions["blocker_recomputation"] = receipt["blockers"] == REQUIRED_PROMOTION_GATES
    assertions["promotion_policy"] = receipt["promotion_policy"] == PROMOTION_POLICY
    assertions["promotion_boolean"] = receipt["physical_promotion_allowed"] is False
    assertions["status_recomputation"] = receipt["status"] == CONDITIONAL_STATUS
    assertions["order_guard"] = receipt["perturbative_order"] == "strict_one_loop"
    assertions["mixing_not_used_at_1L"] = (
        receipt["neutral_matrix"]["strict_1L_root_rule"]
        == STRICT_NEUTRAL_RULE
    )
    assertions["neutral_status"] = receipt["neutral_matrix"]["status"] == (
        "complete_at_declared_evaluation_point" if point_neutral else "missing_AA_AZ_ZA_coefficients"
    )
    assertions["neutral_presence"] = (
        receipt["neutral_matrix"]["diagnostics"] is not None
        if point_neutral
        else receipt["neutral_matrix"]["diagnostics"] is None
    )
    if point_neutral:
        dAA = c(fixture_self_energy["Delta_AA_at_z_GeV2"])
        dAZ = c(fixture_self_energy["Delta_AZ_at_z_GeV2"])
        dZA = c(fixture_self_energy["Delta_ZA_at_z_GeV2"])
        expected_matrix = (
            (complex(z, 0.0) - dAA, -dAZ),
            (-dZA, -dZ),
        )
        expected_det = (
            expected_matrix[0][0] * expected_matrix[1][1]
            - expected_matrix[0][1] * expected_matrix[1][0]
        )
        diagnostics = receipt["neutral_matrix"]["diagnostics"]
        got_matrix = diagnostics["Gamma_N_at_z_GeV2"]
        assertions["neutral_evaluation_point"] = close_complex(
            complex(z, 0.0), c(diagnostics["evaluation_s_GeV2"]), tol2, rel
        )
        for row in range(2):
            for column in range(2):
                assertions[f"neutral_matrix_{row}{column}"] = close_complex(
                    expected_matrix[row][column], c(got_matrix[row][column]), tol2, rel
                )
        assertions["neutral_determinant"] = close_complex(
            expected_det, c(diagnostics["det_Gamma_N_at_z_GeV4"]), tol2, rel
        )
        assertions["neutral_mixing_product"] = close_complex(
            dZA * dAZ,
            c(diagnostics["Delta_ZA_times_Delta_AZ_GeV4"]),
            tol2,
            rel,
        )
        assertions["neutral_mixing_power"] = diagnostics["mixing_loop_power"] == 2
        assertions["neutral_mixing_excluded"] = (
            diagnostics["mixing_used_in_strict_1L_root"] is False
        )
        assertions["neutral_effective_term"] = (
            diagnostics["two_loop_effective_inverse_term"] == NEUTRAL_EFFECTIVE_TERM
        )

    failed = sorted(name for name, ok in assertions.items() if not ok)
    if failed:
        raise SystemExit("FAIL: " + ", ".join(failed))
    return {
        "status": "PASS",
        "receipt": str(receipt_path),
        "assertions": assertions,
        "physical_promotion_allowed": False,
        "classification": receipt["status"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    package = Path(__file__).resolve().parents[1]
    parser.add_argument("--receipt", type=Path, default=package / "outputs" / "conditional_strict_1l_pole_map_receipt.json")
    parser.add_argument("--fixture", type=Path, default=package / "data" / "conditional_smdr_order1_fixture.json")
    parser.add_argument("--schema", type=Path, default=package / "schemas" / "physical_wz_strict_1l_pole_map_receipt.schema.json")
    args = parser.parse_args()
    result = check(args.receipt, args.fixture, args.schema)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
