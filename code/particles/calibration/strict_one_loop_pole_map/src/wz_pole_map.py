#!/usr/bin/env python3
"""Strict-one-loop W/Z complex-pole map.

This module deliberately separates three objects:

1. strict perturbative pole coefficients, used for finite-order gauge checks;
2. the truncated complex pole s_V through one loop; and
3. the conventional energy-pole coordinates M_V and Gamma_V obtained by
   taking the square root of the truncated s_V.

The inverse-propagator convention is

    Gamma_W^T(s) = s - w - Delta_WW^(1)(s) + O(h^2),

    Gamma_N^T(s) = [[s - Delta_AA^(1)(s), -Delta_AZ^(1)(s)],
                    [-Delta_ZA^(1)(s), s-z-Delta_ZZ^(1)(s)]] + O(h^2).

Each Delta^(1) passed here already includes the loop factor, counterterms and
all terms retained by the declared one-loop mask.  Thus

    s_W^[1] = w + Delta_WW^(1)(w),
    s_Z^[1] = z + Delta_ZZ^(1)(z).

The product Delta_ZA^(1) Delta_AZ^(1) is loop order two.  It may be retained
in a neutral-matrix diagnostic but is never inserted into a strict-one-loop
Z root by this module.
"""

from __future__ import annotations

import argparse
import cmath
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA_ID = "physical_wz_strict_1l_pole_map_receipt_v1"
CONVENTION_ID = "Gamma=s-m2-Delta; s_pole=m2+Delta_at_tree"
SELF_ENERGY_TRANSLATION = "Delta=-Pi; Gamma=s-m2+Pi=s-m2-Delta"

EXTERNAL_EVIDENCE_GATES = (
    "full_neutral_matrix_supplied",
    "renormalized_FJ_vev_selected",
    "complete_FJ_tadpole_conversion_checked",
    "target_clean_source_matching",
    "source_law_or_deterministic_exactness_receipt",
    "independent_self_energy_engine",
    "finite_order_BRST_ST_Ward_Nielsen_receipt",
    "physical_clock_closed",
    "no_target_ancestry",
)


class PoleMapError(ValueError):
    """Raised when a receipt would violate the strict-one-loop contract."""


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_complex(value: Sequence[float], *, field: str) -> complex:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 2:
        raise PoleMapError(f"{field} must be [real, imaginary]")
    out = complex(float(value[0]), float(value[1]))
    if not (math.isfinite(out.real) and math.isfinite(out.imag)):
        raise PoleMapError(f"{field} must be finite")
    return out


def complex_pair(value: complex) -> list[float]:
    return [float(value.real), float(value.imag)]


def checked_positive(name: str, value: float) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise PoleMapError(f"{name} must be positive and finite")
    return out


@dataclass(frozen=True)
class TreeMasses:
    w: float
    z: float

    @property
    def mW0(self) -> float:
        return math.sqrt(self.w)

    @property
    def mZ0(self) -> float:
        return math.sqrt(self.z)


def tree_masses(g: float, gp: float, v_f: float) -> TreeMasses:
    """Return w=g^2 v_F^2/4 and z=(g^2+g'^2)v_F^2/4."""

    g = checked_positive("g", g)
    gp = checked_positive("gp", gp)
    v_f = checked_positive("v_f", v_f)
    w = g * g * v_f * v_f / 4.0
    z = (g * g + gp * gp) * v_f * v_f / 4.0
    if not z > w > 0.0:
        raise PoleMapError("tree masses must obey z > w > 0")
    return TreeMasses(w=w, z=z)


def energy_pole_coordinates(s_pole: complex) -> dict[str, Any]:
    """Return the branch M>0, Gamma>=0 for s=(M-i Gamma/2)^2."""

    if not (math.isfinite(s_pole.real) and math.isfinite(s_pole.imag)):
        raise PoleMapError("complex pole must be finite")
    root = cmath.sqrt(s_pole)
    if root.real < 0.0 or (root.real == 0.0 and root.imag > 0.0):
        root = -root
    if root.imag > 0.0:
        # A physical decaying-pole input should be below the real axis.  Do not
        # silently conjugate it, because that would hide a sign-convention bug.
        raise PoleMapError("pole is on the growing-state sheet (Im sqrt(s)>0)")
    mass = float(root.real)
    width = float(-2.0 * root.imag)
    if mass <= 0.0 or width < 0.0:
        raise PoleMapError("invalid energy-pole branch")
    reconstructed = complex(mass, -0.5 * width) ** 2
    return {
        "M_GeV": mass,
        "Gamma_GeV": width,
        "sqrt_s_GeV": complex_pair(root),
        "s_reconstructed_GeV2": complex_pair(reconstructed),
        "coordinate_transform_note": (
            "Exact square root of the one-loop-truncated s pole; the nonlinear "
            "coordinate readout contains kinematic terms beyond strict one loop."
        ),
    }


def strict_mass_width_coefficients(m0: float, delta: complex) -> dict[str, float]:
    """Linearized one-loop mass/width coefficients from s=m0^2+Delta."""

    m0 = checked_positive("tree mass", m0)
    delta_m = float(delta.real / (2.0 * m0))
    gamma = float(-delta.imag / m0)
    if gamma < -1e-14:
        raise PoleMapError("one-loop width coefficient is negative; check the sheet/sign convention")
    gamma = max(0.0, gamma)
    return {
        "M0_GeV": m0,
        "delta_M_1L_GeV": delta_m,
        "M_strict_1L_GeV": m0 + delta_m,
        "Gamma_strict_1L_GeV": gamma,
    }


def neutral_inverse_matrix(
    s: complex,
    *,
    z: float,
    delta_aa: complex,
    delta_az: complex,
    delta_za: complex,
    delta_zz: complex,
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    """Neutral transverse inverse matrix at a declared evaluation point."""

    return (
        (s - delta_aa, -delta_az),
        (-delta_za, s - z - delta_zz),
    )


def det2(matrix: tuple[tuple[complex, complex], tuple[complex, complex]]) -> complex:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def neutral_diagnostics_at_tree(
    *,
    z: float,
    delta_aa: complex,
    delta_az: complex,
    delta_za: complex,
    delta_zz: complex,
) -> dict[str, Any]:
    """Emit the full neutral determinant at s=z without using mixing in the 1L root."""

    matrix = neutral_inverse_matrix(
        complex(z, 0.0),
        z=z,
        delta_aa=delta_aa,
        delta_az=delta_az,
        delta_za=delta_za,
        delta_zz=delta_zz,
    )
    determinant = det2(matrix)
    mixing_product = delta_za * delta_az
    return {
        "evaluation_s_GeV2": [z, 0.0],
        "Gamma_N_at_z_GeV2": [
            [complex_pair(matrix[0][0]), complex_pair(matrix[0][1])],
            [complex_pair(matrix[1][0]), complex_pair(matrix[1][1])],
        ],
        "det_Gamma_N_at_z_GeV4": complex_pair(determinant),
        "Delta_ZA_times_Delta_AZ_GeV4": complex_pair(mixing_product),
        "mixing_loop_power": 2,
        "mixing_used_in_strict_1L_root": False,
        "two_loop_effective_inverse_term": "-Delta_ZA^(1)(s)*Delta_AZ^(1)(s)/s",
    }


def validate_order_contract(fixture: Mapping[str, Any]) -> None:
    order = fixture.get("perturbative_order")
    if order != "strict_one_loop":
        raise PoleMapError(f"strict-one-loop map refuses perturbative_order={order!r}")
    if fixture.get("loop_factor_included_in_delta") is not True:
        raise PoleMapError("Delta coefficients must already include the loop factor")
    mask = fixture.get("contribution_mask")
    if not isinstance(mask, Mapping):
        raise PoleMapError("missing contribution_mask")
    if mask.get("maximum_loop_power") != 1:
        raise PoleMapError("maximum_loop_power must be 1")
    if mask.get("one_loop_complete_for_declared_backend") is not True:
        raise PoleMapError("declared backend one-loop mask is incomplete")


def _unverified_evidence_claims(fixture: Mapping[str, Any]) -> dict[str, bool]:
    evidence = fixture.get("evidence_gates", {})
    return {
        name: bool(evidence.get(name, False))
        for name in EXTERNAL_EVIDENCE_GATES
        if name != "full_neutral_matrix_supplied"
    }


def _gate_summary(neutral_point_supplied: bool) -> dict[str, bool]:
    """Return only gates proved by this algebraic v1 producer.

    External scientific evidence can never self-attest through fixture
    booleans.  A later aggregate checker must resolve and verify hash-bound
    evidence artifacts against this receipt's exact subject digest.
    """

    gates = {
        "strict_one_loop_pole_equation_implemented": True,
        "strict_coefficients_and_sqrt_readout_separated": True,
        "neutral_matrix_point_diagnostic_supplied": neutral_point_supplied,
    }
    gates.update({name: False for name in EXTERNAL_EVIDENCE_GATES})
    return gates


def build_receipt(fixture: Mapping[str, Any], *, fixture_path: Path | None = None) -> dict[str, Any]:
    """Build a fail-closed strict-one-loop pole-map receipt."""

    validate_order_contract(fixture)
    if fixture.get("inverse_propagator_convention") != CONVENTION_ID:
        raise PoleMapError("inverse-propagator convention mismatch")

    params = fixture["renormalized_parameters"]
    masses = tree_masses(float(params["g"]), float(params["gp"]), float(params["v_F_GeV"]))
    self_energy = fixture["renormalized_one_loop_self_energy"]
    delta_w = parse_complex(self_energy["Delta_WW_at_w_GeV2"], field="Delta_WW_at_w_GeV2")
    delta_z = parse_complex(self_energy["Delta_ZZ_at_z_GeV2"], field="Delta_ZZ_at_z_GeV2")

    s_w = complex(masses.w, 0.0) + delta_w
    s_z = complex(masses.z, 0.0) + delta_z
    w_readout = energy_pole_coordinates(s_w)
    z_readout = energy_pole_coordinates(s_z)
    w_strict = strict_mass_width_coefficients(masses.mW0, delta_w)
    z_strict = strict_mass_width_coefficients(masses.mZ0, delta_z)

    neutral_keys = ("Delta_AA_at_z_GeV2", "Delta_AZ_at_z_GeV2", "Delta_ZA_at_z_GeV2")
    neutral_point_supplied = all(self_energy.get(key) is not None for key in neutral_keys)
    neutral = None
    if neutral_point_supplied:
        neutral = neutral_diagnostics_at_tree(
            z=masses.z,
            delta_aa=parse_complex(self_energy["Delta_AA_at_z_GeV2"], field="Delta_AA_at_z_GeV2"),
            delta_az=parse_complex(self_energy["Delta_AZ_at_z_GeV2"], field="Delta_AZ_at_z_GeV2"),
            delta_za=parse_complex(self_energy["Delta_ZA_at_z_GeV2"], field="Delta_ZA_at_z_GeV2"),
            delta_zz=delta_z,
        )

    gates = _gate_summary(neutral_point_supplied)
    unverified_claims = _unverified_evidence_claims(fixture)
    promotion_required = list(EXTERNAL_EVIDENCE_GATES)
    # This v1 producer proves elementary pole algebra only.  It deliberately
    # cannot promote itself from caller-provided booleans.  Promotion belongs
    # to a separate aggregate verifier that resolves hash-bound evidence and
    # binds every receipt to this exact numerical subject.
    promotion_allowed = False

    source = dict(fixture.get("source", {}))
    if fixture_path is not None:
        source["fixture_path"] = fixture_path.name
        source["fixture_sha256"] = sha256_file(fixture_path)

    return {
        "schema": SCHEMA_ID,
        "receipt_id": fixture["receipt_id"],
        "claim_class": fixture["claim_class"],
        "status": "CONDITIONAL_STRICT_1L_POLE_MAP_NOT_OPH_NATIVE_PHYSICAL",
        "map_evaluated": True,
        "physical_promotion_allowed": promotion_allowed,
        "inverse_propagator_convention": CONVENTION_ID,
        "self_energy_sign_translation": SELF_ENERGY_TRANSLATION,
        "pole_convention": "s_V=(M_V-i*Gamma_V/2)^2; Re(sqrt(s_V))>0; Im(sqrt(s_V))<=0",
        "perturbative_order": "strict_one_loop",
        "source": source,
        "renormalized_parameters": dict(params),
        "contribution_mask": dict(fixture["contribution_mask"]),
        "tree": {
            "w_GeV2": masses.w,
            "z_GeV2": masses.z,
            "mW0_GeV": masses.mW0,
            "mZ0_GeV": masses.mZ0,
        },
        "one_loop_coefficients": {
            "Delta_WW_at_w_GeV2": complex_pair(delta_w),
            "Delta_ZZ_at_z_GeV2": complex_pair(delta_z),
            "strict_W": w_strict,
            "strict_Z": z_strict,
        },
        "complex_poles": {
            "W": {
                "s_pole_truncated_1L_GeV2": complex_pair(s_w),
                "energy_pole_readout": w_readout,
            },
            "Z": {
                "s_pole_truncated_1L_GeV2": complex_pair(s_z),
                "energy_pole_readout": z_readout,
            },
        },
        "neutral_matrix": {
            "status": "complete_at_declared_evaluation_point" if neutral_point_supplied else "missing_AA_AZ_ZA_coefficients",
            "diagnostics": neutral,
            "strict_1L_root_rule": "s_Z=z+Delta_ZZ^(1)(z); off-diagonal product excluded as O(h^2)",
        },
        "source_covariance": dict(fixture.get("source_covariance", {"status": "not_supplied"})),
        "evidence_gates": gates,
        "unverified_evidence_claims": unverified_claims,
        "promotion_required_gates": promotion_required,
        "blockers": [name for name in promotion_required if not gates[name]],
        "promotion_policy": (
            "v1 is an algebraic diagnostic and cannot self-promote; a separate "
            "aggregate checker must verify hash-bound evidence artifacts and subject binding"
        ),
        "tolerances": dict(fixture["tolerances"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fixture = json.loads(args.fixture.read_text(encoding="utf-8"))
    receipt = build_receipt(fixture, fixture_path=args.fixture)
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
