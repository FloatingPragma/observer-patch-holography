#!/usr/bin/env python3
"""Certify the target-free leading photon-lepton threshold family for FZ-12.

The frozen FZ-12 packet supplies the photon coefficient.  This producer keeps
the electron and positron coefficients independent, states additive
energy-momentum conservation as a conditional preferred-frame premise, and
derives the leading head-on pair-production threshold.  It also pins the exact
``q <= 1`` spatial-symbol remainder without promoting it to a physical photon
mass shell.

No public comparison value, event row, source fit, detector response, or
future AugerPrime release is read.  The output is a kinematic and degeneracy
certificate, not an OPH score or a physical pair-production prediction.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import subprocess
from typing import Any

import sympy as sp

import a5_multipole_fixed_point_certificate as base


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
OBSERVATION_PATH = RUNTIME / "fz12_observation_map_receipt.json"
REMAINDER_PATH = RUNTIME / "fz12_full_symbol_remainder_receipt.json"
CUSTODY_PATH = RUNTIME / "fz12_custody_projection.json"
LEAN_PATH = REPO_ROOT / "Lean" / "Screen" / "SeamCurrentPhotonLeptonThreshold.lean"
LEAN_DIR = REPO_ROOT / "Lean"
LEAN_BUILD_TARGET = "SeamCurrentPhotonLeptonThreshold"
RECEIPT_PATH = RUNTIME / "fz12_joint_threshold_receipt.json"

LEAN_EXPORTED_THEOREMS = (
    "additive_head_on_implies_leading_threshold",
    "symmetric_share_coefficient",
    "charge_odd_lepton_direction_invisible",
    "common_shift_degeneracy",
    "charge_exchange_degeneracy",
    "fz12_universal_principal_symbol_coefficient",
    "li_lepton_mass_penalty_eq_symmetric_iff",
    "photon_only_residual_le_symmetric_share",
    "photon_only_residual_eq_symmetric_share_iff",
    "critical_threshold_polynomial_zero",
    "critical_threshold_derivative_zero",
    "exact_fz12_symbol_le_radiusSquared",
    "exact_fz12_frequency_le_euclideanMomentumMagnitude",
    "subluminal_photon_decay_impossible_with_li_leptons",
    "subluminal_pair_open_implies_li_pair_open",
)

SCHEMA = "oph.fz12.joint_photon_lepton_threshold.v1"
STATUS = (
    "EXACT_LEADING_JOINT_THRESHOLD__LI_LEPTON_BRANCH_UNIQUE_EQUAL_SHARE__"
    "GENERAL_EQUAL_SHARE_TWO_DIMENSIONAL_COEFFICIENT_FIBER__"
    "PHYSICAL_LEPTON_INTERACTION_SOURCE_AND_POWER_OPEN"
)

PARENT_CONTRACTS = (
    {
        "path": OBSERVATION_PATH,
        "relative_path": (
            "code/a5_fingerprint/runtime/fz12_observation_map_receipt.json"
        ),
        "bytes": 5276,
        "sha256": (
            "0b34218b075a3d51fe0badb00e3bb889743ae19b6dcf161f7877683e25121d17"
        ),
        "schema": "oph.fz12.formal_observation_map.v1",
        "status": (
            "EXACT_FORMAL_SERIES_MAP_FROM_FROZEN_FZ12_RAY__"
            "PHYSICAL_SECTOR_CLOCK_FRAME_REMAINDER_AND_COMPARISON_OPEN"
        ),
        "self_field": "receipt_sha256",
        "role": "target-free FZ-12 photon dispersion convention map",
    },
    {
        "path": REMAINDER_PATH,
        "relative_path": (
            "code/a5_fingerprint/runtime/fz12_full_symbol_remainder_receipt.json"
        ),
        "bytes": 5367,
        "sha256": (
            "ce4052f633c891515ebc319e9d7ff2bc0044bb8bb9e76a890e8e22c8e882dcfa"
        ),
        "schema": "oph.fz12.full_symbol_remainder.v1",
        "status": (
            "EXACT_TARGET_FREE_EDGE_SYMBOL_Q_LE_ONE_REMAINDER__"
            "PHYSICAL_FREQUENCY_AND_COMPARISON_OPEN"
        ),
        "self_field": "receipt_sha256",
        "role": "exact FZ-12 spatial symbol and q at most one remainder",
    },
    {
        "path": CUSTODY_PATH,
        "relative_path": "code/a5_fingerprint/runtime/fz12_custody_projection.json",
        "bytes": 3624,
        "sha256": (
            "dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643"
        ),
        "schema": "oph.fz12.custody_projection.v1",
        "status": "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA",
        "self_field": "projection_sha256",
        "role": "FZ-12-only custody projection without comparison values",
    },
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise base.FingerprintError(message)


def sha256_hex(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


@lru_cache(maxsize=1)
def verify_lean_kernel_build() -> str:
    """Build the pinned module and require axiom-audit output for every export."""
    completed = subprocess.run(
        ["lake", "build", LEAN_BUILD_TARGET],
        cwd=LEAN_DIR,
        check=False,
        capture_output=True,
        text=True,
    )
    output = completed.stdout + completed.stderr
    require(completed.returncode == 0, "threshold Lean kernel build failed")
    require("sorryAx" not in output, "threshold Lean kernel build contains sorryAx")
    namespace = "OPH.SeamCurrentPhotonLeptonThreshold."
    for theorem in LEAN_EXPORTED_THEOREMS:
        require(
            f"{namespace}{theorem}" in output,
            f"missing Lean axiom-audit output: {theorem}",
        )
    return output


def load_parent(contract: dict[str, Any]) -> tuple[bytes, dict[str, Any]]:
    raw = contract["path"].read_bytes()
    require(len(raw) == contract["bytes"], f"parent byte drift: {contract['path'].name}")
    require(
        sha256_hex(raw) == contract["sha256"],
        f"parent hash drift: {contract['path'].name}",
    )
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(
            f"invalid parent JSON: {contract['path'].name}"
        ) from error
    require(isinstance(value, dict), f"parent root drift: {contract['path'].name}")
    require(raw == base.canonical_json_bytes(value), f"noncanonical parent: {contract['path'].name}")
    require(value.get("schema") == contract["schema"], f"schema drift: {contract['path'].name}")
    require(value.get("status") == contract["status"], f"status drift: {contract['path'].name}")
    self_field = contract["self_field"]
    body = {key: item for key, item in value.items() if key != self_field}
    require(
        value.get(self_field) == base.tagged_sha256(base.canonical_json_bytes(body)),
        f"parent self-digest drift: {contract['path'].name}",
    )
    return raw, value


def load_lean_proof() -> bytes:
    verify_lean_kernel_build()
    raw = LEAN_PATH.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise base.FingerprintError("non-UTF-8 threshold Lean proof") from error
    normalized = re.sub(r"\s+", " ", text)
    required = (
        "theorem additive_head_on_implies_leading_threshold",
        "theorem symmetric_share_coefficient",
        "theorem charge_odd_lepton_direction_invisible",
        "theorem common_shift_degeneracy",
        "theorem symmetric_share_cancellation",
        "theorem fz12_universal_principal_symbol_coefficient",
        "theorem li_lepton_mass_penalty_eq_symmetric_iff",
        "theorem photon_only_residual_le_symmetric_share",
        "theorem photon_only_residual_eq_symmetric_share_iff",
        "theorem critical_threshold_polynomial_zero",
        "theorem exact_fz12_symbol_le_radiusSquared",
        "theorem exact_fz12_frequency_le_euclideanMomentumMagnitude",
        "theorem subluminal_photon_decay_impossible_with_li_leptons",
        "theorem subluminal_pair_open_implies_li_pair_open",
        "ordinary product norm on `Fin 3 -> R` is not used",
        "Nothing here selects the electron or positron action",
    )
    for fragment in required:
        require(fragment in normalized, "claim-bearing threshold proof drift")
    require("sorry" not in text and "admit" not in text, "placeholder in threshold proof")
    return raw


def symbolic_derivation() -> dict[str, Any]:
    delta_gamma, delta_plus, delta_minus = sp.symbols(
        "delta_gamma delta_plus delta_minus", real=True
    )
    energy, epsilon, mass, share = sp.symbols(
        "energy epsilon mass share", real=True, nonzero=True
    )
    photon_momentum = energy - delta_gamma * energy**3 / 2
    plus_momentum = (
        share * (energy + epsilon)
        - mass**2 / (2 * share * energy)
        - delta_plus * (share * energy) ** 3 / 2
    )
    minus_momentum = (
        (1 - share) * (energy + epsilon)
        - mass**2 / (2 * (1 - share) * energy)
        - delta_minus * ((1 - share) * energy) ** 3 / 2
    )
    effective = (
        delta_gamma
        - delta_plus * share**3
        - delta_minus * (1 - share) ** 3
    )
    residual = (
        effective * energy**4
        + 4 * epsilon * energy
        - mass**2 / (share * (1 - share))
    )
    balance = photon_momentum - epsilon - plus_momentum - minus_momentum
    require(
        sp.simplify(2 * energy * balance + residual) == 0,
        "additive threshold derivation drift",
    )

    half = sp.Rational(1, 2)
    symmetric = sp.simplify(effective.subs(share, half))
    expected = delta_gamma - (delta_plus + delta_minus) / 8
    require(sp.simplify(symmetric - expected) == 0, "symmetric combination drift")

    observable = sp.Matrix([[1, sp.Rational(-1, 8), sp.Rational(-1, 8)]])
    nullspace = observable.nullspace()
    require(observable.rank() == 1, "threshold observable rank drift")
    require(len(nullspace) == 2, "threshold coefficient-fiber dimension drift")
    charge_odd = sp.Matrix([0, 1, -1])
    common_shift = sp.Matrix([1, 4, 4])
    require(observable * charge_odd == sp.zeros(1, 1), "charge-odd fiber drift")
    require(observable * common_shift == sp.zeros(1, 1), "common-shift fiber drift")

    delta_fz12 = -sp.Symbol("a", real=True) ** 2 / 20
    photon_only = sp.simplify(expected.subs({delta_gamma: delta_fz12, delta_plus: 0, delta_minus: 0}))
    universal = sp.simplify(
        expected.subs(
            {
                delta_gamma: delta_fz12,
                delta_plus: delta_fz12,
                delta_minus: delta_fz12,
            }
        )
    )
    require(photon_only == delta_fz12, "FZ-12 photon-only combination drift")
    require(universal == -3 * sp.Symbol("a", real=True) ** 2 / 80, "universal combination drift")

    li_penalty_excess = mass**2 / (share * (1 - share)) - 4 * mass**2
    li_penalty_square = mass**2 * (2 * share - 1) ** 2 / (
        share * (1 - share)
    )
    require(
        sp.simplify(li_penalty_excess - li_penalty_square) == 0,
        "standard-lepton equal-share optimization identity drift",
    )

    critical_energy = 4 * mass**2 / (3 * epsilon)
    critical_delta = -27 * epsilon**4 / (64 * mass**6)
    polynomial = critical_delta * energy**4 + 4 * epsilon * energy - 4 * mass**2
    derivative = 4 * critical_delta * energy**3 + 4 * epsilon
    require(sp.simplify(polynomial.subs(energy, critical_energy)) == 0, "critical root drift")
    require(sp.simplify(derivative.subs(energy, critical_energy)) == 0, "critical derivative drift")

    wrong_soft_factor = effective * energy**4 + 2 * epsilon * energy - mass**2 / (
        share * (1 - share)
    )
    require(
        sp.simplify(2 * energy * balance + wrong_soft_factor) != 0,
        "wrong head-on soft factor escaped",
    )
    wrong_share_power = delta_gamma - delta_plus * share**2 - delta_minus * (1 - share) ** 2
    require(sp.simplify(wrong_share_power - effective) != 0, "wrong share power escaped")

    return {
        "preferred_frame": "declared frame in which the soft background is isotropic",
        "composition_law": "ordinary additive energy and three-momentum",
        "hard_photon_momentum_expansion": "p_gamma = E - delta_gamma,2 E^3/2",
        "outgoing_energy_partition": "E_plus = x(E+epsilon), E_minus = (1-x)(E+epsilon) at leading order",
        "charged_lepton_momentum_expansion": (
            "p_s = x_s(E+epsilon) - m_e^2/(2 x_s E) "
            "- delta_s,2 (x_s E)^3/2"
        ),
        "threshold_equation": (
            "[delta_gamma,2 - delta_plus,2 x^3 - "
            "delta_minus,2 (1-x)^3] E^4 + 4 epsilon E "
            "- m_e^2/[x(1-x)] = 0"
        ),
        "fixed_share_observable": (
            "delta_eff(x) = delta_gamma,2 - delta_plus,2 x^3 "
            "- delta_minus,2 (1-x)^3"
        ),
        "equal_share_observable": (
            "delta_eff(1/2) = delta_gamma,2 "
            "- (delta_plus,2 + delta_minus,2)/8"
        ),
        "equal_share_linear_map": ["1", "-1/8", "-1/8"],
        "equal_share_linear_rank": 1,
        "equal_share_coefficient_fiber_dimension": 2,
        "explicit_fiber_basis": [
            {
                "coordinates": ["0", "1", "-1"],
                "role": "charge-odd lepton direction",
            },
            {
                "coordinates": ["1", "4", "4"],
                "role": "common photon-lepton shift direction",
            },
        ],
        "charge_exchange_symmetry": (
            "(delta_plus,delta_minus,x) -> "
            "(delta_minus,delta_plus,1-x)"
        ),
        "leading_cancellation_plane": (
            "delta_plus,2 + delta_minus,2 = 8 delta_gamma,2"
        ),
        "standard_lepton_share_optimization": {
            "domain": "0 < x < 1 and m_e > 0",
            "mass_penalty_identity": (
                "m_e^2/[x(1-x)] - 4 m_e^2 = "
                "m_e^2(2x-1)^2/[x(1-x)]"
            ),
            "mass_penalty_result": (
                "the unique global minimum is 4 m_e^2 at x=1/2"
            ),
            "threshold_residual_result": (
                "with delta_plus,2=delta_minus,2=0, the unique global "
                "maximum of the leading head-on collinear residual is at x=1/2"
            ),
            "scope": (
                "Lorentz-invariant charged-lepton leading branch only; general "
                "independent-lepton and full anisotropic minimizations remain open"
            ),
        },
        "critical_double_root": {
            "hard_energy": "4 m_e^2/(3 epsilon)",
            "effective_coefficient": "-27 epsilon^4/(64 m_e^6)",
        },
    }


def verify_parent_semantics(parents: list[dict[str, Any]]) -> None:
    observation, remainder, custody = parents
    coefficients = observation["formal_input"]["frozen_coefficients"]
    require(
        coefficients["C4_over_a2"] == "-1/20",
        "frozen photon coefficient drift",
    )
    require(
        observation["exact_formal_result"]["frequency"].startswith(
            "omega = k - (a^2/40) k^3"
        ),
        "formal positive-frequency branch drift",
    )
    require(
        observation["exposure_boundary"]["comparison_permitted"] is False,
        "observation map exposure drift",
    )
    exact = remainder["exact_symbol"]
    require(
        exact["full_symbol"]
        == "lambda_hat(q,n) = (1/5) sum_e [1 - cos(q (w_e . n))]",
        "full symbol drift",
    )
    require(
        remainder["taylor_remainders"]["R8"]["global_upper"]
        == "(7/34992000) q^10",
        "q at most one remainder drift",
    )
    require(
        remainder["geometry_contract"]["q_domain"]
        == {"minimum": "0", "maximum": "1", "inclusive": True},
        "remainder domain drift",
    )
    scope = custody["projection_scope"]
    for key in (
        "includes_measurement_values",
        "includes_comparison_values",
        "includes_likelihood_values",
        "includes_other_campaign_rows",
    ):
        require(scope[key] is False, f"custody exposure drift: {key}")


def build_receipt() -> dict[str, Any]:
    loaded = [load_parent(contract) for contract in PARENT_CONTRACTS]
    parents = [value for _raw, value in loaded]
    verify_parent_semantics(parents)
    proof_raw = load_lean_proof()
    derivation = symbolic_derivation()

    parent_pins = []
    for contract, (raw, value) in zip(PARENT_CONTRACTS, loaded, strict=True):
        parent_pins.append(
            {
                "path": contract["relative_path"],
                "role": contract["role"],
                "bytes": len(raw),
                "sha256": f"sha256:{sha256_hex(raw)}",
                "schema": value["schema"],
                "status": value["status"],
                "self_digest": value[contract["self_field"]],
            }
        )
    parent_pins.append(
        {
            "path": LEAN_PATH.relative_to(REPO_ROOT).as_posix(),
            "role": (
                "kernel-checked leading threshold and equal-share "
                "coefficient-fiber algebra"
            ),
            "bytes": len(proof_raw),
            "sha256": f"sha256:{sha256_hex(proof_raw)}",
            "sorry_free": True,
            "theorems": list(LEAN_EXPORTED_THEOREMS),
        }
    )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 670,
        "parent_pins": parent_pins,
        "kinematic_scope": {
            "expansion_order": (
                "leading ultrarelativistic order in mass squared, soft energy, "
                "and dimension-six dispersion"
            ),
            "soft_photon_dispersion": "Lorentz-invariant leading order",
            "threshold_geometry": "head-on incoming and collinear outgoing momenta",
            "outgoing_share_domain": "0 < x < 1",
            "conservation_law_source_selected": False,
            "interaction_vertex_proved": False,
            "standard_lepton_collinear_share_optimization_proved": True,
            "general_independent_lepton_share_optimization_proved": False,
            "full_anisotropic_minimization_proved": False,
            "cross_section_proved": False,
            "exact_anisotropic_problem": (
                "minimize E_minus(p)+E_plus(K-p) over the full outgoing "
                "three-momentum; an interior minimizer has equal final group "
                "velocities, not necessarily collinear momenta"
            ),
            "leading_collinear_formula_promoted_to_exact_anisotropic_threshold": False,
        },
        "exact_leading_derivation": derivation,
        "fz12_conditional_branches": {
            "photon_coefficient": "delta_gamma,2 = -a^2/20",
            "standard_lepton_premise": {
                "delta_plus,2": "0",
                "delta_minus,2": "0",
                "visible_equal_share_coefficient": "-a^2/20",
                "source_selected": False,
            },
            "universal_principal_symbol_premise": {
                "delta_plus,2": "-a^2/20",
                "delta_minus,2": "-a^2/20",
                "visible_equal_share_coefficient": "-3a^2/80",
                "source_selected": False,
            },
            "unrestricted_lepton_result": (
                "the leading pair threshold cannot identify or bound the photon "
                "coefficient alone; the equal-share map has a two-dimensional fiber"
            ),
        },
        "exact_photon_symbol_control": {
            "dimensionless_symbol": (
                "lambda_hat(q,n) = (1/5) sum_e [1-cos(q w_e.n)]"
            ),
            "q_domain": "0 <= q <= 1",
            "global_exact_upper_bound": (
                "lambda_a(k) <= |k|_Euclidean^2 for every momentum and a!=0"
            ),
            "P8_remainder": (
                "0 <= lambda_hat-P8 <= (7/34992000) q^10"
            ),
            "mass_shell_remainder_if_frequency_premise_is_added": (
                "0 <= E_gamma^2-a^-2 P8 <= "
                "(7/34992000) a^-2 q^10"
            ),
            "physical_frequency_premise_proved": False,
            "charged_lepton_exact_symbol_available": False,
            "statement": (
                "the exact photon remainder is pinned; a physical threshold "
                "interval additionally needs the photon frequency attachment and "
                "independently controlled charged-lepton symbols"
            ),
        },
        "exact_li_lepton_consequences": {
            "premises": [
                "the spatial symbol is identified with physical photon frequency squared",
                "the exact cosine-symbol inequality lambda_a(k) <= |k|_Euclidean^2 holds for a!=0",
                "electron and positron use standard positive-energy relativistic dispersion",
                "energy and momentum compose additively in the declared frame",
            ],
            "subluminal_frequency": "Omega_gamma(k,n) <= |k|",
            "photon_decay": (
                "kinematically impossible because every LI lepton pair with "
                "total momentum k has energy strictly greater than |k|"
            ),
            "pair_production_domain": (
                "at fixed incoming momenta, Lorentz-invariant incoming energy is "
                "no smaller than FZ incoming energy, so the FZ energy-budget "
                "domain is contained in the Lorentz-invariant energy-budget "
                "domain; no final-state existence, reachability, vertex, cross "
                "section, opacity, or rate is inferred"
            ),
            "EFT_truncation_used": False,
            "rank_six_direction_restriction_used": False,
            "physical_photon_attachment_proved": False,
            "standard_lepton_premise_source_selected": False,
        },
        "negative_controls": {
            "all_fired": True,
            "controls": [
                "replace the head-on 4 epsilon E term by 2 epsilon E",
                "replace the outgoing share cubes by squares",
                "move the standard-lepton residual maximizer away from equal sharing",
                "erase either independent charged-lepton coefficient",
                "promote the additive composition law to an A1-A3 theorem",
                "promote the universal-principal-symbol branch to source-selected",
                "use the exact photon remainder above q=1",
                "interpret the rank-one observable map as three identified coefficients",
            ],
        },
        "exposure_boundary": {
            "comparison_inputs": [],
            "public_measurement_read": False,
            "future_event_rows_read": False,
            "existing_Auger_limit_read": False,
            "comparison_permitted": False,
            "score_emitted": False,
            "evidence_claimed": False,
            "verdict_emitted": False,
        },
        "open_physical_attachments": {
            "FZ12_oscillator_identified_with_photon": False,
            "electron_action_selected": False,
            "positron_action_selected": False,
            "energy_momentum_composition_selected": False,
            "pair_vertex_and_cross_section_selected": False,
            "background_transport_selected": False,
            "source_spectrum_and_composition_selected": False,
            "atmospheric_shower_response_selected": False,
            "detector_classification_selected": False,
            "target_free_power_envelope_complete": False,
            "prospective_contract_armed": False,
        },
    }
    receipt["receipt_sha256"] = base.tagged_sha256(base.canonical_json_bytes(receipt))
    return receipt


def verify_committed_receipt() -> dict[str, Any]:
    raw = RECEIPT_PATH.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid joint-threshold receipt JSON") from error
    require(isinstance(committed, dict), "joint-threshold receipt root drift")
    require(raw == base.canonical_json_bytes(committed), "noncanonical joint-threshold receipt")
    require(committed.get("schema") == SCHEMA, "joint-threshold schema drift")
    require(committed.get("status") == STATUS, "joint-threshold status drift")
    body = {key: item for key, item in committed.items() if key != "receipt_sha256"}
    require(
        committed.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        "joint-threshold self-digest drift",
    )
    rebuilt = build_receipt()
    require(raw == base.canonical_json_bytes(rebuilt), "joint-threshold ancestry or result drift")
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if args.write and args.verify:
        parser.error("choose --write or --verify")
    receipt = verify_committed_receipt() if args.verify else build_receipt()
    if args.write:
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
