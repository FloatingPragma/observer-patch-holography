#!/usr/bin/env python3
"""Independent verifier for the FZ-12 joint photon-lepton threshold packet."""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import subprocess
from typing import Any

import sympy as sp


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
DEFAULT_RECEIPT = HERE / "runtime" / "fz12_joint_threshold_receipt.json"
LEAN_DIR = REPO_ROOT / "Lean"
LEAN_BUILD_TARGET = "SeamCurrentPhotonLeptonThreshold"

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

EXPECTED_TOP_LEVEL_KEYS = {
    "exact_leading_derivation",
    "exact_li_lepton_consequences",
    "exact_photon_symbol_control",
    "exposure_boundary",
    "fz12_conditional_branches",
    "issue",
    "kinematic_scope",
    "negative_controls",
    "open_physical_attachments",
    "parent_pins",
    "receipt_sha256",
    "schema",
    "status",
}

EXPECTED_PARENT_PINS = [
    {
        "bytes": 5276,
        "path": "code/a5_fingerprint/runtime/fz12_observation_map_receipt.json",
        "role": "target-free FZ-12 photon dispersion convention map",
        "schema": "oph.fz12.formal_observation_map.v1",
        "self_digest": "sha256:c426a35029de3f2c25347313153c9cf704827869b184c1a922fb0da7804bc3ee",
        "sha256": "sha256:0b34218b075a3d51fe0badb00e3bb889743ae19b6dcf161f7877683e25121d17",
        "status": (
            "EXACT_FORMAL_SERIES_MAP_FROM_FROZEN_FZ12_RAY__"
            "PHYSICAL_SECTOR_CLOCK_FRAME_REMAINDER_AND_COMPARISON_OPEN"
        ),
    },
    {
        "bytes": 5367,
        "path": "code/a5_fingerprint/runtime/fz12_full_symbol_remainder_receipt.json",
        "role": "exact FZ-12 spatial symbol and q at most one remainder",
        "schema": "oph.fz12.full_symbol_remainder.v1",
        "self_digest": "sha256:78ac6e509a3eb4bcb6d5ae15d191f8ebc473dad1e5554f139d9f3dbab6b5525a",
        "sha256": "sha256:ce4052f633c891515ebc319e9d7ff2bc0044bb8bb9e76a890e8e22c8e882dcfa",
        "status": (
            "EXACT_TARGET_FREE_EDGE_SYMBOL_Q_LE_ONE_REMAINDER__"
            "PHYSICAL_FREQUENCY_AND_COMPARISON_OPEN"
        ),
    },
    {
        "bytes": 3624,
        "path": "code/a5_fingerprint/runtime/fz12_custody_projection.json",
        "role": "FZ-12-only custody projection without comparison values",
        "schema": "oph.fz12.custody_projection.v1",
        "self_digest": "sha256:334d9877023d4eaa2ac64bca51d9e2193f0eaad81afd664b22cada2b890a8527",
        "sha256": "sha256:dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643",
        "status": "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA",
    },
    {
        "bytes": 23814,
        "path": "Lean/Screen/SeamCurrentPhotonLeptonThreshold.lean",
        "role": (
            "kernel-checked leading threshold and equal-share "
            "coefficient-fiber algebra"
        ),
        "sha256": "sha256:ae5b378d1449aee3bfa3cc6c5bd2f070682d690a3619c190f60e54c4e0cd670c",
        "sorry_free": True,
        "theorems": [
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
        ],
    },
]


class VerificationError(ValueError):
    """The independent joint-threshold verification failed closed."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


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
    check(completed.returncode == 0, "threshold Lean kernel build failed")
    check("sorryAx" not in output, "threshold Lean kernel build contains sorryAx")
    namespace = "OPH.SeamCurrentPhotonLeptonThreshold."
    for theorem in LEAN_EXPORTED_THEOREMS:
        check(
            f"{namespace}{theorem}" in output,
            f"missing Lean axiom-audit output: {theorem}",
        )
    return output


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("ascii")


def tagged_digest(raw: bytes) -> str:
    return f"sha256:{hashlib.sha256(raw).hexdigest()}"


def verify_self_digest(value: dict[str, Any], field: str) -> None:
    body = {key: item for key, item in value.items() if key != field}
    check(value.get(field) == tagged_digest(canonical_json_bytes(body)), f"{field} drift")


def load_json(path: Path) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()

    def object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for key, value in pairs:
            check(key not in out, f"duplicate JSON key {key!r}: {path}")
            out[key] = value
        return out

    def reject_constant(value: str) -> Any:
        raise VerificationError(f"non-finite JSON constant {value!r}: {path}")

    try:
        value = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=object_pairs,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise VerificationError(f"invalid JSON: {path}") from error
    check(isinstance(value, dict), f"JSON root is not an object: {path}")
    check(raw == canonical_json_bytes(value), f"noncanonical JSON: {path}")
    return raw, value


def verify_parent(
    pin: dict[str, Any], expected: dict[str, Any]
) -> dict[str, Any] | None:
    check(pin == expected, f"parent pin contract drift: {expected['path']}")
    path = REPO_ROOT / pin["path"]
    raw = path.read_bytes()
    check(len(raw) == pin["bytes"], f"parent byte drift: {path.name}")
    check(tagged_digest(raw) == pin["sha256"], f"parent digest drift: {path.name}")
    if path.suffix == ".lean":
        text = raw.decode("utf-8")
        check("sorry" not in text and "admit" not in text, "Lean placeholder")
        for theorem in pin["theorems"]:
            check(
                re.search(rf"\btheorem\s+{re.escape(theorem)}\b", text) is not None,
                f"missing Lean theorem: {theorem}",
            )
        return None
    parent_raw, parent = load_json(path)
    check(parent_raw == raw, "parent reread drift")
    check(parent.get("schema") == pin["schema"], f"parent schema drift: {path.name}")
    check(parent.get("status") == pin["status"], f"parent status drift: {path.name}")
    self_field = "projection_sha256" if "projection_sha256" in parent else "receipt_sha256"
    verify_self_digest(parent, self_field)
    check(parent[self_field] == pin["self_digest"], f"parent self pin drift: {path.name}")
    return parent


def recompute_algebra() -> dict[str, Any]:
    dg, dp, dm, energy, epsilon, mass, share = sp.symbols(
        "dg dp dm energy epsilon mass share", real=True, nonzero=True
    )
    p_gamma = energy - dg * energy**3 / 2
    p_plus = (
        share * (energy + epsilon)
        - mass**2 / (2 * share * energy)
        - dp * (share * energy) ** 3 / 2
    )
    p_minus = (
        (1 - share) * (energy + epsilon)
        - mass**2 / (2 * (1 - share) * energy)
        - dm * ((1 - share) * energy) ** 3 / 2
    )
    effective = dg - dp * share**3 - dm * (1 - share) ** 3
    residual = effective * energy**4 + 4 * epsilon * energy - mass**2 / (
        share * (1 - share)
    )
    balance = p_gamma - epsilon - p_plus - p_minus
    check(sp.simplify(2 * energy * balance + residual) == 0, "threshold identity failed")

    half_effective = sp.simplify(effective.subs(share, sp.Rational(1, 2)))
    check(
        sp.simplify(half_effective - (dg - (dp + dm) / 8)) == 0,
        "equal-share combination failed",
    )
    row = sp.Matrix([[1, sp.Rational(-1, 8), sp.Rational(-1, 8)]])
    check(row.rank() == 1 and len(row.nullspace()) == 2, "coefficient fiber failed")
    check(row * sp.Matrix([0, 1, -1]) == sp.zeros(1, 1), "charge-odd null failed")
    check(row * sp.Matrix([1, 4, 4]) == sp.zeros(1, 1), "common-shift null failed")

    a = sp.symbols("a", real=True)
    delta_fz12 = -a**2 / 20
    universal = (dg - (dp + dm) / 8).subs(
        {dg: delta_fz12, dp: delta_fz12, dm: delta_fz12}
    )
    check(sp.simplify(universal + 3 * a**2 / 80) == 0, "universal branch failed")

    li_penalty_excess = mass**2 / (share * (1 - share)) - 4 * mass**2
    li_penalty_square = mass**2 * (2 * share - 1) ** 2 / (
        share * (1 - share)
    )
    check(
        sp.simplify(li_penalty_excess - li_penalty_square) == 0,
        "standard-lepton equal-share optimization identity failed",
    )

    ecrit = 4 * mass**2 / (3 * epsilon)
    dcrit = -27 * epsilon**4 / (64 * mass**6)
    polynomial = dcrit * energy**4 + 4 * epsilon * energy - 4 * mass**2
    derivative = 4 * dcrit * energy**3 + 4 * epsilon
    check(sp.simplify(polynomial.subs(energy, ecrit)) == 0, "critical root failed")
    check(sp.simplify(derivative.subs(energy, ecrit)) == 0, "double root failed")

    return {
        "threshold_equation": (
            "[delta_gamma,2 - delta_plus,2 x^3 - "
            "delta_minus,2 (1-x)^3] E^4 + 4 epsilon E "
            "- m_e^2/[x(1-x)] = 0"
        ),
        "equal_share_observable": (
            "delta_eff(1/2) = delta_gamma,2 "
            "- (delta_plus,2 + delta_minus,2)/8"
        ),
        "rank": 1,
        "fiber_dimension": 2,
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
        "critical_energy": "4 m_e^2/(3 epsilon)",
        "critical_delta": "-27 epsilon^4/(64 m_e^6)",
    }


def verify(receipt_path: Path) -> dict[str, Any]:
    verify_lean_kernel_build()
    raw, receipt = load_json(receipt_path)
    check(set(receipt) == EXPECTED_TOP_LEVEL_KEYS, "receipt top-level key drift")
    check(receipt.get("schema") == SCHEMA, "receipt schema drift")
    check(receipt.get("status") == STATUS, "receipt status drift")
    check(receipt.get("issue") == 670, "receipt issue drift")
    verify_self_digest(receipt, "receipt_sha256")

    pins = receipt.get("parent_pins")
    check(isinstance(pins, list) and len(pins) == 4, "parent pin count drift")
    check(pins == EXPECTED_PARENT_PINS, "parent pin inventory drift")
    parents = [
        verify_parent(pin, expected)
        for pin, expected in zip(pins, EXPECTED_PARENT_PINS, strict=True)
    ]
    observation, remainder, custody, lean = parents
    check(lean is None, "Lean pin typed as JSON")
    assert observation is not None and remainder is not None and custody is not None
    check(
        observation["formal_input"]["frozen_coefficients"]["C4_over_a2"] == "-1/20",
        "FZ-12 photon coefficient drift",
    )
    check(
        remainder["taylor_remainders"]["R8"]["global_upper"]
        == "(7/34992000) q^10",
        "exact remainder drift",
    )
    for key in (
        "includes_measurement_values",
        "includes_comparison_values",
        "includes_likelihood_values",
        "includes_other_campaign_rows",
    ):
        check(custody["projection_scope"][key] is False, f"custody exposure: {key}")

    algebra = recompute_algebra()
    exact = receipt["exact_leading_derivation"]
    check(
        set(exact)
        == {
            "charge_exchange_symmetry",
            "charged_lepton_momentum_expansion",
            "composition_law",
            "critical_double_root",
            "equal_share_coefficient_fiber_dimension",
            "equal_share_linear_map",
            "equal_share_linear_rank",
            "equal_share_observable",
            "explicit_fiber_basis",
            "fixed_share_observable",
            "hard_photon_momentum_expansion",
            "leading_cancellation_plane",
            "outgoing_energy_partition",
            "preferred_frame",
            "standard_lepton_share_optimization",
            "threshold_equation",
        },
        "leading derivation key drift",
    )
    check(exact["threshold_equation"] == algebra["threshold_equation"], "threshold text drift")
    check(exact["equal_share_observable"] == algebra["equal_share_observable"], "observable text drift")
    check(exact["equal_share_linear_rank"] == algebra["rank"], "rank drift")
    check(
        exact["equal_share_coefficient_fiber_dimension"] == algebra["fiber_dimension"],
        "fiber dimension drift",
    )
    check(
        exact["standard_lepton_share_optimization"]
        == algebra["standard_lepton_share_optimization"],
        "standard-lepton share optimization drift",
    )
    check(
        exact["critical_double_root"]
        == {
            "hard_energy": algebra["critical_energy"],
            "effective_coefficient": algebra["critical_delta"],
        },
        "critical boundary drift",
    )
    check(
        exact["composition_law"] == "ordinary additive energy and three-momentum"
        and exact["preferred_frame"]
        == "declared frame in which the soft background is isotropic"
        and exact["equal_share_linear_map"] == ["1", "-1/8", "-1/8"]
        and exact["explicit_fiber_basis"]
        == [
            {
                "coordinates": ["0", "1", "-1"],
                "role": "charge-odd lepton direction",
            },
            {
                "coordinates": ["1", "4", "4"],
                "role": "common photon-lepton shift direction",
            },
        ],
        "leading derivation contract drift",
    )

    check(
        receipt["kinematic_scope"]
        == {
            "conservation_law_source_selected": False,
            "cross_section_proved": False,
            "exact_anisotropic_problem": (
                "minimize E_minus(p)+E_plus(K-p) over the full outgoing "
                "three-momentum; an interior minimizer has equal final group "
                "velocities, not necessarily collinear momenta"
            ),
            "expansion_order": (
                "leading ultrarelativistic order in mass squared, soft energy, "
                "and dimension-six dispersion"
            ),
            "interaction_vertex_proved": False,
            "leading_collinear_formula_promoted_to_exact_anisotropic_threshold": False,
            "outgoing_share_domain": "0 < x < 1",
            "soft_photon_dispersion": "Lorentz-invariant leading order",
            "standard_lepton_collinear_share_optimization_proved": True,
            "general_independent_lepton_share_optimization_proved": False,
            "full_anisotropic_minimization_proved": False,
            "threshold_geometry": "head-on incoming and collinear outgoing momenta",
        },
        "kinematic-scope promotion drift",
    )

    branches = receipt["fz12_conditional_branches"]
    check(
        branches
        == {
            "photon_coefficient": "delta_gamma,2 = -a^2/20",
            "standard_lepton_premise": {
                "delta_minus,2": "0",
                "delta_plus,2": "0",
                "source_selected": False,
                "visible_equal_share_coefficient": "-a^2/20",
            },
            "universal_principal_symbol_premise": {
                "delta_minus,2": "-a^2/20",
                "delta_plus,2": "-a^2/20",
                "source_selected": False,
                "visible_equal_share_coefficient": "-3a^2/80",
            },
            "unrestricted_lepton_result": (
                "the leading pair threshold cannot identify or bound the photon "
                "coefficient alone; the equal-share map has a two-dimensional fiber"
            ),
        },
        "conditional-branch contract drift",
    )
    check(
        receipt["exact_photon_symbol_control"]
        == {
            "P8_remainder": "0 <= lambda_hat-P8 <= (7/34992000) q^10",
            "charged_lepton_exact_symbol_available": False,
            "dimensionless_symbol": (
                "lambda_hat(q,n) = (1/5) sum_e [1-cos(q w_e.n)]"
            ),
            "global_exact_upper_bound": (
                "lambda_a(k) <= |k|_Euclidean^2 for every momentum and a!=0"
            ),
            "mass_shell_remainder_if_frequency_premise_is_added": (
                "0 <= E_gamma^2-a^-2 P8 <= (7/34992000) a^-2 q^10"
            ),
            "physical_frequency_premise_proved": False,
            "q_domain": "0 <= q <= 1",
            "statement": (
                "the exact photon remainder is pinned; a physical threshold "
                "interval additionally needs the photon frequency attachment and "
                "independently controlled charged-lepton symbols"
            ),
        },
        "photon-symbol control drift",
    )
    exact_li = receipt["exact_li_lepton_consequences"]
    check(
        exact_li
        == {
            "EFT_truncation_used": False,
            "pair_production_domain": (
                "at fixed incoming momenta, Lorentz-invariant incoming energy is "
                "no smaller than FZ incoming energy, so the FZ energy-budget "
                "domain is contained in the Lorentz-invariant energy-budget "
                "domain; no final-state existence, reachability, vertex, cross "
                "section, opacity, or rate is inferred"
            ),
            "photon_decay": (
                "kinematically impossible because every LI lepton pair with total "
                "momentum k has energy strictly greater than |k|"
            ),
            "physical_photon_attachment_proved": False,
            "premises": [
                "the spatial symbol is identified with physical photon frequency squared",
                (
                    "the exact cosine-symbol inequality lambda_a(k) <= "
                    "|k|_Euclidean^2 holds for a!=0"
                ),
                (
                    "electron and positron use standard positive-energy "
                    "relativistic dispersion"
                ),
                "energy and momentum compose additively in the declared frame",
            ],
            "rank_six_direction_restriction_used": False,
            "standard_lepton_premise_source_selected": False,
            "subluminal_frequency": "Omega_gamma(k,n) <= |k|",
        },
        "exact Lorentz-invariant-lepton boundary drift",
    )
    exposure = receipt["exposure_boundary"]
    check(
        exposure
        == {
            "comparison_inputs": [],
            "comparison_permitted": False,
            "evidence_claimed": False,
            "existing_Auger_limit_read": False,
            "future_event_rows_read": False,
            "public_measurement_read": False,
            "score_emitted": False,
            "verdict_emitted": False,
        },
        "exposure boundary drift",
    )
    controls = receipt["negative_controls"]
    check(
        controls
        == {
            "all_fired": True,
            "controls": [
                "replace the head-on 4 epsilon E term by 2 epsilon E",
                "replace the outgoing share cubes by squares",
                "move the standard-lepton residual maximizer away from equal sharing",
                "erase either independent charged-lepton coefficient",
                "promote the additive composition law to an A1-A3 theorem",
                (
                    "promote the universal-principal-symbol branch to "
                    "source-selected"
                ),
                "use the exact photon remainder above q=1",
                "interpret the rank-one observable map as three identified coefficients",
            ],
        },
        "negative-control contract drift",
    )
    check(
        receipt["open_physical_attachments"]
        == {
            "FZ12_oscillator_identified_with_photon": False,
            "atmospheric_shower_response_selected": False,
            "background_transport_selected": False,
            "detector_classification_selected": False,
            "electron_action_selected": False,
            "energy_momentum_composition_selected": False,
            "pair_vertex_and_cross_section_selected": False,
            "positron_action_selected": False,
            "prospective_contract_armed": False,
            "source_spectrum_and_composition_selected": False,
            "target_free_power_envelope_complete": False,
        },
        "open-attachment promotion drift",
    )
    check(raw == canonical_json_bytes(receipt), "receipt bytes changed during verification")
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args(argv)
    try:
        receipt = verify(args.receipt)
    except (OSError, KeyError, TypeError, VerificationError) as error:
        print(f"FZ12_JOINT_THRESHOLD_INVALID: {error}")
        return 1
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
