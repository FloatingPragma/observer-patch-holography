#!/usr/bin/env python3
"""Exact global cubic certificate for the entropic charged branch.

Conditioning the twelve-port record on a nonzero W5 norm leaves the
quadratic repair cost flat on the W5 sphere, so orbit selection falls to the
universal Taylor coefficients of the record entropy,
``-sum q ln q = const - |w|^2/2 + S3(w)/6 - S4(w)/12 + ...`` with the
port-wise power sums ``S_k``.  At leading order in the band amplitude the
selection functional is the cubic ``S3`` alone, with no chosen coefficient.

Exact statement.  A vector in W5 has six antipodal-pair values ``a_i``, each
repeated twice, with ``sum a_i=0`` and ``2 sum a_i^2=1``.  At a constrained
extremum of ``S3=2 sum a_i^3``, the Lagrange equation makes every ``a_i`` a
root of one quadratic, hence there are at most two distinct values.  Exact
enumeration of their multiplicity ``m=1,...,5`` gives stationary values
``(3-m)/sqrt(3m(6-m))``.  The unique maximum up to permutation is therefore
``2/sqrt(15)`` at multiplicity one, the C5-axis orbit; the minimum is its
antipode.  The icosahedral two-design identity then gives quadrupole spectrum
``(-2/sqrt(15),-2/sqrt(15),4/sqrt(15))`` at the maximum.  Under the separately
declared quadrupole-to-physical-log-mass attachment, that degeneracy would give
two equal charged masses at leading cubic order.  The old forty-seed projected
ascent is retained only as a redundant numerical replay of the exact result.

Scope.  The separate quartic certificate
(``flavor/entropy_w5_shape_certificate.py``) now gives an exact global-
minimizer classification on both the strict full-support domain and its
closed probability simplex.  The strict minima are C5-degenerate or the
compare-only-excluded golden orbit, and above the golden positivity bound the
strict infimum is unattained.  The closed simplex instead has a viable
zero-weight boundary branch with a continuous simple spectrum.  Thus the
quartic result is a conditional strict-support no-go, never a global packet
exclusion; regularization, higher orders, full entropy, source-amplitude, and
physical-attachment routes remain open.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

try:
    from leptons.derive_charged_w5_orbit_decision_geometry import (
        P5,
        spectrum_report,
    )
except ModuleNotFoundError:
    from derive_charged_w5_orbit_decision_geometry import P5, spectrum_report

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_entropic_branch_no_go.json"
)
QUARTIC_RECEIPT = (
    ROOT / "particles" / "runs" / "flavor" / "entropy_w5_shape_certificate.json"
)

Q5 = tuple[Fraction, Fraction]


def _q5_add(left: Q5, right: Q5) -> Q5:
    return left[0] + right[0], left[1] + right[1]


def _q5_neg(value: Q5) -> Q5:
    return -value[0], -value[1]


def _q5_mul(left: Q5, right: Q5) -> Q5:
    return (
        left[0] * right[0] + 5 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def _q5_inv(value: Q5) -> Q5:
    norm = value[0] * value[0] - 5 * value[1] * value[1]
    if norm == 0:
        raise ZeroDivisionError("zero in Q(sqrt(5))")
    return value[0] / norm, -value[1] / norm


def _q5_div(left: Q5, right: Q5) -> Q5:
    return _q5_mul(left, _q5_inv(right))


def _exact_icosahedral_axis_frame_check() -> dict[str, Any]:
    """Verify sum p_i p_i^T=2I for six normalized axes in Q(sqrt(5))."""

    zero: Q5 = (Fraction(0), Fraction(0))
    one: Q5 = (Fraction(1), Fraction(0))
    phi: Q5 = (Fraction(1, 2), Fraction(1, 2))
    minus_phi = _q5_neg(phi)
    minus_one = _q5_neg(one)
    axes: list[tuple[Q5, Q5, Q5]] = [
        (zero, one, phi),
        (zero, one, minus_phi),
        (one, phi, zero),
        (one, minus_phi, zero),
        (phi, zero, one),
        (phi, zero, minus_one),
    ]
    norms = [
        _q5_add(_q5_add(_q5_mul(x, x), _q5_mul(y, y)), _q5_mul(z, z))
        for x, y, z in axes
    ]
    common_norm = norms[0]
    equal_norms = all(norm == common_norm for norm in norms)
    outer_sum: list[list[Q5]] = [[zero for _ in range(3)] for _ in range(3)]
    for axis in axes:
        for row in range(3):
            for column in range(3):
                outer_sum[row][column] = _q5_add(
                    outer_sum[row][column],
                    _q5_mul(axis[row], axis[column]),
                )
    normalized = [
        [_q5_div(value, common_norm) for value in row]
        for row in outer_sum
    ]
    target = [
        [(Fraction(2), Fraction(0)) if row == column else zero for column in range(3)]
        for row in range(3)
    ]
    identity_verified = equal_norms and normalized == target
    return {
        "field": "Q(sqrt(5)) exact arithmetic",
        "axis_representatives": (
            "(0,1,+/-phi),(1,+/-phi,0),(phi,0,+/-1), normalized by "
            "sqrt(1+phi^2)"
        ),
        "common_unnormalized_norm_squared": "1+phi^2",
        "all_six_norms_equal": equal_norms,
        "normalized_outer_product_sum": [
            ["2" if row == column else "0" for column in range(3)]
            for row in range(3)
        ],
        "sum_over_six_axes_pipt_equals_2I": identity_verified,
    }


def maximize_s3(seeds: int = 40, iters: int = 4000) -> tuple[float, np.ndarray]:
    best: tuple[float, np.ndarray] | None = None
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        w = P5 @ rng.standard_normal(12)
        w /= np.linalg.norm(w)
        for _ in range(iters):
            grad = P5 @ (3.0 * w**2)
            grad -= (grad @ w) * w
            if np.linalg.norm(grad) < 1.0e-13:
                break
            w = w + 0.05 * grad
            w /= np.linalg.norm(w)
        value = float(np.sum(w**3))
        if best is None or value > best[0]:
            best = (value, w.copy())
    assert best is not None
    return best


def exact_cubic_certificate() -> dict[str, Any]:
    """Return the exhaustive two-value stationary-point enumeration."""

    stationary = []
    for multiplicity in range(1, 6):
        other_multiplicity = 6 - multiplicity
        b_over_a = Fraction(-multiplicity, other_multiplicity)
        a_squared = Fraction(other_multiplicity, 12 * multiplicity)
        linear_constraint = (
            Fraction(multiplicity)
            + Fraction(other_multiplicity) * b_over_a
        )
        norm_constraint = (
            2
            * (
                Fraction(multiplicity)
                + Fraction(other_multiplicity) * b_over_a * b_over_a
            )
            * a_squared
        )
        objective_coefficient_in_a_cubed = 2 * (
            Fraction(multiplicity)
            + Fraction(other_multiplicity) * b_over_a**3
        )
        derived_objective_squared = (
            objective_coefficient_in_a_cubed**2 * a_squared**3
        )
        numerator = 3 - multiplicity
        denominator_squared = 3 * multiplicity * (6 - multiplicity)
        objective_squared = Fraction(numerator * numerator, denominator_squared)
        derived_sign = (
            "positive"
            if objective_coefficient_in_a_cubed > 0
            else "negative"
            if objective_coefficient_in_a_cubed < 0
            else "zero"
        )
        formula_sign = "positive" if numerator > 0 else "negative" if numerator < 0 else "zero"
        stationary.append(
            {
                "positive_root_multiplicity": multiplicity,
                "other_root_multiplicity": other_multiplicity,
                "other_root_over_positive_root": {
                    "numerator": b_over_a.numerator,
                    "denominator": b_over_a.denominator,
                },
                "positive_root_squared": {
                    "numerator": a_squared.numerator,
                    "denominator": a_squared.denominator,
                },
                "objective_exact": (
                    "0"
                    if numerator == 0
                    else f"{numerator}/sqrt({denominator_squared})"
                ),
                "objective_sign": derived_sign,
                "objective_squared": {
                    "numerator": derived_objective_squared.numerator,
                    "denominator": derived_objective_squared.denominator,
                },
                "objective_approximation": float(
                    numerator / np.sqrt(denominator_squared)
                ),
                "exact_checks": {
                    "zero_sum_constraint": linear_constraint == 0,
                    "unit_norm_constraint": norm_constraint == 1,
                    "stationary_formula_squared_matches": (
                        derived_objective_squared == objective_squared
                    ),
                    "stationary_formula_sign_matches": derived_sign == formula_sign,
                },
            }
        )

    positive_rows = [row for row in stationary if row["objective_sign"] == "positive"]
    maximum = max(
        positive_rows,
        key=lambda row: Fraction(
            row["objective_squared"]["numerator"],
            row["objective_squared"]["denominator"],
        ),
    )
    exact_maximum_is_m1 = maximum["positive_root_multiplicity"] == 1
    frame_check = _exact_icosahedral_axis_frame_check()
    scale_squared = Fraction(1, 60)
    high_coefficient = Fraction(5)
    low_coefficient = Fraction(-1)
    m1_zero_sum = high_coefficient + 5 * low_coefficient
    m1_norm = 2 * (
        high_coefficient**2 + 5 * low_coefficient**2
    ) * scale_squared
    transverse_coefficient = 4 * low_coefficient
    axial_coefficient = -8 * low_coefficient
    transverse_squared = transverse_coefficient**2 * scale_squared
    axial_squared = axial_coefficient**2 * scale_squared
    spectrum_matches_claim = (
        transverse_coefficient < 0
        and axial_coefficient > 0
        and transverse_squared == Fraction(4, 15)
        and axial_squared == Fraction(16, 15)
        and 2 * transverse_coefficient + axial_coefficient == 0
    )
    quadrupole_double_eigenvalue_exact = (
        frame_check["sum_over_six_axes_pipt_equals_2I"]
        and m1_zero_sum == 0
        and m1_norm == 1
        and spectrum_matches_claim
    )
    return {
        "proof_kind": "exact_closed_form_global_extremum_certificate",
        "w5_pair_reduction": {
            "coordinates": "w=(a1,a1,...,a6,a6) after ordering antipodal pairs",
            "linear_constraint": "sum_i a_i = 0",
            "unit_sphere_constraint": "2*sum_i a_i^2 = 1",
            "objective": "S3(w)=2*sum_i a_i^3",
            "justification": (
                "range(P5) is the antipodal-even subspace with its constant "
                "line removed"
            ),
        },
        "completeness_argument": {
            "compactness": "the constrained W5 unit sphere is compact",
            "regularity": (
                "the gradients of sum_i a_i and sum_i a_i^2 are independent "
                "on the nonzero zero-sum sphere"
            ),
            "lagrange_equation": "6*a_i^2 = lambda + 4*mu*a_i",
            "root_bound": (
                "all six a_i are roots of one quadratic, hence every "
                "stationary point has at most two distinct values"
            ),
            "multiplicities_exhausted": [1, 2, 3, 4, 5],
        },
        "stationary_value_formula": "S3_m=(3-m)/sqrt(3*m*(6-m))",
        "stationary_enumeration": stationary,
        "global_maximum": {
            "multiplicity": 1,
            "value_exact": "2/sqrt(15)",
            "value_squared_exact": {"numerator": 4, "denominator": 15},
            "unique_up_to_pair_permutation": True,
            "a_high_squared": {"numerator": 5, "denominator": 12},
            "a_low_squared": {"numerator": 1, "denominator": 60},
            "root_relation": "a_high=-5*a_low with a_high>0",
            "orbit": "one icosahedral antipodal vertex pair versus the other five; C5 axis",
        },
        "global_minimum": {
            "multiplicity": 5,
            "value_exact": "-2/sqrt(15)",
            "relation": "antipode of the global maximum",
        },
        "quadrupole_degeneracy": {
            "axis_frame_check": frame_check,
            "m1_pair_coefficients_in_units_of_1_over_sqrt60": {
                "a_high": int(high_coefficient),
                "a_low": int(low_coefficient),
                "zero_sum": m1_zero_sum == 0,
                "unit_norm": m1_norm == 1,
            },
            "derivation": (
                "for pair values A and B with A+5B=0, Q has eigenvalues "
                "4B,4B,-8B"
            ),
            "spectrum_coefficients_in_units_of_1_over_sqrt60": [
                int(transverse_coefficient),
                int(transverse_coefficient),
                int(axial_coefficient),
            ],
            "spectrum_squared_checks": {
                "transverse_equals_4_over_15": transverse_squared == Fraction(4, 15),
                "axial_equals_16_over_15": axial_squared == Fraction(16, 15),
                "signs_and_trace_match": spectrum_matches_claim,
            },
            "maximum_spectrum_exact": [
                "-2/sqrt(15)",
                "-2/sqrt(15)",
                "4/sqrt(15)",
            ],
            "double_eigenvalue_exact": quadrupole_double_eigenvalue_exact,
        },
        "checks": {
            "all_stationary_multiplicities_enumerated": len(stationary) == 5,
            "all_pair_constraints_and_stationary_formulas_verified": all(
                all(row["exact_checks"].values()) for row in stationary
            ),
            "global_maximum_is_multiplicity_one": exact_maximum_is_m1,
            "m1_value_squared_is_4_over_15": Fraction(
                maximum["objective_squared"]["numerator"],
                maximum["objective_squared"]["denominator"],
            )
            == Fraction(4, 15),
            "every_other_positive_stationary_value_is_smaller": all(
                Fraction(
                    row["objective_squared"]["numerator"],
                    row["objective_squared"]["denominator"],
                )
                < Fraction(4, 15)
                for row in positive_rows
                if row["positive_root_multiplicity"] != 1
            ),
            "icosahedral_axis_frame_identity_exact": frame_check[
                "sum_over_six_axes_pipt_equals_2I"
            ],
            "quadrupole_double_eigenvalue_exact": quadrupole_double_eigenvalue_exact,
        },
    }


def build() -> dict[str, Any]:
    quartic = json.loads(QUARTIC_RECEIPT.read_text(encoding="utf-8"))
    quartic_boundary = quartic.get("exhaustiveness_boundary", {})
    quartic_global_minimizers_proved = bool(
        quartic.get("schema") == "oph.entropy_w5_shape_certificate.v3"
        and quartic_boundary.get("global_minimizer_classification_proved") is True
        and quartic_boundary.get("quartic_packet_globally_excluded") is False
        and quartic.get("finite_seed_or_numerical_optimization_used_as_proof") is False
    )
    if not quartic_global_minimizers_proved:
        raise ValueError("the exact quartic W5 certificate is absent or fail-closed")
    exact = exact_cubic_certificate()
    value, w = maximize_s3()
    report = spectrum_report(w)
    eigenvalues = report["eigenvalues"]
    min_gap = min(
        eigenvalues[1] - eigenvalues[0], eigenvalues[2] - eigenvalues[1]
    )
    degenerate = min_gap < 1.0e-8
    # S3 is odd, so the antipode has the negated candidate objective value.
    antipodal_eigenvalues = spectrum_report(-w)["eigenvalues"]
    antipodal_value = float(np.sum((-w) ** 3))
    antipodal_gap = min(
        antipodal_eigenvalues[1] - antipodal_eigenvalues[0],
        antipodal_eigenvalues[2] - antipodal_eigenvalues[1],
    )
    numerical_checks = {
        "extremum_found": value > 0.5,
        "extremal_spectrum_degenerate": bool(degenerate),
        "antipodal_value_is_negated": abs(antipodal_value + value) < 1.0e-12,
        "antipodal_spectrum_degenerate": bool(antipodal_gap < 1.0e-8),
        "c5_candidate_reproduced": bool(
            abs(value - 2.0 / np.sqrt(15.0)) < 1.0e-10
        ),
    }
    return {
        "artifact": "oph_charged_entropic_branch_no_go",
        "schema_version": 3,
        "status": "ENTROPIC_CUBIC_AND_QUARTIC_EXACT_GLOBAL_CLASSIFICATION",
        "row_class": "exact_cubic_extremum_and_quartic_global_minimizer_theorems",
        "promotion_allowed": False,
        "cubic_truncation_theorem_established": True,
        "quartic_global_minimizer_theorem_established": quartic_global_minimizers_proved,
        "quartic_packet_globally_excluded": False,
        "legacy_filename_note": (
            "the filename predates the exact certificates; the cubic no-go "
            "is exact, while the quartic no-go is conditional on strict full "
            "support and the physical mass attachment"
        ),
        "selection_functional": (
            "maximize the universal entropy cubic S3 on the W5 unit sphere, "
            "the leading shape-dependent term in the band amplitude; no "
            "chosen coefficient"
        ),
        "exact_cubic_certificate": exact,
        "global_maximum_value_s3_exact": "2/sqrt(15)",
        "global_maximum_spectrum_exact": [
            "-2/sqrt(15)",
            "-2/sqrt(15)",
            "4/sqrt(15)",
        ],
        "oddness": (
            "S3(-w) = -S3(w); after the exhaustive stationary-point "
            "enumeration proves the global maximum, oddness proves that its "
            "antipode is the global minimum"
        ),
        "consequence": (
            "the exact global maximizer of the leading cubic entropy "
            "truncation is C5-axis and has a double quadrupole eigenvalue; "
            "two equal physical charged masses follow only conditional on "
            "the quadrupole-to-physical-log-mass attachment. The exact "
            "quartic theorem gives only a strict-full-support conditional "
            "no-go because its closed-simplex boundary branch remains viable; "
            "no higher-order or full-mechanism no-go follows"
        ),
        "epistemic_gates": {
            "cubic_global_extremum_proved": True,
            "cubic_no_go_certified": True,
            "finite_seed_search_can_close_route": False,
            "quartic_global_optimality_proved": quartic_global_minimizers_proved,
            "strict_full_support_quartic_no_go_conditional": (
                quartic_global_minimizers_proved
            ),
            "quartic_packet_globally_excluded": False,
            "full_entropic_mechanism_no_go_certified": False,
            "quadrupole_to_physical_log_mass_attachment_established": False,
        },
        "numerical_replay": {
            "role": "redundant_diagnostic_only_not_used_in_proof",
            "seeds": 40,
            "projected_ascent_iterations_per_seed": 4000,
            "value_s3": value,
            "spectrum": eigenvalues,
            "minimum_gap": min_gap,
            "antipodal_value_s3": antipodal_value,
            "antipodal_spectrum": antipodal_eigenvalues,
            "checks": numerical_checks,
        },
        "scope": {
            "truncation_order": "leading cubic term at small band amplitude",
            "quartic_order": (
                "the exact companion certificate classifies global minima on "
                "the strict and closed probability domains. The strict branch "
                "is conditionally excluded after the separate mass attachment, "
                "but a zero-weight closed-simplex boundary branch with a "
                "continuous simple spectrum remains viable"
            ),
            "quartic_certificate": "runs/flavor/entropy_w5_shape_certificate.json",
            "beyond_quartic": "open",
            "mechanism_status": (
                "exact cubic and quartic minimizer theorems; quartic strict-"
                "support no-go remains conditional, the packet is not globally "
                "excluded, and boundary regularization, higher orders, source "
                "amplitude, physical attachment, and full entropy remain open"
            ),
            "universal_impossibility_claimed": False,
            "physical_mass_attachment": (
                "open: the exact statement is eigenvalue degeneracy on the "
                "W5 quadrupole; identification with physical charged "
                "log-masses is a separate conditional attachment"
            ),
        },
        "checks": {
            "exact_certificate_checks_pass": all(
                bool(v) for v in exact["checks"].values()
            ),
            "numerical_replay_checks_pass": all(
                bool(v) for v in numerical_checks.values()
            ),
            "quartic_boundary_and_regularization_routes_kept_open": True,
            "quartic_global_minimizer_receipt_consumed": (
                quartic_global_minimizers_proved
            ),
            "full_mechanism_no_go_not_claimed": True,
        },
        "checks_pass": (
            all(bool(v) for v in exact["checks"].values())
            and all(bool(v) for v in numerical_checks.values())
            and quartic_global_minimizers_proved
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "numerical_replay_spectrum": artifact["numerical_replay"]["spectrum"],
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
