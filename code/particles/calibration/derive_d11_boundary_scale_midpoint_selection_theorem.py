#!/usr/bin/env python3
"""Conditional selection theorem for the criticality boundary scale.

The boundary-scale selection audit closes the flow-internal route and leaves
the scale to the source structure.  This lane states and proves the strongest
selection theorem available: a rigorous implication from three named
premises to the unique boundary scale, with the premises carried as
fail-closed gates.

Premises (anchor-reconciliation package):

- AR1 (reconciliation): the criticality boundary is a record that reconciles
  the model's two pre-existing high-scale anchor records, the
  gauge-unification record at ``mu_U = E_star exp(-2 pi) P^(1/6)`` and the
  transmutation record at ``E_cell = E_star P^(-1/2)``.
- AR2 (quadratic cost): the reconciliation cost is quadratic in the
  renormalization time ``t = ln mu`` offset from each anchor.  This is the
  mismatch form of the repair calculus and the leading form of relative
  entropy around a reference record.
- AR3 (equal capacity): the two anchor records carry equal capacity, so the
  two quadratic weights are equal.

Theorem (proved here, exactly): under AR1, AR2, AR3 the cost
``F(t) = w (t - t1)^2 + w (t - t2)^2`` is strictly convex with the unique
minimizer at the log-midpoint, so the boundary scale is
``mu_b = sqrt(mu_U E_cell) = E_star exp(-pi) P^(-1/6)``, and with unequal
capacities the minimizer is the capacity-weighted geometric mean.

The implication is exact.  The premises are a proposed selection principle:
AR1 and AR3 have no derivation from the five axioms in the current corpus,
and AR2 is motivated by, without being derived from, the quadratic mismatch
form.  The theorem therefore carries the same claim class as every completed
conditional law in this repository, and the registered three-loop
discriminating test of the selection audit applies to it unchanged.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SCAN = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_boundary_scan.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_boundary_scale_midpoint_selection_theorem.json"
)


def weighted_minimizer_exact(
    t1: Fraction, t2: Fraction, w1: Fraction, w2: Fraction
) -> Fraction:
    """Unique minimizer of w1 (t-t1)^2 + w2 (t-t2)^2, from the exact gradient."""

    return (w1 * t1 + w2 * t2) / (w1 + w2)


def implication_checks(sample_count: int = 40) -> dict[str, Any]:
    """Exact rational certificates for the implication."""

    rng = random.Random(20260713)
    stationarity_zero = 0
    convexity_positive = 0
    dominance_holds = 0
    for _ in range(sample_count):
        t1 = Fraction(rng.randint(-999, 999), rng.randint(1, 97))
        t2 = Fraction(rng.randint(-999, 999), rng.randint(1, 97))
        w1 = Fraction(rng.randint(1, 999), rng.randint(1, 97))
        w2 = Fraction(rng.randint(1, 999), rng.randint(1, 97))
        t_star = weighted_minimizer_exact(t1, t2, w1, w2)
        gradient = 2 * w1 * (t_star - t1) + 2 * w2 * (t_star - t2)
        if gradient == 0:
            stationarity_zero += 1
        if 2 * (w1 + w2) > 0:
            convexity_positive += 1
        probe = t_star + Fraction(rng.randint(1, 99), rng.randint(1, 97))
        cost_star = w1 * (t_star - t1) ** 2 + w2 * (t_star - t2) ** 2
        cost_probe = w1 * (probe - t1) ** 2 + w2 * (probe - t2) ** 2
        if cost_probe > cost_star:
            dominance_holds += 1
    equal_weight_midpoint = weighted_minimizer_exact(
        Fraction(3), Fraction(11), Fraction(5), Fraction(5)
    ) == Fraction(7)
    return {
        "sampled_rational_points": sample_count,
        "gradient_zero_at_minimizer": stationarity_zero,
        "second_derivative_positive": convexity_positive,
        "strict_dominance_off_minimizer": dominance_holds,
        "equal_weights_give_exact_midpoint": bool(equal_weight_midpoint),
        "implication_holds": stationarity_zero == sample_count
        and convexity_positive == sample_count
        and dominance_holds == sample_count
        and equal_weight_midpoint,
    }


def closed_form_identity() -> dict[str, Any]:
    """Exact exponent algebra: sqrt(mu_U E_cell) = E_star exp(-pi) P^(-1/6)."""

    # Exponents of (e^{-pi}) and P in each anchor: mu_U = E* e^{-2pi} P^{1/6},
    # E_cell = E* P^{-1/2}.  The log-midpoint halves the exponent sums.
    e_exponent = Fraction(-2, 1) * Fraction(1, 2)  # -> -1 (in units of pi)
    p_exponent = (Fraction(1, 6) + Fraction(-1, 2)) * Fraction(1, 2)
    return {
        "pi_exponent_halved": str(e_exponent),
        "p_exponent_halved": str(p_exponent),
        "identity": "sqrt(mu_U * E_cell) = E_star * exp(-pi) * P^(-1/6)",
        "exact": e_exponent == Fraction(-1) and p_exponent == Fraction(-1, 6),
    }


def capacity_asymmetry_band(scan: dict[str, Any]) -> dict[str, Any]:
    """Map a capacity-weight asymmetry onto a boundary-scale and mass band.

    With weights w1, w2 and r = w2/w1, the minimizer offset from the midpoint
    is (delta_t/2) (r - 1)/(r + 1), with delta_t = ln(E_cell/mu_U).  The
    two-loop curve moves the Higgs coordinate by about 2.1 GeV per e-fold, so
    a stated asymmetry bound translates into a mass band.
    """

    scales = scan["source_scales_gev"]
    delta_t = math.log(scales["E_cell_pixel_energy"] / scales["mu_U_gauge_unification"])
    rows = []
    for r in (1.0, 1.25, 2.0):
        offset = (delta_t / 2.0) * (r - 1.0) / (r + 1.0)
        rows.append(
            {
                "capacity_ratio": r,
                "offset_e_folds": offset,
                "approx_mh_shift_gev": 2.1 * offset,
            }
        )
    return {
        "delta_t_e_folds": delta_t,
        "rows": rows,
        "statement": (
            "AR3 is testable: a capacity asymmetry moves the selected scale "
            "off the midpoint by a computable amount, and the three-loop "
            "implied scale measures it."
        ),
    }


def build_artifact(scan: dict[str, Any]) -> dict[str, Any]:
    implication = implication_checks()
    identity = closed_form_identity()
    band = capacity_asymmetry_band(scan)

    midpoint = scan["source_scales_gev"]["log_midpoint_half_turn"]
    consequence = scan["two_loop_named_boundaries"].get("log_midpoint_half_turn", {})

    premises = {
        "AR1_reconciliation": {
            "status": "open",
            "statement": (
                "the criticality boundary reconciles the gauge-unification "
                "record and the transmutation record"
            ),
            "needed": "derivation from the record/repair structure of the axioms",
        },
        "AR2_quadratic_cost": {
            "status": "open",
            "statement": (
                "reconciliation cost is quadratic in renormalization time "
                "offset"
            ),
            "needed": (
                "promotion of the repair-calculus quadratic mismatch form to "
                "this record class"
            ),
        },
        "AR3_equal_capacity": {
            "status": "open",
            "statement": "the two anchor records carry equal capacity",
            "needed": (
                "a capacity computation for both anchor records; the "
                "three-loop implied scale measures any asymmetry"
            ),
        },
    }

    checks = {
        "implication_proved_exactly": bool(implication["implication_holds"]),
        "closed_form_identity_exact": bool(identity["exact"]),
        "consequence_rows_present": "mt_pole_gev" in consequence,
        "premises_carried_fail_closed": all(
            gate["status"] == "open" for gate in premises.values()
        ),
    }

    return {
        "artifact": "oph_d11_boundary_scale_midpoint_selection_theorem",
        "schema_version": 1,
        "status": "CONDITIONAL_SELECTION_THEOREM_PREMISES_OPEN",
        "row_class": "conditional_theorem_with_named_premises",
        "promotion_allowed": False,
        "theorem": {
            "name": "ANCHOR_RECONCILIATION_MIDPOINT_SELECTION_THEOREM",
            "hypotheses": ["AR1_reconciliation", "AR2_quadratic_cost", "AR3_equal_capacity"],
            "conclusion": (
                "the boundary scale is the log-midpoint "
                "sqrt(mu_U * E_cell) = E_star exp(-pi) P^(-1/6); with unequal "
                "capacities it is the capacity-weighted geometric mean"
            ),
            "proof": (
                "the cost is strictly convex in renormalization time with "
                "positive second derivative 2(w1+w2); the exact gradient "
                "vanishes only at the weighted mean; equal weights give the "
                "midpoint; exponent algebra gives the closed form"
            ),
        },
        "implication_certificates": implication,
        "closed_form_identity": identity,
        "premises": premises,
        "capacity_asymmetry_band": band,
        "selected_scale_gev": midpoint,
        "two_loop_consequence": {
            "mt_pole_gev": consequence.get("mt_pole_gev"),
            "mh_tree_gev": consequence.get("mh_tree_gev"),
        },
        "claim_boundary": (
            "The implication from AR1, AR2, AR3 to the midpoint scale is "
            "proved exactly. The premises are a proposed selection principle "
            "without a derivation from the five axioms; the theorem shares "
            "the claim class of the completed conditional laws and inherits "
            "the registered three-loop discriminating test."
        ),
        "checks": checks,
        "checks_pass": all(bool(v) for v in checks.values()),
    }


def build() -> dict[str, Any]:
    scan = json.loads(SCAN.read_text(encoding="utf-8"))
    return build_artifact(scan)


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
                "selected_scale_gev": artifact["selected_scale_gev"],
                "two_loop_consequence": artifact["two_loop_consequence"],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
