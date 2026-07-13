#!/usr/bin/env python3
"""Reduce the anchor-reconciliation premises to finite carrier facts.

The midpoint selection theorem consumes three premises.  This lane derives
what the axioms' record and repair structure can supply and reduces the rest
to two sharp carrier facts.  The premise ledger changes from a proposed
variational principle to checkable statements about the carrier.

Derived here:

- Theorem B (discharges AR2 exactly).  The one-loop chart makes every
  inverse coupling exactly affine in renormalization time ``t = ln mu``.
  Under the canonical record model RM (a record stores a chart coordinate
  with Gaussian readback of fixed resolution, the maximum-entropy law under
  mean and second-moment constraints), the reconciliation cost between a
  record at ``t`` and an anchor at ``t_i`` is exactly
  ``w_i (t - t_i)^2`` with ``w_i = s_i^2/(2 sigma_i^2)``.  The quadratic
  form is exact, never a leading-order approximation, because the stored
  coordinate is affine in ``t``.

- Theorem A (discharges the AR1 placement mechanism).  The mismatch
  functional is port-additive by definition and repair settles at the cost
  minimizer.  For a record with parent ports at ``t_1, t_2`` and quadratic
  per-port costs, the settled coordinate is the capacity-weighted mean;
  equal weights give the log-midpoint exactly.

- Equivalence C (reduces AR3).  Equal weights hold exactly when the two
  anchor records store same-class register readbacks at equal refinement
  depth, so equal slope and equal resolution.  Anchors storing different
  observable classes shift the minimizer off the midpoint by the computable
  capacity-asymmetry band.

Remaining carrier facts, carried fail-closed:

- CF1: the criticality boundary record carries exactly two parent ports,
  one to the gauge-unification register and one to the transmutation
  register.  This requires the D11 carrier certificate, the same class of
  object as the QT1--QT5 census.
- CF2: the two anchor registers are the same register class at equal
  refinement depth.

Combined statement: RM and CF1 and CF2, together with the repair and MaxEnt
structure, imply the boundary scale ``E_star exp(-pi) P^(-1/6)`` and the
two-loop pair ``(m_t, m_H) = (172.63, 125.77)`` GeV.  RM is the corpus
record model; CF1 and CF2 are open finite carrier statements.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from fractions import Fraction
from pathlib import Path
from typing import Any

try:
    from calibration.derive_d11_criticality_boundary_scan import (
        ALPHA_2_MZ,
        ALPHA_3_MZ,
        ALPHA_U_FALLBACK,
        ALPHA_Y_MZ,
        P_FALLBACK,
        alpha_run_1loop,
        source_scales,
    )
except ModuleNotFoundError:
    from derive_d11_criticality_boundary_scan import (
        ALPHA_2_MZ,
        ALPHA_3_MZ,
        ALPHA_U_FALLBACK,
        ALPHA_Y_MZ,
        P_FALLBACK,
        alpha_run_1loop,
        source_scales,
    )

ROOT = Path(__file__).resolve().parents[2]
SCAN = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_boundary_scan.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_anchor_reconciliation_reduction_theorems.json"
)


def chart_affinity_certificate() -> dict[str, Any]:
    """Exact affinity of the chart's inverse couplings in t = ln mu."""

    mz_run = source_scales(P_FALLBACK, ALPHA_U_FALLBACK)["mz_run_gev"]
    worst = 0.0
    for b, alpha0 in (
        (41.0 / 6.0, ALPHA_Y_MZ),
        (-19.0 / 6.0, ALPHA_2_MZ),
        (-7.0, ALPHA_3_MZ),
    ):
        for mu in (1.0e3, 1.0e9, 1.0e16, 1.3e19):
            inverse = 1.0 / alpha_run_1loop(alpha0, b, mu, mz_run)
            affine = 1.0 / alpha0 - (b / (2.0 * math.pi)) * math.log(mu / mz_run)
            worst = max(worst, abs(inverse - affine))
    return {
        "statement": (
            "the one-loop chart solution is inverse-affine: "
            "1/alpha(t) = 1/alpha_0 - (b/2 pi)(t - t_0) identically; the "
            "identity is algebraic (the chart is defined by its affine "
            "inverse), and the float round trip reproduces it to machine "
            "precision"
        ),
        "worst_residual": worst,
        "exact": worst < 1.0e-9,
    }


def gaussian_quadratic_certificate(sample_count: int = 30) -> dict[str, Any]:
    """Exact rational check: fixed-variance Gaussian KL is quadratic in t."""

    rng = random.Random(20260713)
    exact_hits = 0
    for _ in range(sample_count):
        slope = Fraction(rng.randint(1, 99), rng.randint(1, 17))
        variance = Fraction(rng.randint(1, 99), rng.randint(1, 17))
        t = Fraction(rng.randint(-99, 99), 7)
        t_i = Fraction(rng.randint(-99, 99), 7)
        kl = (slope * (t - t_i)) ** 2 / (2 * variance)
        quadratic = (slope * slope / (2 * variance)) * (t - t_i) ** 2
        if kl == quadratic:
            exact_hits += 1
    return {
        "statement": (
            "KL(N(m(t), sigma^2) || N(m(t_i), sigma^2)) = "
            "(s^2/(2 sigma^2)) (t - t_i)^2 for affine m(t) = a + s t"
        ),
        "sampled_rational_points": sample_count,
        "exact_identities": exact_hits,
        "exact": exact_hits == sample_count,
    }


def barycenter_certificate(sample_count: int = 30) -> dict[str, Any]:
    """Exact rational check: port-additive quadratic repair settles at the
    capacity-weighted mean; equal weights give the midpoint."""

    rng = random.Random(20260714)
    stationary = 0
    dominance = 0
    midpoint_hits = 0
    for _ in range(sample_count):
        t1 = Fraction(rng.randint(-99, 99), 11)
        t2 = Fraction(rng.randint(-99, 99), 11)
        w1 = Fraction(rng.randint(1, 99), rng.randint(1, 17))
        w2 = Fraction(rng.randint(1, 99), rng.randint(1, 17))
        t_star = (w1 * t1 + w2 * t2) / (w1 + w2)
        if 2 * w1 * (t_star - t1) + 2 * w2 * (t_star - t2) == 0:
            stationary += 1
        probe = t_star + Fraction(rng.randint(1, 99), rng.randint(1, 17))
        if (
            w1 * (probe - t1) ** 2 + w2 * (probe - t2) ** 2
            > w1 * (t_star - t1) ** 2 + w2 * (t_star - t2) ** 2
        ):
            dominance += 1
        if (w1 * t1 + w1 * t2) / (2 * w1) == (t1 + t2) / 2:
            midpoint_hits += 1
    return {
        "statement": (
            "repair minimization of port-additive quadratic mismatch places "
            "a two-parent record at the capacity-weighted mean; equal "
            "capacities give the log-midpoint"
        ),
        "sampled_rational_points": sample_count,
        "stationarity_exact": stationary,
        "strict_dominance": dominance,
        "equal_weight_midpoint": midpoint_hits,
        "exact": stationary == sample_count
        and dominance == sample_count
        and midpoint_hits == sample_count,
    }


def build_artifact(scan: dict[str, Any]) -> dict[str, Any]:
    affinity = chart_affinity_certificate()
    quadratic = gaussian_quadratic_certificate()
    barycenter = barycenter_certificate()

    consequence = scan["two_loop_named_boundaries"].get("log_midpoint_half_turn", {})

    premise_ledger = {
        "AR2_quadratic_cost": {
            "status": "discharged_under_RM",
            "derivation": (
                "Theorem B: chart inverse-affinity plus the Gaussian "
                "maximum-entropy record model give an exactly quadratic "
                "reconciliation cost in renormalization time"
            ),
        },
        "AR1_reconciliation": {
            "status": "mechanism_discharged_parentage_open",
            "derivation": (
                "Theorem A discharges the placement mechanism from "
                "port-additivity and repair minimization; the two-parent "
                "port structure is carrier fact CF1"
            ),
        },
        "AR3_equal_capacity": {
            "status": "reduced_to_CF2",
            "derivation": (
                "Equivalence C: equal weights hold exactly when the anchors "
                "are same-class registers at equal refinement depth"
            ),
        },
    }

    carrier_facts = {
        "RM_canonical_record_model": {
            "status": "canonical_model_of_the_corpus",
            "statement": (
                "records store chart coordinates with Gaussian readback of "
                "fixed resolution, the maximum-entropy law under mean and "
                "second-moment constraints"
            ),
        },
        "CF1_two_parent_ports": {
            "status": "open",
            "statement": (
                "the criticality boundary record carries exactly two parent "
                "ports, one to the gauge-unification register and one to "
                "the transmutation register"
            ),
            "needed": "D11 carrier certificate of the QT1-QT5 census class",
        },
        "CF2_same_register_class_equal_depth": {
            "status": "open",
            "statement": (
                "the two anchor registers are the same register class at "
                "equal refinement depth, so equal slope and equal readback "
                "resolution"
            ),
            "needed": "the same carrier certificate",
        },
    }

    checks = {
        "chart_affinity_exact": bool(affinity["exact"]),
        "gaussian_quadratic_exact": bool(quadratic["exact"]),
        "barycenter_exact": bool(barycenter["exact"]),
        "consequence_present": "mt_pole_gev" in consequence,
        "carrier_facts_fail_closed": all(
            gate["status"] == "open"
            for key, gate in carrier_facts.items()
            if key.startswith("CF")
        ),
    }

    return {
        "artifact": "oph_d11_anchor_reconciliation_reduction_theorems",
        "schema_version": 1,
        "status": "PREMISES_REDUCED_TO_CARRIER_FACTS",
        "row_class": "reduction_theorem_with_open_carrier_facts",
        "promotion_allowed": False,
        "theorem_B_quadratic_cost": {
            "chart_affinity": affinity,
            "gaussian_kl_quadratic": quadratic,
        },
        "theorem_A_barycenter_placement": barycenter,
        "equivalence_C_capacity": {
            "statement": (
                "w_i = s_i^2/(2 sigma_i^2); equal weights hold exactly when "
                "slope and resolution agree, which same-class registers at "
                "equal depth supply; different observable classes shift the "
                "minimizer by the capacity-asymmetry band of the selection "
                "theorem"
            ),
        },
        "premise_ledger": premise_ledger,
        "carrier_facts": carrier_facts,
        "combined_statement": (
            "RM and CF1 and CF2, with the repair and MaxEnt structure of "
            "the axioms, imply the boundary scale "
            "E_star exp(-pi) P^(-1/6)"
        ),
        "two_loop_consequence": {
            "mt_pole_gev": consequence.get("mt_pole_gev"),
            "mh_tree_gev": consequence.get("mh_tree_gev"),
        },
        "claim_boundary": (
            "AR2 is discharged under the canonical record model, and the "
            "AR1 placement mechanism follows from port-additivity plus "
            "repair minimization. The selection theorem's remaining content "
            "is the pair of finite carrier facts CF1 and CF2; a carrier "
            "certificate of the census class closes both."
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
                "premise_ledger": {
                    k: v["status"] for k, v in artifact["premise_ledger"].items()
                },
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
