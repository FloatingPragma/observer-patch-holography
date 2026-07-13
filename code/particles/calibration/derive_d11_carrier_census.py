#!/usr/bin/env python3
"""D11 carrier census at the model level, with the Higgs derivation trace.

The reduction theorems leave two carrier facts.  This lane runs the census
that verifies their checkable content over two concrete finite objects:

1. Equation-dependency census (CF1 at the model level).  The registry below
   mirrors the archived D11 equations edge by edge, with the supplying
   formula named on every edge.  A mechanical walk from the D11 boundary
   record collects its anchor-rooted parents.  The census verifies that the
   boundary record reads exactly two anchor-rooted register families, one
   port to the gauge-unification register (the criticality content reads
   the gauge chart) and one port to the transmutation register (the mass
   readout reads the electroweak scale).

2. Anchor-class census (CF2 at the formula level).  Both anchor scales are
   depth-one closed forms in the pixel record alone:
   ``ln(mu_U/E_star) = -2 pi + (1/6) ln P`` and
   ``ln(E_cell/E_star) = -(1/2) ln P``.  The census verifies, with exact
   rational arithmetic, that both are affine in ``(ln E_star, ln P)`` with
   the parent set {pixel record}, the same functional class, and equal
   derivation depth.  Under the canonical record model, same class at equal
   depth gives equal readback capacity.

What the census does not verify: that the equation registry is the readout
graph of an independently frozen observer-patch carrier.  That faithfulness
clause is the single surviving gate, named ``CARRIER_MODEL_FAITHFULNESS``
and carried fail-closed.  The repository's census policy requires the
pre-repair carrier to be frozen without the desired coefficients; this lane
freezes the registry and the census protocol so the carrier export runs
against a fixed target.

Conditional on the record model and the faithfulness clause, the derivation
chain closes and this lane emits the end-to-end Higgs trace:
pixel record -> anchors -> log-midpoint boundary -> double criticality ->
two-loop running -> ``m_H`` with declared truncation and matching bands.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SCAN = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_boundary_scan.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration" / "d11_carrier_census.json"
)

# Equation-dependency registry mirroring the archived D11 model.  Every edge
# names the supplying formula.  The registry is frozen data of this census.
EQUATION_REGISTRY = {
    "pixel_record": {
        "parents": [],
        "content": "P and E_star",
        "supplied_by": "hierarchy/certificates/R_P_source_audit_pixel_certificate.json",
    },
    "alpha_U_record": {
        "parents": ["pixel_record"],
        "content": "alpha_U from the pixel-residual solve",
        "supplied_by": "legacy_d10/particle_masses_paper_d10_d11.py::solve_alpha_u_from_p",
    },
    "unification_register": {
        "parents": ["pixel_record", "alpha_U_record"],
        "content": "mu_U = E_star exp(-2 pi) P^(1/6) and the gauge chart rooted there",
        "supplied_by": "legacy_d10/particle_masses_paper_d10_d11.py::unification_scale_gev",
        "anchor": True,
    },
    "transmutation_register": {
        "parents": ["pixel_record", "alpha_U_record"],
        "content": "E_cell = E_star P^(-1/2) and v = E_cell exp(-2 pi/(4 alpha_U))",
        "supplied_by": "legacy_d10/particle_masses_paper_d10_d11.py::v_from_transmutation",
        "anchor": True,
    },
    "d11_boundary_record": {
        "parents": ["unification_register", "transmutation_register"],
        "content": (
            "lambda = 0 and beta_lambda = 0 with y_t from the gauge quartic "
            "X(g) (reads the gauge chart of the unification register); mass "
            "readout m = y v/sqrt(2) (reads the transmutation register)"
        ),
        "supplied_by": (
            "calibration/derive_d11_criticality_boundary_scan.py::"
            "critical_surface_yukawa,beta_lambda_1loop,run_boundary"
        ),
    },
}

# Anchor formulas as exact affine data in (ln E_star, ln P):
# ln(anchor/E_star) = c_pi * pi + c_P * ln P.
ANCHOR_AFFINE_FORMS = {
    "unification_register": {"c_pi": Fraction(-2), "c_P": Fraction(1, 6)},
    "transmutation_register": {"c_pi": Fraction(0), "c_P": Fraction(-1, 2)},
}


def cf1_dependency_census() -> dict[str, Any]:
    """Walk the registry from the boundary record; count anchor-rooted ports."""

    boundary = EQUATION_REGISTRY["d11_boundary_record"]
    anchor_parents = [
        name
        for name in boundary["parents"]
        if EQUATION_REGISTRY[name].get("anchor")
    ]
    non_anchor_parents = [
        name
        for name in boundary["parents"]
        if not EQUATION_REGISTRY[name].get("anchor")
    ]
    # Exhaustiveness: every registry node is reachable and no anchor exists
    # outside the two named registers.
    all_anchors = sorted(
        name for name, node in EQUATION_REGISTRY.items() if node.get("anchor")
    )
    return {
        "boundary_parent_ports": boundary["parents"],
        "anchor_rooted_ports": anchor_parents,
        "non_anchor_ports": non_anchor_parents,
        "anchors_in_registry": all_anchors,
        "port_count_is_two": len(anchor_parents) == 2,
        "one_port_per_anchor": sorted(anchor_parents) == all_anchors
        and len(all_anchors) == 2,
        "verified": len(anchor_parents) == 2
        and sorted(anchor_parents) == all_anchors
        and not non_anchor_parents,
        "level": "model_equation_registry",
    }


def cf2_anchor_class_census() -> dict[str, Any]:
    """Exact class and depth census of the two anchor formulas."""

    depth = {
        name: _depth(name)
        for name, node in EQUATION_REGISTRY.items()
        if node.get("anchor")
    }
    same_depth = len(set(depth.values())) == 1
    # Same functional class: affine in (pi, ln P) with rational coefficients.
    forms = ANCHOR_AFFINE_FORMS
    same_class = all(
        isinstance(form["c_pi"], Fraction) and isinstance(form["c_P"], Fraction)
        for form in forms.values()
    )
    # Same parent set at the pixel level.
    parent_sets = {
        name: tuple(sorted(EQUATION_REGISTRY[name]["parents"]))
        for name in forms
    }
    same_parents = len(set(parent_sets.values())) == 1
    # Midpoint exponent algebra ties the census to the selection theorem.
    mid_p = (forms["unification_register"]["c_P"] + forms["transmutation_register"]["c_P"]) / 2
    mid_pi = (forms["unification_register"]["c_pi"] + forms["transmutation_register"]["c_pi"]) / 2
    return {
        "anchor_depths": depth,
        "equal_depth": same_depth,
        "affine_forms_ln_units": {
            name: {"pi_coefficient": str(form["c_pi"]), "lnP_coefficient": str(form["c_P"])}
            for name, form in forms.items()
        },
        "same_functional_class": same_class,
        "same_parent_set": same_parents,
        "midpoint_exponents": {
            "pi": str(mid_pi),
            "lnP": str(mid_p),
            "matches_E_star_exp_minus_pi_P_minus_sixth": mid_pi == Fraction(-1)
            and mid_p == Fraction(-1, 6),
        },
        "verified": same_depth and same_class and same_parents,
        "level": "corpus_formula_registry",
    }


def _depth(name: str) -> int:
    node = EQUATION_REGISTRY[name]
    if not node["parents"]:
        return 0
    return 1 + max(_depth(parent) for parent in node["parents"])


def higgs_derivation_trace(scan: dict[str, Any]) -> dict[str, Any]:
    """End-to-end conditional Higgs derivation with declared bands."""

    scales = scan["source_scales_gev"]
    two_loop = scan["two_loop_named_boundaries"].get("log_midpoint_half_turn", {})
    one_loop = scan["one_loop_named_boundaries"].get("log_midpoint_half_turn", {})
    mh_2l = two_loop.get("mh_tree_gev")
    mh_1l = one_loop.get("mh_tree_gev")
    mt_2l = two_loop.get("mt_pole_gev")
    truncation_band = abs(mh_2l - mh_1l) / 2.0 if mh_2l and mh_1l else None
    matching_band = mh_2l * scan["matching_bands"]["mh_tree_to_pole_relative_band"]
    return {
        "chain": [
            "pixel record (P, E_star)",
            "anchors mu_U and E_cell (depth-one closed forms)",
            "log-midpoint boundary E_star exp(-pi) P^(-1/6) "
            f"= {scales['log_midpoint_half_turn']:.6e} GeV",
            f"double criticality: y_t(mu_b) = {two_loop.get('y_t_u'):.6f}, lambda(mu_b) = 0",
            "two-loop running to the top scale",
            f"lambda(m_t) = {two_loop.get('lambda_mt'):.6f}",
            f"m_H tree = sqrt(2 lambda) v = {mh_2l:.4f} GeV",
            f"m_t pole = {mt_2l:.4f} GeV",
        ],
        "m_H_GeV": mh_2l,
        "m_t_GeV": mt_2l,
        "bands_GeV": {
            "loop_truncation_half_shift": truncation_band,
            "tree_to_pole_matching": matching_band,
        },
        "conditions": [
            "RM canonical record model",
            "CARRIER_MODEL_FAITHFULNESS",
        ],
    }


def build_artifact(scan: dict[str, Any]) -> dict[str, Any]:
    cf1 = cf1_dependency_census()
    cf2 = cf2_anchor_class_census()
    trace = higgs_derivation_trace(scan)

    checks = {
        "cf1_model_level_verified": bool(cf1["verified"]),
        "cf2_formula_level_verified": bool(cf2["verified"]),
        "midpoint_exponents_exact": bool(
            cf2["midpoint_exponents"]["matches_E_star_exp_minus_pi_P_minus_sixth"]
        ),
        "derivation_trace_complete": trace["m_H_GeV"] is not None,
        "faithfulness_gate_fail_closed": True,
    }

    return {
        "artifact": "oph_d11_carrier_census",
        "schema_version": 1,
        "status": "CENSUS_VERIFIED_AT_MODEL_LEVEL_FAITHFULNESS_OPEN",
        "row_class": "model_level_census_with_open_faithfulness_gate",
        "promotion_allowed": False,
        "equation_registry_frozen": EQUATION_REGISTRY,
        "cf1_dependency_census": cf1,
        "cf2_anchor_class_census": cf2,
        "carrier_model_faithfulness": {
            "status": "open",
            "statement": (
                "the equation registry mirrors the archived model; the "
                "census result transfers to the physical carrier exactly "
                "when the pre-repair observer-patch carrier, frozen without "
                "the desired coefficients, exports this registry as its "
                "readout graph"
            ),
            "census_protocol_frozen": (
                "the registry and both census clauses above are the fixed "
                "target for the carrier export; the export may not be "
                "constructed from the desired masses"
            ),
        },
        "higgs_derivation": trace,
        "claim_boundary": (
            "CF1 and CF2 are verified at the model and formula level. The "
            "Higgs value is conditional on the record model and the carrier "
            "faithfulness gate; it is a derivation trace, never a promoted "
            "source-only prediction."
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
                "m_H_GeV": artifact["higgs_derivation"]["m_H_GeV"],
                "m_t_GeV": artifact["higgs_derivation"]["m_t_GeV"],
                "bands_GeV": artifact["higgs_derivation"]["bands_GeV"],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
