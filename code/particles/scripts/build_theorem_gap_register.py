#!/usr/bin/env python3
"""Build the theorem/proof gap register for the particle program.

One artifact enumerates every named open object across the particle
families, classifies it, and records the closures of the current pass with
their artifact pointers.  Classes:

- closed: proved or certified, artifact named.
- research_open: needs a new theorem; no known mechanical route.
- closable_by_standard_computation: standard mathematics or engineering of
  known form; effort-bound, no new principle needed.
- computationally_blocked: needs computation beyond current resources.
- awaiting_discriminating_test: frozen candidates with a registered test.

The register also carries the massless-carrier statement receipts: exact
photon masslessness from unbroken electromagnetic gauge invariance of the
recovered quotient, and carrier-level gluon masslessness under unbroken
color with the confinement caveat.  Both are statement receipts conditional
on the recovered-core gauge structure; neither is a numerical lane.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve()
PARTICLES = HERE.parents[1]
RUNS = PARTICLES / "runs"
DEFAULT_JSON_OUT = RUNS / "status" / "theorem_gap_register.json"
DEFAULT_MD_OUT = PARTICLES / "THEOREM_GAP_REGISTER.md"


def closed_entries() -> list[dict[str, Any]]:
    return [
        {
            "id": "HIGGS_DEFICIT_ATTRIBUTION",
            "family": "Higgs/top",
            "class": "closed",
            "statement": (
                "the archived (115.1, 164.1) deficit decomposes exactly into "
                "the boundary-scale choice plus loop truncation; both "
                "boundaries derive from the gauge sector through the "
                "double-criticality law"
            ),
            "artifact": "runs/calibration/d11_criticality_boundary_scan.json",
        },
        {
            "id": "FLOW_INTERNAL_SELECTION_NO_GO",
            "family": "Higgs/top",
            "class": "closed",
            "statement": (
                "triple criticality has no root: the flow derivative stays "
                "above 4.4e-5 across the window, so the boundary scale is a "
                "source-side object"
            ),
            "artifact": "runs/calibration/d11_boundary_scale_selection_audit.json",
        },
        {
            "id": "MT_MH_RELATION",
            "family": "Higgs/top",
            "class": "closed",
            "statement": (
                "the fit-free criticality curve reproduces the Higgs "
                "coordinate at the measured top to 0.47 percent at two "
                "loops, inside the declared matching band"
            ),
            "artifact": "runs/calibration/d11_criticality_comparison.json",
        },
        {
            "id": "ANCHOR_RECONCILIATION_MIDPOINT_IMPLICATION",
            "family": "Higgs/top",
            "class": "closed",
            "statement": (
                "AR1+AR2+AR3 imply the unique boundary scale "
                "E_star exp(-pi) P^(-1/6); the implication is proved "
                "exactly, the premises are open"
            ),
            "artifact": (
                "runs/calibration/"
                "d11_boundary_scale_midpoint_selection_theorem.json"
            ),
        },
        {
            "id": "AR_PREMISE_REDUCTION",
            "family": "Higgs/top",
            "class": "closed",
            "statement": (
                "AR2 is discharged exactly under the canonical record model "
                "(chart inverse-affinity plus Gaussian MaxEnt readback), and "
                "the AR1 placement mechanism follows from port-additivity "
                "plus repair minimization; AR1 parentage and AR3 reduce to "
                "the finite carrier facts CF1 and CF2"
            ),
            "artifact": (
                "runs/calibration/"
                "d11_anchor_reconciliation_reduction_theorems.json"
            ),
        },
        {
            "id": "ENTROPIC_CONDITIONED_BRANCH_NO_GO",
            "family": "charged leptons",
            "class": "closed",
            "statement": (
                "the parameter-free entropic conditioned branch selects the "
                "C5-axis orbit with an exactly degenerate spectrum: two "
                "equal charged masses; the second of the three candidate "
                "mechanisms is closed, and the surviving route is a "
                "source-emitted charged interaction off the entropic ray"
            ),
            "artifact": "runs/leptons/charged_entropic_branch_no_go.json",
        },
        {
            "id": "W5_ORBIT_DECISION_GEOMETRY",
            "family": "charged leptons",
            "class": "closed",
            "statement": (
                "stratum classification certified (C5 and C3 axes exactly "
                "degenerate, C2 simple: three distinct masses force the "
                "orbit off the symmetry axes); the MCPR target shape lies "
                "inside the simple-spectrum region, witnessed by a frozen "
                "coefficient point with minimizing-orbit gap ratio within "
                "0.006 of the target; the harness maps any coefficient "
                "packet to the family shape with no further choices"
            ),
            "artifact": "runs/leptons/charged_w5_orbit_decision_geometry.json",
        },
        {
            "id": "UP_TYPE_INTEGER_POWER_LAW_REMOVED",
            "family": "quarks",
            "class": "closed",
            "statement": (
                "the frozen four-base integer exponent scan is negative in "
                "every base with both channels reported; the law family is "
                "removed from the candidate space prospectively"
            ),
            "artifact": "runs/flavor/up_type_register_exponent_scan.json",
        },
        {
            "id": "DOWN_TYPE_CLEBSCH_CONDITIONAL_SECTOR",
            "family": "quarks",
            "class": "closed",
            "statement": (
                "conditional on the MCPR leptons and the register Clebsch "
                "pattern (1, 1/3, 3), the down-type sector is emitted: "
                "ratios and the Gatto-Sartori-Tonin angle land at the "
                "ten-percent scale, the absolute normalization carries the "
                "named third-generation tension"
            ),
            "artifact": "runs/flavor/down_type_register_clebsch_lane.json",
        },
        {
            "id": "CF1_CF2_MODEL_LEVEL_CENSUS",
            "family": "Higgs/top",
            "class": "closed",
            "statement": (
                "the equation-dependency census verifies exactly two "
                "anchor-rooted parent ports on the D11 boundary record, and "
                "the anchor-class census certifies same functional class, "
                "same pixel parent set, and equal depth for the two anchors, "
                "with exact midpoint exponent algebra"
            ),
            "artifact": "runs/calibration/d11_carrier_census.json",
        },
        {
            "id": "CHARGED_HOMOGENEOUS_SHAPE_SILENCE",
            "family": "charged leptons",
            "class": "closed",
            "statement": (
                "a unique invariant MaxEnt state has zero non-singlet "
                "expectation and the homogeneous twelve-port branch has the "
                "uniform minimizer: family shape requires a selected W5 "
                "orbit"
            ),
            "artifact": "runs/leptons/charged_mcpr_completion_conditional.json",
        },
        {
            "id": "PHOTON_EXACT_MASSLESSNESS_RECEIPT",
            "family": "massless carriers",
            "class": "closed",
            "statement": (
                "unbroken electromagnetic gauge invariance of the recovered "
                "quotient forces a transverse self-energy with vanishing "
                "zero-momentum part; the photon mass is exactly zero, "
                "conditional on the recovered-core gauge structure"
            ),
            "artifact": "runs/status/theorem_gap_register.json",
        },
        {
            "id": "GLUON_CARRIER_MASSLESSNESS_RECEIPT",
            "family": "massless carriers",
            "class": "closed",
            "statement": (
                "unbroken color forces carrier-level gluon masslessness; "
                "confinement removes asymptotic gluon states, so the "
                "statement concerns the gauge field, never a physical pole"
            ),
            "artifact": "runs/status/theorem_gap_register.json",
        },
    ]


def open_entries() -> list[dict[str, Any]]:
    return [
        {
            "id": "CARRIER_MODEL_FAITHFULNESS",
            "family": "Higgs/top",
            "class": "research_open",
            "statement": (
                "export the pre-repair observer-patch carrier, frozen "
                "without the desired coefficients, and show that its "
                "readout graph is the frozen equation registry of the "
                "census; this transfers CF1 and CF2 from the model level "
                "to the carrier level and completes the Higgs selection "
                "theorem"
            ),
        },
        {
            "id": "CRITICALITY_BOUNDARY_SCALE_SELECTION",
            "family": "Higgs/top",
            "class": "awaiting_discriminating_test",
            "statement": (
                "four frozen candidates; the three-loop implied scale "
                "discriminates"
            ),
            "test": "FROZEN_THREE_LOOP_RG_MATCHING_PACKET",
        },
        {
            "id": "FROZEN_THREE_LOOP_RG_MATCHING_PACKET",
            "family": "Higgs/top",
            "class": "closable_by_standard_computation",
            "statement": (
                "three-loop running plus full NNLO matching, frozen "
                "target-free; standard published mathematics, "
                "implementation-bound"
            ),
        },
        {
            "id": "D10_DISCRETE_LAW_SELECTION",
            "family": "W/Z",
            "class": "research_open",
            "statement": (
                "select between the zero-selector law and the carrier value "
                "law; the QME mechanism applies once the electroweak moment "
                "vector is emitted"
            ),
        },
        {
            "id": "SM_MOMENT_VECTOR_EMISSION",
            "family": "cross-cutting",
            "class": "research_open",
            "statement": (
                "the complete quotient-local operator basis and moment "
                "vector c_r(P, N) for the source-action rigidity mechanism"
            ),
        },
        {
            "id": "W5_ORBIT_EFFECTIVE_ACTION",
            "family": "charged leptons",
            "class": "research_open",
            "statement": (
                "SOURCE_W5_COEFFICIENT_EMISSION: the charge-conditioned "
                "repair branch must emit the invariant coefficient packet "
                "(h5, a, b, c, e); the decision geometry and harness are "
                "closed, so emission alone completes the charged family "
                "shape"
            ),
        },
        {
            "id": "A5_FAMILY_LIFT_AND_DET_CAN",
            "family": "charged leptons",
            "class": "research_open",
            "statement": (
                "the physical A5 family lift with the multiplicity-one "
                "attachment, and the normed determinant-line descent with "
                "kinetic factors"
            ),
        },
        {
            "id": "QF1_QF9_FLAVOR_CARRIER",
            "family": "quarks",
            "class": "research_open",
            "statement": (
                "the physical flavor carrier certificate and per-sector "
                "spread-fiber elimination; after the Clebsch lane the open "
                "quark objects are the Clebsch selection theorem, the "
                "third-generation register factor, and the charm/up "
                "selectors (the integer power-law family is removed by a "
                "negative frozen scan)"
            ),
        },
        {
            "id": "NEUTRINO_MECHANISM_SELECTION",
            "family": "neutrinos",
            "class": "research_open",
            "statement": (
                "Majorana or Dirac from the quotient content, then the "
                "absolute-scale mechanism and the neutral-lane family shape"
            ),
        },
        {
            "id": "WARD_ENDPOINT_AND_HADRONIC_TRANSPORT",
            "family": "cross-cutting (P root)",
            "class": "computationally_blocked",
            "statement": (
                "the source hadronic spectral transport (issue 425) and the "
                "scheme bridge (issue 545); the production route is the "
                "unquenched lattice backend with a vector correlator"
            ),
        },
        {
            "id": "OPERATIONAL_CLOCK_CHAIN",
            "family": "scale",
            "class": "mixed",
            "statement": (
                "R_alpha follows the P root; R_e_abs follows the charged "
                "completion; the cesium nuclear packet is computationally "
                "blocked; the Feshbach scalar evaluators are "
                "closable_by_standard_computation with an interval backend"
            ),
        },
        {
            "id": "BRST_POLE_KERNEL_PACKET",
            "family": "cross-cutting",
            "class": "research_open",
            "statement": (
                "BRST-complete two-point kernels with sheets, residues, and "
                "widths; every mass above is a tree/chart or threshold "
                "coordinate until this closes"
            ),
        },
        {
            "id": "NONPERTURBATIVE_HADRON_FACTOR",
            "family": "hadrons",
            "class": "computationally_blocked",
            "statement": (
                "the source-only m_hadron/Lambda factor; the conditional "
                "route through the published lattice ratio carries the "
                "spectrum in the interim"
            ),
        },
    ]


def build() -> dict[str, Any]:
    closed = closed_entries()
    open_items = open_entries()
    return {
        "artifact": "oph_theorem_gap_register",
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "status": "register_current",
        "promotion_allowed": False,
        "closed_this_program": closed,
        "open_register": open_items,
        "counts": {
            "closed": len(closed),
            "research_open": sum(
                1 for e in open_items if e["class"] == "research_open"
            ),
            "closable_by_standard_computation": sum(
                1
                for e in open_items
                if e["class"] == "closable_by_standard_computation"
            ),
            "computationally_blocked": sum(
                1 for e in open_items if e["class"] == "computationally_blocked"
            ),
            "awaiting_discriminating_test": sum(
                1
                for e in open_items
                if e["class"] == "awaiting_discriminating_test"
            ),
        },
    }


def render_markdown(register: dict[str, Any]) -> str:
    lines = [
        "# Theorem Gap Register",
        "",
        "Every named open object of the particle program with its class, "
        "plus the closures recorded by the current pass. Generated by "
        "`scripts/build_theorem_gap_register.py`.",
        "",
        "## Closed",
        "",
        "| Id | Family | Statement | Artifact |",
        "|---|---|---|---|",
    ]
    for entry in register["closed_this_program"]:
        lines.append(
            f"| {entry['id']} | {entry['family']} | {entry['statement']} "
            f"| `{entry['artifact']}` |"
        )
    lines += [
        "",
        "## Open",
        "",
        "| Id | Family | Class | Statement |",
        "|---|---|---|---|",
    ]
    for entry in register["open_register"]:
        lines.append(
            f"| {entry['id']} | {entry['family']} | {entry['class']} "
            f"| {entry['statement']} |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    parser.add_argument("--md-out", type=Path, default=DEFAULT_MD_OUT)
    args = parser.parse_args()
    register = build()
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(
        json.dumps(register, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.md_out.write_text(render_markdown(register), encoding="utf-8")
    print(
        json.dumps(
            {"status": register["status"], "counts": register["counts"]},
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
