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

The register also carries three conditional classical carrier-mode receipts.
They record the zero hard quadratic mass parameter on separately declared
Maxwell, perturbative Yang--Mills, and pure-Einstein branches.  The recovered
symmetry group does not supply those action, phase, background, or positivity
premises.  The receipts do not construct quantum particles or physical poles.
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
            "id": "KOIDE_CIRCULANT_IDENTITY",
            "family": "charged leptons",
            "class": "closed",
            "statement": (
                "inside the positive-eigenvalue chamber, the Hermitian C3 "
                "circulant has Q = 1/3 + (2/3)(|b|/a)^2 and hence Q = 2/3 "
                "iff |b|/a = 1/sqrt(2); the finite tracial-GNS result fixes "
                "that modulus only under its declared attachment premises, "
                "while the phase, two mass ratios, and physical attachment "
                "remain open"
            ),
            "artifact": "runs/leptons/koide_circulant_identity.json",
        },
        {
            "id": "W5_STABILISER_SPECTRUM_BOUND",
            "family": "charged leptons",
            "class": "closed",
            "statement": (
                "in W5 = Sym^2_0(R^3), C3 and C5 fixed loci are "
                "one-dimensional and force a double eigenvalue, whereas the "
                "C2 fixed locus is three-dimensional (projective dimension "
                "two) and admits a simple spectrum; stabilizer geometry "
                "therefore does not select the two charged-lepton mass "
                "ratios, which require a screen-derived potential"
            ),
            "artifact": "runs/flavor/w5_stabiliser_spectrum_bound.json",
        },
        {
            "id": "ICOSAHEDRAL_DIRECT_AXIS_CABIBBO_NO_GO",
            "family": "quark mixing",
            "class": "closed",
            "statement": (
                "the exhaustive exact spectrum of the 31 unoriented real "
                "three-dimensional icosahedral residual axes has minimum "
                "nonzero acute angle 20.905157 degrees, excluding only a "
                "direct identification of the Cabibbo angle with one such "
                "axis-pair angle; other A5 representations, symmetry "
                "breaking, and general overlap geometry are not excluded"
            ),
            "artifact": "runs/flavor/icosahedral_axis_angle_spectrum.json",
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
            "id": "DOWN_TYPE_REGISTER_CLEBSCH_ROUTE_REJECTED",
            "family": "quarks",
            "class": "closed",
            "statement": (
                "the channel-compatibility theorem and target-free F1/F2 "
                "enumeration leave the conditional unordered multiset "
                "{1/3, 1, 3}, but every one of its six generation "
                "assignments is retrospectively rejected by the conservative "
                "FLAG gate; the least-discrepant assignment is target-informed "
                "and does not select a physical generation order"
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
            "id": "MAXWELL_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT",
            "family": "conditional classical carrier modes",
            "class": "closed",
            "statement": (
                "on the separately declared unbroken, deconfined Maxwell "
                "branch with positive kinetic coefficient, the transverse "
                "classical modes have zero hard quadratic mass parameter; "
                "the OPH quantum-particle and positive-pole receipts are open"
            ),
            "artifact": "runs/status/carrier_mode_acceptance.json",
        },
        {
            "id": "YANG_MILLS_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT",
            "family": "conditional classical carrier modes",
            "class": "closed",
            "statement": (
                "on the separately declared perturbative Yang--Mills branch "
                "with positive kinetic coefficient, the free transverse "
                "classical modes have zero hard quadratic mass parameter; "
                "confinement blocks promotion to an asymptotic gluon particle"
            ),
            "artifact": "runs/status/carrier_mode_acceptance.json",
        },
        {
            "id": "EINSTEIN_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT",
            "family": "conditional classical carrier modes",
            "class": "closed",
            "statement": (
                "on the separately declared pure-Einstein branch about a "
                "suitable Ricci-flat background with positive kinetic "
                "coefficient, the transverse-traceless classical modes have "
                "zero hard quadratic mass parameter; no OPH graviton "
                "Hilbert space or positive-pole receipt is supplied"
            ),
            "artifact": "runs/status/carrier_mode_acceptance.json",
        },
        {
            "id": "SM_QFT_Q1_Q4_CONDITIONAL_IMPLICATIONS",
            "family": "W/Z and quantum-field landing",
            "class": "closed",
            "statement": (
                "the finite local G6 action implication, two exact finite "
                "QFT-Q2 criteria, formal QFT-Q3 BV/ST restoration, strict "
                "finite-order W/Z algebra, and separate QFT-Q4 OS/resonance "
                "implications are stated and checked without promoting their "
                "antecedents"
            ),
            "artifact": "calibration/q1_q4_wz_theorem_oracle/",
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
                "derive a screen-selected A5-invariant potential and prove "
                "that its global orbit selection emits the charged-family "
                "shape. Stabilizer symmetry alone cannot do this: C3/C5 "
                "force degeneracy and the C2 fixed locus retains exactly two "
                "projective shape parameters"
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
                "spread-fiber elimination. The complete declared "
                "register-Clebsch common-transport assignment family is "
                "retrospectively rejected inside its tested class. Direct "
                "Cabibbo identification with the canonical real-axis menu is "
                "excluded inside that separate finite class. "
                "channel compatibility and the conditional F1/F2 unordered "
                "multiset do not equate Yukawa coefficients or select a "
                "generation order. A source-derived mass and "
                "mixing mechanism, including charm/up selectors, is required"
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
            "id": "SM_QFT_NATIVE_PRODUCER_STACK",
            "family": "W/Z and quantum-field landing",
            "class": "research_open",
            "statement": (
                "construct the source-selected action/normalization, a "
                "QFT-Q2-E measure or QFT-Q2-H Hamiltonian, the "
                "regulator-specific QFT-Q3 restoration/matching/identity "
                "packet with dressed current amplitudes and numerical freeze, "
                "and the separate QFT-Q4 observable tower and resonance-sheet "
                "continuation packet"
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
