#!/usr/bin/env python3
"""Build the per-family source-only mass prediction surface.

This status surface collects, for every particle family, the strongest
numerical mass coordinate the corpus emits, its claim tier, and the named
objects that block the next tier.  It reads generated lane artifacts and
writes one JSON artifact plus the repository document
``SOURCE_ONLY_MASS_PREDICTIONS.md``.

Tier ladder used by every row:

- T0 comparison-only: measured value or target-anchored postdiction; never an
  ancestor of a prediction.
- T1 empirical closure: OPH structure plus a declared empirical transport
  packet; certified intervals.
- T2 declared-model conditional: source coordinates propagated through a
  declared architecture whose selection is unproved.
- T3 source-only dimensionless: unique consequence of the source records plus
  proved selection theorems.
- T4 source-only absolute: T3 plus a closed operational clock ratio.

GeV and MeV displays that use the unclosed clock candidate inherit its
calibration-checksum status and are annotated per row.  This builder reads no
measured particle mass.
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

INPUTS = {
    "mcpr": RUNS / "leptons" / "charged_mcpr_completion_conditional.json",
    "kappa": RUNS / "leptons" / "charged_kappa_interval_from_alpha_transport.json",
    "ward": RUNS / "leptons" / "charged_ward_determinant_line.json",
    "wzh": RUNS / "calibration" / "wzh_residual_elimination_boundary.json",
    "ew": RUNS / "calibration" / "conditional_ew_predictions_current.json",
    "qme": RUNS / "calibration" / "qme_source_action_rigidity_mechanism.json",
    "clock_audit": RUNS / "clock" / "clock_source_energy_closure_audit.json",
    "feshbach": RUNS / "clock" / "cs133_feshbach_scalarization.json",
    "criticality": RUNS / "calibration" / "d11_criticality_boundary_scan.json",
    "criticality_compare": RUNS / "calibration" / "d11_criticality_comparison.json",
    "selection_audit": RUNS / "calibration" / "d11_boundary_scale_selection_audit.json",
    "lambda_qcd": RUNS / "qcd" / "lambda_qcd_source_transmutation.json",
    "nucleon": RUNS / "hadron" / "nucleon_mass_external_qcd_ratio.json",
    "clebsch": RUNS / "flavor" / "down_type_register_clebsch_lane.json",
    "up_scan": RUNS / "flavor" / "up_type_register_exponent_scan.json",
}

DEFAULT_JSON_OUT = RUNS / "status" / "source_only_mass_prediction_surface.json"
DEFAULT_MD_OUT = PARTICLES / "SOURCE_ONLY_MASS_PREDICTIONS.md"

DISPLAY_NOTE = (
    "GeV/MeV displays use the unclosed clock candidate; the clock audit "
    "classifies that decimal as a calibration checksum."
)


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _lepton_rows(mcpr: dict, kappa: dict, ward: dict) -> dict[str, Any]:
    pred = mcpr["conditional_prediction"]
    display = mcpr["optional_scale_display"]
    interval_rows = kappa["conditional_mass_rows"]
    return {
        "family": "charged leptons",
        "rows": [
            {
                "lane": "MCPR completed-response candidate",
                "tier": "T2",
                "explanation": (
                    "The eight-register architecture reproduces the charged "
                    "triple at 84 ppm with zero runtime charged reference, "
                    "and the architecture itself is a declared model input. "
                    "The source-only completion is the A5/W5 orbit program; "
                    "the homogeneous branch is proven shape-silent."
                ),
                "row_class": mcpr["row_class"],
                "masses_over_E_star": pred["masses_over_E_star"],
                "masses_MeV_display": display["masses_MeV"],
                "ratios": pred["ratios"],
                "artifact": "runs/leptons/charged_mcpr_completion_conditional.json",
                "blocking_objects": sorted(
                    mcpr["audit_boundary"]["completion_route"].keys()
                ),
            },
            {
                "lane": "alpha-transport kappa interval (empirical closure)",
                "tier": "T1",
                "explanation": (
                    "Certified intervals contain the physical triple; the "
                    "width reduces to the open hadronic transport and the "
                    "anchor bridge."
                ),
                "row_class": kappa["row_class"],
                "mass_rows_GeV": [
                    {
                        "particle": row["particle"],
                        "central": row["mass_central"],
                        "interval": row["mass_interval"],
                    }
                    for row in interval_rows
                ],
                "artifact": (
                    "runs/leptons/charged_kappa_interval_from_alpha_transport.json"
                ),
                "blocking_objects": [
                    "source hadronic spectral measure (issue 425)",
                    "anchor scheme bridge source branch (issue 545)",
                ],
            },
            {
                "lane": "Ward determinant-line inversion (affine-scale producer)",
                "tier": "theorem complete, physical parents open",
                "row_class": "fail_closed_source_packet_gate",
                "artifact": "runs/leptons/charged_ward_determinant_line.json",
                "blocking_objects": [
                    entry["id"] if isinstance(entry, dict) else str(entry)
                    for entry in ward["physical_source_packet"]["open_parents"]
                ],
            },
        ],
    }


def _boson_rows(wzh: dict, ew: dict) -> dict[str, Any]:
    two_law = wzh["strict_branch_two_law_boundary"]
    display = wzh["two_law_display"]
    bench = ew["rows_at_calibration_P"]["running_tree_value_law"]
    return {
        "family": "electroweak bosons",
        "rows": [
            {
                "lane": "strict source-audit branch, zero-selector law",
                "tier": "T2 (T3 once the discrete law selection is a theorem)",
                "explanation": (
                    "Tree/chart coordinates within 0.06 to 0.08 percent of "
                    "measurement with no target input; the residual is the "
                    "unmodeled pole packet plus the witness-level P root."
                ),
                "row_class": wzh["row_class"],
                "MW_over_E_star": two_law["zero_selector_law"]["MW_over_E_star"],
                "MZ_over_E_star": two_law["zero_selector_law"]["MZ_over_E_star"],
                "MW_GeV_display": display["zero_selector_law"]["MW_GeV"],
                "MZ_GeV_display": display["zero_selector_law"]["MZ_GeV"],
                "artifact": "runs/calibration/wzh_residual_elimination_boundary.json",
                "blocking_objects": [
                    "SOURCE_LAW_SELECTION_PRINCIPLE",
                    "BRST-complete pole kernels",
                    "operational clock",
                ],
            },
            {
                "lane": "strict source-audit branch, carrier value law",
                "tier": "T2 (T3 once the discrete law selection is a theorem)",
                "row_class": wzh["row_class"],
                "MW_over_E_star": two_law["nonzero_carrier_law"]["MW_over_E_star"],
                "MZ_over_E_star": two_law["nonzero_carrier_law"]["MZ_over_E_star"],
                "MW_GeV_display": display["nonzero_carrier_law"]["MW_GeV"],
                "MZ_GeV_display": display["nonzero_carrier_law"]["MZ_GeV"],
                "artifact": "runs/calibration/wzh_residual_elimination_boundary.json",
                "blocking_objects": [
                    "SOURCE_LAW_SELECTION_PRINCIPLE",
                    "C10 carrier prospective selection",
                    "BRST-complete pole kernels",
                    "operational clock",
                ],
            },
            {
                "lane": "declared benchmark surface (public endpoint branch continuation)",
                "tier": "T2; the public branch P consumes the measured Thomson endpoint",
                "explanation": (
                    "Target-free continuation on the calibration pixel; the "
                    "pixel itself sits on the endpoint branch, so this row "
                    "carries endpoint ancestry and validates the formula "
                    "stack rather than adding an independent prediction."
                ),
                "row_class": ew["row_class"],
                "MW_GeV": bench["MW_pole_gev"],
                "MZ_GeV": bench["MZ_pole_gev"],
                "artifact": "runs/calibration/conditional_ew_predictions_current.json",
                "blocking_objects": [
                    "same gates as the strict branch plus the declared surface"
                ],
            },
        ],
        "discrete_ambiguity_width_GeV": wzh["remaining_discrete_selection"][
            "certified_ambiguity_width_GeV"
        ],
    }


def _higgs_top_rows(
    wzh: dict, ew: dict, crit: dict, crit_compare: dict, audit: dict
) -> dict[str, Any]:
    ivp = wzh["literal_source_ivp"]
    bench = ew["rows_at_calibration_P"]["running_tree_value_law"]
    named_1l = crit["one_loop_named_boundaries"]
    named_2l = crit["two_loop_named_boundaries"]
    relation = crit_compare.get("mt_mh_relation_test", {}).get("two_loop", {})

    def pair(row: dict) -> dict[str, float]:
        return {
            "mt_pole_GeV_display": row["mt_pole_gev"],
            "mH_tree_GeV_display": row["mh_tree_gev"],
        }

    adopted_name = audit["adopted_working_conditional_branch"]["name"]
    adopted_two_loop = named_2l.get(adopted_name, {})

    return {
        "family": "Higgs boson and top quark",
        "rows": [
            {
                "lane": (
                    "adopted conditional branch: criticality at the "
                    "log-midpoint anchor E_star exp(-pi) P^(-1/6)"
                ),
                "tier": (
                    "T2, conditional on the frozen post-exposure "
                    "boundary-scale candidate"
                ),
                "explanation": (
                    "Two-loop values of the adopted branch; the candidate "
                    "was identified after the implied scale was computed "
                    "and is frozen prospectively against the three-loop "
                    "discriminating test. The flow-internal selection "
                    "route is closed by a no-go."
                ),
                "mt_pole_GeV_display": adopted_two_loop.get("mt_pole_gev"),
                "mH_tree_GeV_display": adopted_two_loop.get("mh_tree_gev"),
                "artifact": (
                    "runs/calibration/d11_boundary_scale_selection_audit.json"
                ),
                "blocking_objects": [
                    "CRITICALITY_BOUNDARY_SCALE_SELECTION_THEOREM "
                    "(candidates frozen; three-loop packet discriminates)",
                ],
            },
            {
                "lane": (
                    "double-criticality boundary family (lambda = 0 and "
                    "beta_lambda = 0 at one source scale; zero continuous "
                    "parameters)"
                ),
                "tier": "T2, conditional on the boundary-scale selection theorem",
                "explanation": (
                    "Both Yukawa-sector boundaries derive from the gauge "
                    "sector through the criticality law. The family over "
                    "named source scales brackets the measured pair; the "
                    "archived mu_U branch is the low edge and its deficit "
                    "decomposes into the boundary-scale choice plus loop "
                    "truncation, both named and quantified."
                ),
                "one_loop": {
                    name: pair(row) for name, row in named_1l.items()
                },
                "two_loop": {
                    name: pair(row)
                    for name, row in named_2l.items()
                    if isinstance(row, dict) and "mt_pole_gev" in row
                },
                "artifact": "runs/calibration/d11_criticality_boundary_scan.json",
                "blocking_objects": [
                    "CRITICALITY_BOUNDARY_SCALE_SELECTION_THEOREM",
                    "frozen three-loop RG and matching packet",
                    "BRST-complete pole kernels",
                ],
            },
            {
                "lane": "m_t to m_H criticality relation (fit-free curve, compare-only test)",
                "tier": "relation verified within the declared matching band",
                "explanation": (
                    "Along the fit-free criticality curve, the Higgs "
                    "coordinate at the measured top mass tests the relation "
                    "independently of the boundary-scale selection."
                ),
                "mH_on_curve_at_measured_mt_GeV": relation.get(
                    "mh_on_curve_at_measured_mt_gev"
                ),
                "mh_relative_residual": relation.get("mh_relative"),
                "implied_boundary_scale_GeV": relation.get(
                    "implied_boundary_scale_gev"
                ),
                "artifact": "runs/calibration/d11_criticality_comparison.json",
                "blocking_objects": [
                    "CRITICALITY_BOUNDARY_SCALE_SELECTION_THEOREM"
                ],
            },
            {
                "lane": "archived literal one-loop core at the mu_U boundary",
                "tier": "T2 (target-free tree coordinates; low edge of the family)",
                "explanation": (
                    "The 115.1 GeV Higgs coordinate is the archived branch "
                    "of the same family: criticality imposed at the gauge "
                    "unification scale with one-loop running. Its deficit "
                    "against measurement is the boundary-scale choice plus "
                    "truncation, quantified by the scan."
                ),
                "mH_tree_GeV_display": ivp["display_GeV_using_unclosed_clock"][
                    "mH_tree"
                ],
                "mt_QCD_converted_GeV_display": ivp[
                    "display_GeV_using_unclosed_clock"
                ]["mt_QCD_converted"],
                "artifact": "runs/calibration/wzh_residual_elimination_boundary.json",
                "blocking_objects": [
                    "CRITICALITY_BOUNDARY_SCALE_SELECTION_THEOREM",
                    "frozen three-loop RG and matching packet",
                ],
            },
            {
                "lane": "declared calibration surface",
                "tier": "target-anchored fit; not a prediction",
                "explanation": (
                    "These values are back-solved from the measured Higgs "
                    "and top through the synchronization-scale scan and the "
                    "live-exact adapters. They validate the formula stack "
                    "and carry zero predictive content."
                ),
                "mH_GeV": bench["mH_gev"],
                "mt_GeV": bench["mt_pole_gev"],
                "artifact": "runs/calibration/conditional_ew_predictions_current.json",
                "blocking_objects": [
                    "synchronization scale is target-ancestral",
                ],
            },
        ],
    }


def _static_families(status_rows: list[dict]) -> list[dict[str, Any]]:
    def statuses(group: str, particles: list[str]) -> dict[str, str]:
        return {
            row["particle"]: row["status"]
            for row in status_rows
            if row.get("particle") in particles
        }

    return [
        {
            "family": "quarks (u, d, s, c, b)",
            "rows": [
                {
                    "lane": "public theorem lane",
                    "tier": "T0 boundary (no forward sextet)",
                    "explanation": (
                        "The current axioms leave a continuous spread fiber "
                        "per sector: no Yukawa selector is emitted on the "
                        "strict lane. The top row lives in the Higgs family "
                        "through the criticality law; the down-type sector "
                        "carries a conditional Clebsch lane below."
                    ),
                    "statuses": statuses("Quarks", ["u", "d", "s", "c", "b"]),
                    "blocking_objects": [
                        "PHYSICAL_QF1_TO_QF9_FLAVOR_CARRIER_CERTIFICATE",
                        "CLEBSCH_REGISTER_SELECTION_THEOREM",
                        "charm/up selectors (integer power-law family removed)",
                    ],
                }
            ],
        },
        {
            "family": "neutrinos",
            "rows": [
                {
                    "lane": "public theorem lane",
                    "tier": "T0 (rejected target-informed candidate only)",
                    "explanation": (
                        "Mechanism selection (Majorana or Dirac) and the "
                        "absolute-scale mechanism are open structural "
                        "questions; no candidate survives the correlated "
                        "oscillation test."
                    ),
                    "statuses": statuses(
                        "Leptons", ["nu_e", "nu_mu", "nu_tau"]
                    ),
                    "blocking_objects": [
                        "Majorana versus Dirac mechanism selection",
                        "absolute scale mechanism",
                        "neutral-lane family shape (W5 orbit program)",
                    ],
                }
            ],
        },
        {
            "family": "massless carriers (photon, gluon)",
            "rows": [
                {
                    "lane": "structural masslessness",
                    "tier": "structural (receipt pending)",
                    "explanation": (
                        "Exact masslessness follows from unbroken gauge "
                        "invariance in the recovered core; the receipt is "
                        "bookkeeping, never a numerical program."
                    ),
                    "statuses": statuses("Bosons", ["gamma", "g (8 color states)"]),
                    "blocking_objects": [
                        "exact masslessness theorem receipts on the carrier lane"
                    ],
                }
            ],
        },
    ]


def _down_type_rows(clebsch: dict, up_scan: dict) -> dict[str, Any]:
    predictions = clebsch["predictions"]
    compare = clebsch["compare_only"]
    return {
        "family": "quarks, conditional lanes",
        "rows": [
            {
                "lane": "down-type sector from MCPR leptons via register Clebsch (1, 1/3, 3)",
                "tier": "T2, conditional on the Clebsch selection theorem",
                "explanation": (
                    "Ratios and the Gatto-Sartori-Tonin Cabibbo angle land "
                    "at the ten-percent scale; the absolute normalization "
                    "carries the named third-generation tension."
                ),
                "mb_mb_GeV": predictions["mb_mb_gev"],
                "ms_2GeV_GeV": predictions["ms_2gev_gev"],
                "md_2GeV_GeV": predictions["md_2gev_gev"],
                "cabibbo_gst": predictions["cabibbo_gst_sqrt_md_over_ms"],
                "cabibbo_relative_compare_only": compare["cabibbo_relative"],
                "artifact": "runs/flavor/down_type_register_clebsch_lane.json",
                "blocking_objects": clebsch["normalization_tension"]["open_objects"],
            },
            {
                "lane": "up-type integer exponent scan (frozen, negative)",
                "tier": "compare-only scan; law family removed",
                "explanation": (
                    "No frozen source-constant base gives integer exponents "
                    "for charm and up; the verdict removes the family "
                    "prospectively and charm/up stay research-open."
                ),
                "verdict": up_scan["status"],
                "artifact": "runs/flavor/up_type_register_exponent_scan.json",
                "blocking_objects": ["charm/up Yukawa selectors"],
            },
        ],
    }


def _hadron_rows(lambda_qcd: dict, nucleon: dict) -> dict[str, Any]:
    return {
        "family": "hadrons",
        "rows": [
            {
                "lane": "Lambda_QCD by dimensional transmutation of the source coupling",
                "tier": "T2 (conditional on declared threshold inputs)",
                "explanation": (
                    "The source strong coupling transmutes perturbatively "
                    "into the QCD scale with no hadronic input. This is the "
                    "perturbative half of every light-hadron mass."
                ),
                "lambda3_GeV_display": lambda_qcd["central"]["lambda3_gev"],
                "lambda3_interval_GeV": lambda_qcd["lambda3_interval_gev"],
                "lambda5_GeV_display": lambda_qcd["central"]["lambda5_gev"],
                "artifact": "runs/qcd/lambda_qcd_source_transmutation.json",
                "blocking_objects": [
                    "declared threshold inputs (bracket-swept)",
                    "strict no-target P root",
                ],
            },
            {
                "lane": "nucleon mass from source Lambda times external lattice ratio",
                "tier": "conditional: oph_plus_external_qcd_theory",
                "explanation": (
                    "The nonperturbative factor m_N/Lambda is a published "
                    "lattice-theory constant, declared as the external "
                    "condition. A source-only nucleon mass requires the "
                    "production hadron backend; this row is the realistic "
                    "conditional route."
                ),
                "m_nucleon_GeV_display": nucleon["prediction"][
                    "m_nucleon_gev_display"
                ],
                "m_nucleon_interval_GeV": nucleon["prediction"][
                    "m_nucleon_interval_gev_display"
                ],
                "compare_only_relative": nucleon["compare_only"][
                    "central_relative_difference"
                ],
                "artifact": "runs/hadron/nucleon_mass_external_qcd_ratio.json",
                "blocking_objects": [
                    "production hadron backend (source-only nonperturbative factor)",
                ],
            },
            {
                "lane": "other hadrons (pions, kaons, hyperons)",
                "tier": "T0",
                "explanation": (
                    "Chiral-sector masses additionally require the light "
                    "quark Yukawa selectors, which are open in the quark "
                    "family; the same external-ratio route extends once "
                    "those close."
                ),
                "blocking_objects": [
                    "light-quark Yukawa selectors (quark family)",
                    "nonperturbative ratios per channel",
                ],
            },
        ],
    }


def _scale_row(clock_audit: dict, feshbach: dict) -> dict[str, Any]:
    return {
        "family": "operational scale E_star",
        "rows": [
            {
                "lane": "clock-to-source-energy conversion",
                "tier": "conversion theorem proved; epsilon_Cs is a checksum",
                "explanation": (
                    "Absolute GeV output is equivalent to a source-only "
                    "operational clock ratio; four of the five clock "
                    "component producers are absent, so every GeV column is "
                    "a display through the checksum decimal."
                ),
                "candidate_E_star_GeV": clock_audit.get(
                    "clock_certificate_receipt", {}
                ).get("candidate_E_star_GeV"),
                "epsilon_role": "calibration_checksum_reproducing_displayed_G",
                "artifact": "runs/clock/clock_source_energy_closure_audit.json",
                "blocking_objects": clock_audit.get(
                    "required_next_source_objects", []
                ),
            },
            {
                "lane": "cesium Feshbach scalarization",
                "tier": "theorem checks pass on a synthetic fixture",
                "artifact": "runs/clock/cs133_feshbach_scalarization.json",
                "blocking_objects": feshbach.get("open_source_packets", []),
            },
        ],
    }


def build() -> dict[str, Any]:
    data = {key: _load(path) for key, path in INPUTS.items()}
    status_rows = _load(PARTICLES / "results_status.json")["rows"]

    families = [
        _lepton_rows(data["mcpr"], data["kappa"], data["ward"]),
        _boson_rows(data["wzh"], data["ew"]),
        _higgs_top_rows(
            data["wzh"],
            data["ew"],
            data["criticality"],
            data["criticality_compare"],
            data["selection_audit"],
        ),
        *_static_families(status_rows),
        _down_type_rows(data["clebsch"], data["up_scan"]),
        _hadron_rows(data["lambda_qcd"], data["nucleon"]),
        _scale_row(data["clock_audit"], data["feshbach"]),
    ]

    input_receipts = {
        key: {"path": str(path.relative_to(PARTICLES)), "exists": path.exists()}
        for key, path in INPUTS.items()
    }

    return {
        "artifact": "oph_source_only_mass_prediction_surface",
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "status": "per_family_claim_tier_surface_no_public_promotion",
        "promotion_allowed": False,
        "display_note": DISPLAY_NOTE,
        "tier_ladder": {
            "T0": "comparison-only; never an ancestor of a prediction",
            "T1": "empirical closure; certified intervals",
            "T2": "declared-model conditional",
            "T3": "source-only dimensionless",
            "T4": "source-only absolute (closed operational clock ratio)",
        },
        "selection_mechanism_receipt": (
            "runs/calibration/qme_source_action_rigidity_mechanism.json"
        ),
        "families": families,
        "inputs": input_receipts,
        "global_blocking_objects": {
            "X1_strict_no_target_P_alpha_root": (
                "source-audit branch witness; full endpoint proof requires the "
                "hadronic spectral transport, the scheme bridge, and the "
                "interval fixed-point proof"
            ),
            "X2_operational_clock_chain": data["clock_audit"].get(
                "required_next_source_objects", []
            ),
            "X3_source_law_selection": (
                "mechanism closed by the source-action rigidity theorem; the "
                "physical Standard-Model moment vector is open"
            ),
            "X4_provenance_freeze": (
                "canonical taint classes, transitive hashes, prospective freeze"
            ),
            "X5_interacting_pole_packet": (
                "BRST-complete two-point kernels, sheets, residues, widths"
            ),
        },
    }


def render_markdown(surface: dict[str, Any]) -> str:
    lines: list[str] = []
    add = lines.append
    add("# Source-Only Mass Prediction Surface")
    add("")
    add(
        "Per-family best numerical mass coordinates with claim tiers and the "
        "named objects blocking the next tier. Generated by "
        "`scripts/build_source_only_mass_prediction_surface.py`; the JSON "
        "artifact is `runs/status/source_only_mass_prediction_surface.json`. "
        "No row is a public promotion."
    )
    add("")
    add("| Tier | Meaning |")
    add("|---|---|")
    for tier, meaning in surface["tier_ladder"].items():
        add(f"| {tier} | {meaning} |")
    add("")
    add(surface["display_note"])
    add("")

    for family in surface["families"]:
        add(f"## {family['family']}")
        add("")
        for row in family["rows"]:
            add(f"- Lane: {row['lane']}")
            add(f"  - Tier: {row['tier']}")
            for key, value in row.items():
                if key in {"lane", "tier", "blocking_objects", "artifact"}:
                    continue
                if isinstance(value, (dict, list)):
                    add(f"  - {key}: `{json.dumps(value)}`")
                else:
                    add(f"  - {key}: {value}")
            if "artifact" in row:
                add(f"  - Artifact: `{row['artifact']}`")
            blockers = row.get("blocking_objects", [])
            if blockers:
                add(f"  - Blocking objects: {', '.join(str(b) for b in blockers)}")
            add("")
        if "discrete_ambiguity_width_GeV" in family:
            width = family["discrete_ambiguity_width_GeV"]
            add(
                f"Discrete two-law ambiguity width: "
                f"{1000 * width['MW']:.1f} MeV on MW, "
                f"{1000 * width['MZ']:.1f} MeV on MZ."
            )
            add("")

    add("## Global blocking objects")
    add("")
    for key, value in surface["global_blocking_objects"].items():
        if isinstance(value, list):
            add(f"- {key}: {', '.join(str(v) for v in value)}")
        else:
            add(f"- {key}: {value}")
    add("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    parser.add_argument("--md-out", type=Path, default=DEFAULT_MD_OUT)
    args = parser.parse_args()
    surface = build()
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(
        json.dumps(surface, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.md_out.write_text(render_markdown(surface), encoding="utf-8")
    print(
        json.dumps(
            {
                "status": surface["status"],
                "families": [f["family"] for f in surface["families"]],
                "json": str(args.json_out),
                "markdown": str(args.md_out),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
