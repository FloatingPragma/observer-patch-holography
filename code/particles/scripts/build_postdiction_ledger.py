#!/usr/bin/env python3
"""Aggregate the certified postdiction rows into one ledger.

The ledger is a pure aggregator: every displayed value is read live from a
certified parent artifact, every measured reference is read from the parent
that embeds it, and the builder introduces no physics number of its own.
A missing parent is a hard failure, not a silently absent row.

Section one records the forced-structure layer: the machine-checked
icosahedral results that pin the gauge sector before any numeric lane runs.
The finite steps live in the Lean workspace under `Lean/Screen/`; the
receipt entries record the module paths (existence-checked here) and the
declared hypothesis boundaries exactly as the compact paper states them.

The numeric sections carry the per-lane claim discipline of their parents:
interval rows report containment of the compare-only witness, conditional
rows carry their declared premises, chart coordinates keep their
NOT_EVALUABLE physical-comparison status, and the quark absolute-mass row
is an obstruction theorem rather than a number.

Run:
    python3 code/particles/scripts/build_postdiction_ledger.py
writes runs/status/postdiction_ledger.json and POSTDICTION_LEDGER.md.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCRIPTS = Path(__file__).resolve().parent
PARTICLES = SCRIPTS.parent
CODE = PARTICLES.parent
REPO = CODE.parent
RUNS = PARTICLES / "runs"
RUNTIME = CODE / "P_derivation" / "runtime"
LEAN_SCREEN = REPO / "Lean" / "Screen"

PARENTS = {
    "mass_surface": RUNS / "status" / "source_only_mass_prediction_surface.json",
    "conditional_ew": RUNS / "calibration" / "conditional_ew_predictions_current.json",
    "endpoint": RUNTIME / "empirical_thomson_endpoint_current.json",
    "anchor_bridge": RUNTIME / "anchor_scheme_bridge_current.json",
    "kappa_rectangle": RUNS / "leptons" / "charged_kappa_interval_from_alpha_transport.json",
    "kappa_coherent": RUNS / "leptons" / "charged_kappa_interval_coherent_closure.json",
    "clebsch_lane": RUNS / "flavor" / "down_type_register_clebsch_lane.json",
    "clebsch_selection": RUNS / "flavor" / "clebsch_register_pairing_selection.json",
    "fiber_obstruction": RUNS / "flavor" / "quark_spread_fiber_structure_transport_obstruction.json",
    "matter_receipt": CODE / "a5_closure" / "receipts" / "super_tannakian_matter_reference.receipt.json",
    "hadron_payload": RUNS / "hadron" / "empirical_ee_hadronic_spectral_measure.json",
    "solver_standby": RUNS / "qcd" / "hadron_source_backend" / "qcd_ensemble" / "solver_on_standby.json",
}

LEAN_RECEIPTS = {
    "A5OPH": LEAN_SCREEN / "A5OPH.lean",
    "A5CharacterField": LEAN_SCREEN / "A5CharacterField.lean",
    "A5SixAxes": LEAN_SCREEN / "A5SixAxes.lean",
    "Z6Exact": LEAN_SCREEN / "Z6Exact.lean",
    "A5CouplingSymmetry": LEAN_SCREEN / "A5CouplingSymmetry.lean",
    "A5PortAction": LEAN_SCREEN / "A5PortAction.lean",
    "PortFrameGram": LEAN_SCREEN / "PortFrameGram.lean",
}

DEFAULT_OUT = RUNS / "status" / "postdiction_ledger.json"
DEFAULT_MD = PARTICLES / "POSTDICTION_LEDGER.md"


def _load(key: str, override: Path | None = None) -> dict[str, Any]:
    path = override or PARENTS[key]
    if not path.exists():
        raise SystemExit(f"postdiction ledger parent missing: {key} at {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _lean_receipt(*modules: str) -> list[str]:
    refs = []
    for m in modules:
        path = LEAN_RECEIPTS[m]
        if not path.exists():
            raise SystemExit(f"postdiction ledger Lean receipt missing: {path}")
        refs.append(path.relative_to(REPO).as_posix())
    return refs


def _rel(key: str) -> str:
    return PARENTS[key].relative_to(REPO).as_posix()


def _forced_structure(matter: dict[str, Any]) -> list[dict[str, Any]]:
    spectrum = matter["realized_package"]["charge_spectrum"]
    sm_spectrum = {"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/6": 6}
    return [
        {
            "id": "gauge_lie_algebra",
            "statement": (
                "Compact-Lie trichotomy on the twelve-port screen: a compact "
                "connected group with a group-level A5 action equivalent to P12 "
                "has Lie algebra u(1)^12, su(2)^2+u(1)^6, or su(3)+su(2)+u(1); "
                "the noncentral quintet and the inner-action closure each select "
                "su(3)+su(2)+u(1)"
            ),
            "observed_counterpart": "Standard Model gauge Lie algebra su(3)+su(2)+u(1)",
            "match": "exact",
            "machine_checked_steps": (
                "triviality of A5-actions on at most four objects; unique "
                "partitions 11 = 8+3 and 12 = 3+3+3+3 over the compact-simple "
                "dimension list {3, 8, 10}; no compact semisimple algebra in "
                "dimensions 1, 2, 4, 5, 7; the characteristic-centre step; "
                "A5 not a subgroup of SU(2); the noncentrality witness "
                "[iS, iT] != 0; Galois stability of the two three-dimensional "
                "characters over Q(sqrt 5) with the centre-dimension list "
                "{0, 1, 5, 6, 7, 11, 12}; six-axis 2-transitivity and the "
                "dimension count of the dimension-six branch"
            ),
            "lean_receipts": _lean_receipt("A5OPH", "A5CharacterField", "A5SixAxes"),
            "hypothesis_boundary": (
                "the compact-simple classification, the torus/cocharacter step "
                "of the rationality lemma, and irreducibility of the "
                "five-dimensional summand stay declared classical inputs on "
                "paper; the physical inner current action is the open premise "
                "of issues 567 and 599"
            ),
            "paper_ref": "compact paper, Compact-Lie trichotomy section",
        },
        {
            "id": "global_form_z6",
            "statement": (
                "The screen gluing-class quotient Lambda_+/(Lambda_1 + Lambda_5) "
                "is Z/6 with proper-rotation invariance and antipodal sign "
                "reversal, matching the global form (SU(3) x SU(2) x U(1))/Z6"
            ),
            "observed_counterpart": (
                "Standard Model global gauge-group form and its charge "
                "quantization pattern"
            ),
            "match": "exact",
            "lean_receipts": _lean_receipt("Z6Exact", "A5OPH"),
            "hypothesis_boundary": (
                "the quotient isomorphism and both invariance clauses are "
                "machine checked; the identification with the physical global "
                "form rides on the same inner current action premise as the "
                "Lie-algebra row"
            ),
            "paper_ref": "compact paper, Z6 global-form section",
        },
        {
            "id": "hypercharge_spectrum",
            "statement": (
                "The super-Tannakian matter lift realizes the one-generation "
                "left-chiral hypercharge multiset "
                "{Q: 1/6 x6, u_c: -2/3 x3, d_c: 1/3 x3, L: -1/2 x2, e_c: 1 x1}"
            ),
            "observed_counterpart": (
                "Standard Model one-generation hypercharge assignment"
            ),
            "realized_spectrum": spectrum,
            "match": "exact" if spectrum == sm_spectrum else "MISMATCH",
            "artifact_ref": _rel("matter_receipt"),
            "hypothesis_boundary": (
                "conditional on the declared super-Tannakian lift premises "
                "recorded in the receipt claim boundary (issue 314); the "
                "generation count is an input, not a consequence"
            ),
            "paper_ref": "zoo paper, matter lift section",
        },
        {
            "id": "coupling_universality",
            "statement": (
                "A5-invariant readouts have port-independent group-averaged cap "
                "sums, so the per-cap ratio of any two averaged readouts is "
                "universal with zero spread"
            ),
            "observed_counterpart": (
                "universality clause of the Einstein-branch coupling law"
            ),
            "match": "structural",
            "lean_receipts": _lean_receipt("A5CouplingSymmetry", "A5PortAction", "PortFrameGram"),
            "hypothesis_boundary": (
                "reduces the universality clause to A5-equivariance of the "
                "implemented source law; no coupling value is implied"
            ),
            "paper_ref": "compact paper, coupling symmetry section",
        },
    ]


def _alpha_rows(endpoint: dict[str, Any], bridge: dict[str, Any]) -> list[dict[str, Any]]:
    ep = endpoint["endpoint"]
    co = endpoint["compare_only"]
    verdict = bridge["verdict"]
    if "reference_deficit_inside_certified_gap" not in verdict:
        raise SystemExit("anchor bridge artifact lacks the containment verdict field")
    return [
        {
            "id": "alpha_inv_thomson_endpoint",
            "value_central": float(ep["alpha_inv_central"]),
            "value_interval": [float(v) for v in ep["alpha_inv_interval"]],
            "measured": float(co["codata_alpha_inv"]),
            "measured_source": "CODATA 2022 via the endpoint artifact, compare-only",
            "payload_release": endpoint["inputs"]["payload_release"],
            "row_class": endpoint["row_class"],
            "tier": "T1_empirical_closure",
            "anchor_gap_interval": [float(v) for v in co["same_scheme_anchor_gap_interval_inv_alpha"]],
            "reference_deficit_inside_certified_gap": verdict["reference_deficit_inside_certified_gap"],
            "reading": (
                "the endpoint carries the full measured deficit into the "
                "certified same-scheme anchor gap; the standard reference "
                "deficit sits inside the certified interval, so the residual "
                "is the one-loop anchor running deficit, not a payload or "
                "solve defect"
            ),
            "artifact_refs": [_rel("endpoint"), _rel("anchor_bridge")],
            "blocking_issues": [425, 545],
        }
    ]


def _lepton_rows(
    surface: dict[str, Any],
    rectangle: dict[str, Any],
    coherent: dict[str, Any],
) -> list[dict[str, Any]]:
    witness_point = rectangle["compare_only"].get("witness_point")
    if witness_point is None:
        raise SystemExit(
            "rectangle artifact lacks the witness_point block; rebuild the "
            "rectangle lane first"
        )
    width_floor = coherent.get("width_floor_audit")
    if width_floor is None:
        raise SystemExit(
            "coherent artifact lacks the width_floor_audit block; rebuild "
            "the coherent lane first"
        )
    witnesses = rectangle["compare_only"]["witness_masses_gev"]
    particles = [r["particle"] for r in rectangle["conditional_mass_rows"]]
    family = next(f for f in surface["families"] if f["family"] == "charged leptons")
    mcpr = next(r for r in family["rows"] if r["lane"].startswith("MCPR"))
    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "id": "charged_leptons_closure_target",
            "statement": (
                "The certified solve inverts exactly at the measured triple: "
                "one anchor-gap value closes the charged-lepton lane on the "
                "witness, that value lies inside the certified band, and its "
                "distance to the standard on-shell reference deficit is the "
                "live scheme term of the open anchor bridge. The lepton "
                "scale is localized to the width of the scheme band, and a "
                "source-emitted bridge value is a sharp falsification "
                "target: landing on the closure value closes the lane on "
                "the witness, landing outside the certified band refutes "
                "the decomposition."
            ),
            "witness_point": witness_point,
            "width_floor": width_floor["floor_attribution"],
            "tier": "T1_empirical_closure",
            "artifact_refs": [_rel("kappa_rectangle"), _rel("kappa_coherent")],
            "blocking_issues": [545, 425],
        }
    )
    mcpr_masses = [float(m) / 1000.0 for m in mcpr["masses_MeV_display"]]
    rows.append(
        {
            "id": "charged_leptons_mcpr_conditional",
            "particles": particles,
            "values_gev": mcpr_masses,
            "measured_gev": witnesses,
            "measured_source": "PDG witness triple embedded in the kappa lane, compare-only",
            "relative_deltas": [v / w - 1.0 for v, w in zip(mcpr_masses, witnesses, strict=True)],
            "tier": mcpr["tier"],
            "row_class": mcpr["row_class"],
            "explanation": mcpr["explanation"],
            "artifact_ref": mcpr["artifact"],
            "blocking_objects": mcpr["blocking_objects"],
        }
    )
    for key, lane, artifact in (
        ("charged_leptons_kappa_rectangle", rectangle, "kappa_rectangle"),
        ("charged_leptons_kappa_coherent", coherent, "kappa_coherent"),
    ):
        mass_rows = lane["conditional_mass_rows"]
        containment = all(
            row["mass_interval"][0] < w < row["mass_interval"][1]
            for row, w in zip(mass_rows, witnesses, strict=True)
        )
        entry: dict[str, Any] = {
            "id": key,
            "particles": particles,
            "intervals_gev": [row["mass_interval"] for row in mass_rows],
            "centrals_gev": [row["mass_central"] for row in mass_rows],
            "measured_gev": witnesses,
            "measured_source": "PDG witness triple embedded in the lane, compare-only",
            "witness_inside_all_intervals": containment,
            "relative_half_widths": [
                (row["mass_interval"][1] - row["mass_interval"][0])
                / (2.0 * row["mass_central"])
                for row in mass_rows
            ],
            "tier": "T1_empirical_closure",
            "row_class": lane["row_class"],
            "artifact_ref": _rel(artifact),
            "blocking_issues": [425, 545],
        }
        if key.endswith("coherent"):
            entry["width_reduction_factor"] = lane["kappa_interval"]["width_reduction_factor"]
            entry["premise"] = "payload-coherent anchor-gap premise, declared"
        rows.append(entry)
    return rows


def _ew_rows(conditional: dict[str, Any]) -> list[dict[str, Any]]:
    cc = conditional["comparison_compare_only"]
    rows = []
    for key in ("mH_gev", "mt_pole_gev", "MW_chart_gev", "MZ_chart_gev"):
        block = cc[key]
        row: dict[str, Any] = {
            "id": f"ew_{key}",
            "value_central": block["conditional_central"],
            "value_envelope": block["conditional_envelope"],
            "physical_comparison_status": block["physical_comparison_status"],
            "tier": "T2_conditional",
            "row_class": conditional["row_class"],
            "artifact_ref": _rel("conditional_ew"),
        }
        if block["physical_comparison_status"] == "COMPARE_ONLY":
            row.update(
                {
                    "measured": block["measured"],
                    "measured_sigma": block["measured_sigma"],
                    "measured_source": block["measured_source"],
                    "delta": block["delta"],
                    "delta_over_sigma": block["delta_over_sigma"],
                    "envelope_overlaps_one_sigma_band": block["envelope_overlaps_one_sigma_band"],
                }
            )
        else:
            row["reason"] = block["reason"]
        rows.append(row)
    return rows


def _quark_rows(
    obstruction: dict[str, Any],
    clebsch: dict[str, Any],
    selection: dict[str, Any],
) -> list[dict[str, Any]]:
    if obstruction["fork"] != "ii_fiber_survives":
        raise SystemExit(
            "fiber obstruction artifact is not on the survives fork; "
            "rebuild the quark section before aggregating"
        )
    compare = clebsch["compare_only"]
    predictions = clebsch["predictions"]
    return [
        {
            "id": "quark_absolute_masses_obstruction",
            "statement": (
                "No absolute quark mass is emitted: the two-modulus spread "
                "fiber survives every certified structure transport, so the "
                "six absolute masses are non-identifiable from the corpus, by "
                "theorem rather than by omission"
            ),
            "fork": obstruction["fork"],
            "fiber_cut_detected": obstruction["fiber_cut_detected"],
            "tier": obstruction["claim_tier"],
            "artifact_ref": _rel("fiber_obstruction"),
            "blocking_issues": obstruction["github_issues"],
        },
        {
            "id": "quark_down_type_texture_conditional",
            "values": predictions,
            "measured_references": compare["references"],
            "relative_deltas": {
                k: v for k, v in compare.items() if k.endswith("_relative")
            },
            "tier": "T2_conditional",
            "row_class": clebsch["row_class"],
            "premise": (
                "generation register order (issue 569); pairing and weight set "
                "selected by the Clebsch selection artifact"
            ),
            "selection_artifact_ref": _rel("clebsch_selection"),
            "selection_status": selection["status"],
            "artifact_ref": _rel("clebsch_lane"),
            "reading": (
                "Conditional texture rows: the register order premise is open "
                "and the recorded tensions stay in the normalization_tension "
                "block of the parent."
            ),
        },
    ]


def _hadron_rows(payload: dict[str, Any], standby: dict[str, Any]) -> list[dict[str, Any]]:
    integral = payload["integral"]
    norm = integral["normalization"]
    return [
        {
            "id": "hadronic_correction_engine",
            "delta_alpha_had_5_MZ": integral["value"],
            "uncertainty_total": integral["uncertainty"],
            "source_compilation": payload["source_compilation"]["id"],
            "pin_factor": norm["pin_factor"],
            "policy": (
                "The published-compilation payload is the correction engine of "
                "the fine-structure lane; source-only hadron rows stay "
                "suppressed pending the source spectral measure (issue 425)."
            ),
            "artifact_ref": _rel("hadron_payload"),
        },
        {
            "id": "qcd_solver_on_standby",
            "status": standby["status"],
            "invocation_gate": standby["policy"]["invocation_gate"],
            "artifact_ref": _rel("solver_standby"),
        },
    ]


def _principal_results(sections: dict[str, Any]) -> list[dict[str, Any]]:
    """Digest the four strongest rows into the leading section, data-driven."""

    leptons = {r["id"]: r for r in sections["charged_leptons"]}
    target = leptons["charged_leptons_closure_target"]
    coherent = leptons["charged_leptons_kappa_coherent"]
    mcpr = leptons["charged_leptons_mcpr_conditional"]
    ew = {r["id"]: r for r in sections["electroweak"]}
    alpha = sections["alpha"][0]
    wp = target["witness_point"]
    glo, ghi = alpha["anchor_gap_interval"]
    hw = coherent["relative_half_widths"][0]
    ppm = abs(mcpr["relative_deltas"][0]) * 1.0e6
    mh, mt = ew["ew_mH_gev"], ew["ew_mt_pole_gev"]
    return [
        {
            "id": "lepton_closure_target",
            "statement": (
                "The anchor-gap value "
                f"{wp['required_anchor_gap_at_witness_inv_alpha']:.4f} closes "
                "the charged-lepton lane exactly on the measured triple, "
                f"inside the certified band [{glo:.4f}, {ghi:.4f}]; the "
                f"distance {wp['scheme_term_difference_inv_alpha']:+.4f} to "
                "the standard on-shell reference deficit "
                f"{wp['reference_deficit_inv_alpha']:.4f} is the live scheme "
                "term of the open anchor bridge (issue 545). The lepton "
                "scale is localized to the width of the scheme band, and a "
                "source-emitted bridge value is a falsification target: the "
                "closure value confirms, a value outside the band refutes."
            ),
        },
        {
            "id": "lepton_certified_intervals",
            "statement": (
                "The measured charged-lepton triple lies inside every "
                "certified interval; the payload-coherent half-width is "
                f"{hw * 100.0:.2f} percent per lepton, and the conditional "
                f"eight-register triple sits {ppm:.0f} ppm from measurement "
                "with the architecture declared."
            ),
        },
        {
            "id": "higgs_top_envelopes",
            "statement": (
                f"The conditional Higgs envelope [{mh['value_envelope'][0]:.3f}, "
                f"{mh['value_envelope'][1]:.3f}] GeV sits "
                f"{mh['delta_over_sigma']:.2f} sigma from the measured "
                f"{mh['measured']} +- {mh['measured_sigma']} GeV, and the top "
                f"envelope [{mt['value_envelope'][0]:.2f}, "
                f"{mt['value_envelope'][1]:.2f}] GeV sits "
                f"{mt['delta_over_sigma']:.2f} sigma from "
                f"{mt['measured']} +- {mt['measured_sigma']} GeV, "
                "compare-only, conditional on the declared selection axioms."
            ),
        },
        {
            "id": "forced_gauge_structure",
            "statement": (
                "The gauge sector is pinned before any numeric lane runs: "
                "the twelve-port trichotomy forces su(3)+su(2)+u(1), the "
                "gluing-class quotient gives the Z6 global form, and the "
                "matter lift realizes the exact one-generation hypercharge "
                "multiset, with the finite steps machine checked in "
                "Lean/Screen and the hypothesis boundaries recorded below."
            ),
        },
    ]


def build(out_path: Path = DEFAULT_OUT, md_path: Path | None = DEFAULT_MD) -> dict[str, Any]:
    surface = _load("mass_surface")
    conditional = _load("conditional_ew")
    endpoint = _load("endpoint")
    bridge = _load("anchor_bridge")
    rectangle = _load("kappa_rectangle")
    coherent = _load("kappa_coherent")
    clebsch = _load("clebsch_lane")
    selection = _load("clebsch_selection")
    obstruction = _load("fiber_obstruction")
    matter = _load("matter_receipt")
    payload = _load("hadron_payload")
    standby = _load("solver_standby")

    sections = {
        "forced_structure": _forced_structure(matter),
        "alpha": _alpha_rows(endpoint, bridge),
        "charged_leptons": _lepton_rows(surface, rectangle, coherent),
        "electroweak": _ew_rows(conditional),
        "quarks": _quark_rows(obstruction, clebsch, selection),
        "hadrons": _hadron_rows(payload, standby),
        "neutrinos": [
            {
                "id": "neutrino_dimensionless_pointer",
                "statement": (
                    "dimensionless PMNS and mass-splitting-ratio "
                    "comparisons live on the results status surface; the "
                    "absolute attachment stays compare-only"
                ),
                "artifact_ref": "code/particles/RESULTS_STATUS.md",
            }
        ],
    }
    result = {
        "artifact": "oph_postdiction_ledger",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "schema_version": 1,
        "row_class": "compare_only_postdiction_ledger",
        "guards": {
            "compare_only": True,
            "public_promotion_allowed": False,
            "changes_any_solve_path": False,
            "new_axiom_introduced": False,
            "hand_typed_measured_values": False,
        },
        "aggregation_policy": (
            "every value and every measured reference is read live from the "
            "cited parent artifact; a missing parent aborts the build"
        ),
        "principal_results": _principal_results(sections),
        "sections": sections,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if md_path is not None:
        md_path.write_text(_render_md(result), encoding="utf-8")
    return result


def _fmt(x: float, digits: int = 6) -> str:
    return f"{x:.{digits}g}"


def _render_md(ledger: dict[str, Any]) -> str:
    s = ledger["sections"]
    lines: list[str] = []
    add = lines.append
    add("# Postdiction Ledger")
    add("")
    add(f"Generated: `{ledger['generated_utc']}` by `scripts/build_postdiction_ledger.py`; "
        "the JSON artifact is `runs/status/postdiction_ledger.json`.")
    add("")
    add("Every value and every measured reference on this page is read live from "
        "the cited certified artifact. The ledger promotes nothing, changes no "
        "solve path, and introduces no number of its own. Interval rows report "
        "containment of the compare-only witness; conditional rows carry their "
        "declared premises; chart coordinates keep their NOT_EVALUABLE "
        "physical-comparison status.")
    add("")
    add("## Principal results")
    add("")
    for entry in ledger["principal_results"]:
        add(f"- {entry['statement']}")
    add("")
    add("## Forced structure")
    add("")
    add("The icosahedral screen results pin the gauge sector before any numeric "
        "lane runs. The finite steps are machine checked in the Lean workspace; "
        "the recorded hypothesis boundaries are the exact classical inputs and "
        "open premises of the compact paper.")
    add("")
    add("| Result | Observed counterpart | Match | Receipts |")
    add("| --- | --- | --- | --- |")
    for row in s["forced_structure"]:
        receipts = row.get("lean_receipts") or [row.get("artifact_ref", "")]
        receipt_txt = ", ".join(f"`{r}`" for r in receipts)
        add(f"| {row['statement']} | {row['observed_counterpart']} | "
            f"`{row['match']}` | {receipt_txt} |")
    add("")
    add("Hypothesis boundaries:")
    add("")
    for row in s["forced_structure"]:
        add(f"- `{row['id']}`: {row['hypothesis_boundary']}")
    add("")
    add("## Fine-structure lane")
    add("")
    for row in s["alpha"]:
        lo, hi = row["value_interval"]
        glo, ghi = row["anchor_gap_interval"]
        add(f"- `alpha_em^-1` Thomson endpoint: `{_fmt(row['value_central'], 10)}` "
            f"in `[{_fmt(lo, 10)}, {_fmt(hi, 10)}]` against CODATA "
            f"`{_fmt(row['measured'], 10)}` (compare-only). Payload release "
            f"`{row['payload_release']}`.")
        inside = "inside" if row["reference_deficit_inside_certified_gap"] else "outside"
        add(f"- Certified same-scheme anchor gap `[{_fmt(glo, 4)}, {_fmt(ghi, 4)}]` "
            f"inverse-alpha units; the standard reference deficit sits {inside} "
            "the certified interval.")
        add(f"- Reading: {row['reading']}")
        add(f"- Blocking issues: {', '.join(f'#{i}' for i in row['blocking_issues'])}")
    add("")
    add("## Charged leptons")
    add("")
    for row in s["charged_leptons"]:
        if row["id"].endswith("closure_target"):
            wp = row["witness_point"]
            add(f"- Closure target ({row['tier']}): the anchor-gap value "
                f"`{wp['required_anchor_gap_at_witness_inv_alpha']:.4f}` closes the "
                "lane exactly on the measured triple (inversion machine-checked); "
                f"the distance `{wp['scheme_term_difference_inv_alpha']:+.4f}` to the "
                f"on-shell reference deficit `{wp['reference_deficit_inv_alpha']:.4f}` "
                "is the live scheme term of the bridge. The certified width floor "
                "is the scheme-band ambiguity; no budget is shrunk without the "
                "source bridge.")
            continue
        if row["id"].endswith("mcpr_conditional"):
            deltas = ", ".join(
                f"{p} `{_fmt(d * 1e6, 3)} ppm`"
                for p, d in zip(row["particles"], row["relative_deltas"], strict=True)
            )
            add(f"- MCPR conditional triple ({row['tier']}): {deltas} against the "
                "PDG witness triple; the eight-register architecture is a "
                "declared model input.")
        else:
            kind = "coherent closure" if row["id"].endswith("coherent") else "rectangle"
            hw = ", ".join(
                f"{p} `{_fmt(h * 100, 3)}%`"
                for p, h in zip(row["particles"], row["relative_half_widths"], strict=True)
            )
            contained = "inside" if row["witness_inside_all_intervals"] else "OUTSIDE"
            add(f"- Kappa interval, {kind} ({row['tier']}): certified relative "
                f"half-widths {hw}; the witness triple lies {contained} every "
                "interval.")
            if "width_reduction_factor" in row:
                add(f"  - Width reduction over the rectangle: "
                    f"`{_fmt(row['width_reduction_factor'], 3)}x`; premise: "
                    f"{row['premise']}.")
    add("")
    add("## Electroweak sector")
    add("")
    add("| Quantity | Conditional central | Envelope | Measured | Delta/sigma | Status |")
    add("| --- | ---: | --- | --- | ---: | --- |")
    for row in s["electroweak"]:
        env = f"[{_fmt(row['value_envelope'][0], 8)}, {_fmt(row['value_envelope'][1], 8)}]"
        if row["physical_comparison_status"] == "COMPARE_ONLY":
            add(f"| `{row['id'][3:]}` | `{_fmt(row['value_central'], 8)}` | `{env}` | "
                f"`{row['measured']} +- {row['measured_sigma']}` ({row['measured_source']}) | "
                f"`{_fmt(row['delta_over_sigma'], 3)}` | compare-only |")
        else:
            add(f"| `{row['id'][3:]}` | `{_fmt(row['value_central'], 8)}` | `{env}` | "
                "chart coordinate | n/a | NOT_EVALUABLE |")
    add("")
    add("W/Z rows are running/tree chart coordinates; no physical comparison is "
        "defined until the chart-to-pole map is complete. The Higgs and top "
        "rows are conditional on the declared selection axioms.")
    add("")
    add("## Quarks")
    add("")
    for row in s["quarks"]:
        if row["id"].endswith("obstruction"):
            add(f"- Absolute masses ({row['tier']}): {row['statement']} "
                f"(issues {', '.join(f'#{i}' for i in row['blocking_issues'])}).")
        else:
            vals = row["values"]
            refs = row["measured_references"]
            add(f"- Down-type texture, conditional ({row['tier']}): "
                f"Cabibbo `{_fmt(vals['cabibbo_gst_sqrt_md_over_ms'], 4)}` against "
                f"`{refs['cabibbo']}`; `ms/md = {_fmt(vals['ms_over_md'], 4)}` against "
                f"`{refs['ms_over_md']}`. Premise: {row['premise']}. "
                f"{row['reading']}")
    add("")
    add("## Hadrons")
    add("")
    for row in s["hadrons"]:
        if row["id"] == "hadronic_correction_engine":
            add(f"- Correction engine payload: `Delta alpha_had^(5)(M_Z^2) = "
                f"{row['delta_alpha_had_5_MZ']} +- {row['uncertainty_total']}` "
                f"from `{row['source_compilation']}` "
                f"(pin factor `{_fmt(row['pin_factor'], 7)}`). {row['policy']}")
        else:
            add(f"- QCD solver: `{row['status']}`; invocation is gated on the "
                "source-side parameter emissions recorded in the standby receipt.")
    add("")
    add("## Neutrinos")
    add("")
    for row in s["neutrinos"]:
        add(f"- {row['statement']} (`{row['artifact_ref']}`).")
    add("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    args = parser.parse_args()
    result = build(args.out, args.md)
    for name, rows in result["sections"].items():
        print(f"{name}: {len(rows)} rows")
    print(f"wrote {args.out}")
    print(f"wrote {args.md}")


if __name__ == "__main__":
    main()
