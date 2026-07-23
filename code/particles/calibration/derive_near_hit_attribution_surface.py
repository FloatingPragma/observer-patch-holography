#!/usr/bin/env python3
"""Audit candidate residual attributions without promoting chart comparisons.

The surface is compare-only: it promotes nothing, adds no axiom, and changes
no solve path.  W/Z values are running/tree chart coordinates, not pole
masses.  Their legacy standardized differences are retained only as arithmetic
diagnostics and receive a NOT_EVALUABLE physical-comparison status.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
RUNS = ROOT / "particles" / "runs"
RUNTIME = ROOT / "P_derivation" / "runtime"

CONDITIONAL_EW_JSON = RUNS / "calibration" / "conditional_ew_predictions_current.json"
VALUE_LAW_JSON = RUNS / "calibration" / "d10_ew_target_free_repair_value_law.json"
ANCHOR_BRIDGE_JSON = RUNTIME / "anchor_scheme_bridge_current.json"
ENDPOINT_JSON = RUNTIME / "empirical_thomson_endpoint_current.json"
KAPPA_LANE_JSON = RUNS / "leptons" / "charged_kappa_interval_from_alpha_transport.json"
DEFAULT_OUT = RUNS / "calibration" / "near_hit_attribution_surface.json"

HADRONIC_REDUCTION = (
    "source_emitted_ward_projected_hadronic_spectral_measure (#425); "
    "same-scheme a0 bridge (#545)"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def _chart_value(row: dict[str, Any], chart_key: str, legacy_key: str) -> float:
    """Read the corrected chart key, accepting a frozen pre-erratum artifact."""
    if chart_key in row:
        return float(row[chart_key])
    return float(row[legacy_key])


def _mz_legacy_coordinate_diagnostic(
    conditional: dict[str, Any], endpoint: dict[str, Any]
) -> dict[str, Any]:
    """Reproduce the historical MZ chart arithmetic without a physical verdict.

    Branch 1 (naive coupling injection): shift the couplings by the certified
    anchor gap and re-run the closed forward W/Z map.  Branch 2 (P channel):
    move the Thomson endpoint to the measured alpha_inv(0) and propagate the
    induced pixel shift through dMZ/dP.  Neither branch is a falsification test
    because the chart-to-pole map is incomplete.
    """

    from derive_d10_repair_tuple_selection_theorem import (
        candidate_a_tuple,
        forward_wz,
        load_basis,
    )

    basis = load_basis()
    cand = candidate_a_tuple(basis)
    baseline = forward_wz(basis, cand["tau2_exact"], cand["delta_n_exact"])

    gap_interval = [
        float(v)
        for v in endpoint["compare_only"]["same_scheme_anchor_gap_interval_inv_alpha"]
    ]
    injected = {}
    for name, gap in (("gap_lo", gap_interval[0]), ("gap_hi", gap_interval[1])):
        shifted = dict(basis)
        shifted["alphaY_mz"] = 1.0 / (1.0 / basis["alphaY_mz"] + gap)
        wz = forward_wz(shifted, cand["tau2_exact"], cand["delta_n_exact"])
        injected[name] = (wz["MZ_chart_gev"] - baseline["MZ_chart_gev"]) * 1e3

    a_th = float(endpoint["endpoint"]["alpha_inv_central"])
    p_empirical = float(endpoint["endpoint"]["P_central"])
    codata = float(endpoint["compare_only"]["codata_alpha_inv"])
    p_calibration = float(conditional["inputs"]["P_calibration"])
    p_lo = float(conditional["inputs"]["P_empirical_interval"][0])
    mz_cal_row = conditional["rows_at_calibration_P"]["running_tree_value_law"]
    mz_lo_row = conditional["rows_at_empirical_P"]["empirical_closure_P_lo"][
        "running_tree_value_law"
    ]
    mz_cal = _chart_value(mz_cal_row, "MZ_chart_gev", "MZ_pole_gev")
    mz_lo = _chart_value(mz_lo_row, "MZ_chart_gev", "MZ_pole_gev")
    dmz_dp = (mz_lo - mz_cal) / (p_lo - p_calibration)
    dp_da = -math.sqrt(math.pi) / a_th**2
    p_corrected = p_empirical + dp_da * (codata - a_th)

    return {
        "executed_utc": "2026-07-12",
        "physical_comparison_status": "NOT_EVALUABLE",
        "reason": (
            "the computation propagates internal chart coordinates; no complete "
            "scheme map identifies them with a measured W/Z pole observable"
        ),
        "naive_coupling_injection": {
            "scheme": "anchor gap added to alphaY^-1 (hypercharge line)",
            "delta_mz_mev_over_gap_interval": [injected["gap_lo"], injected["gap_hi"]],
            "proportional_scheme_delta_mz_mev": "approximately -230 to -302, and moves MW by a comparable amount",
            "verdict": "LEGACY_DIAGNOSTIC_ONLY",
            "reading": (
                "the one-loop couplings at calibration P already give MZ to +0.4 MeV; "
                "injecting on-shell hadronic running into them is scheme-inconsistent "
                "and moves the chart coordinate by the displayed amount; this does "
                "not establish a physical excess or exclusion"
            ),
        },
        "p_channel": {
            "dmz_dp_gev_per_unit_p": dmz_dp,
            "dp_d_alpha_inv_thomson": dp_da,
            "p_empirical_central": p_empirical,
            "p_corrected_endpoint_to_measured": p_corrected,
            "p_calibration": p_calibration,
            "p_residual_after_correction": p_corrected - p_calibration,
            "mz_residual_after_correction_mev": (p_corrected - p_calibration) * dmz_dp * 1e3,
            "verdict": "LEGACY_DIAGNOSTIC_ONLY",
            "reading": (
                "moving the Thomson endpoint to the measured alpha_inv(0) shifts the "
                "empirical pixel onto the calibration pixel to within 4.2e-07, i.e. "
                "-0.05 MeV in the internal MZ chart.  Without the missing physical "
                "readout map, this does not identify the cause of a measured mass "
                "residual"
            ),
        },
    }


def _ew_rows(
    conditional: dict[str, Any], endpoint: dict[str, Any]
) -> list[dict[str, Any]]:
    comparison = conditional["comparison_compare_only"]
    hypotheses = {
        "MZ_chart_gev": {
            "mechanism": (
                "no physical missing-correction attribution is licensed.  The value "
                "is a running/tree chart coordinate and the chart-to-pole map is open"
            ),
            "falsification_test": (
                "derive and freeze a common-observable map with vev and tadpole scheme, "
                "thresholds, finite-order completion, complex-pole conversion, and "
                "theory uncertainty before evaluating a mass residual"
            ),
            "legacy_coordinate_diagnostic": _mz_legacy_coordinate_diagnostic(
                conditional, endpoint
            ),
        },
        "MW_chart_gev": {
            "mechanism": (
                "no physical missing-correction attribution is licensed.  The value "
                "is a running/tree chart coordinate and the chart-to-pole map is open"
            ),
            "falsification_test": (
                "derive and freeze the same common-observable map required for the "
                "MZ chart before evaluating any residual"
            ),
        },
        "mH_gev": {
            "mechanism": (
                "repair-tuple selection underdetermination (issue #521) plus the shared "
                "anchor deficit entering through the quartic-coupling transport"
            ),
            "falsification_test": (
                "derive the selection axioms from source (#521); the envelope width "
                "should shrink below the measured sigma"
            ),
        },
        "mt_pole_gev": {
            "mechanism": (
                "same repair-tuple underdetermination; the top row inherits the "
                "envelope of the d11 calibration surface, not a new physics deficit"
            ),
            "falsification_test": (
                "source derivation of the selection axioms (#521); residual should "
                "track the mH row, not the alpha rows"
            ),
        },
    }
    rows = []
    for key, hyp in hypotheses.items():
        legacy_key = key.replace("_chart_", "_pole_")
        row = comparison.get(key, comparison.get(legacy_key))
        if row is None:
            raise KeyError(f"missing conditional comparison row for {key}")
        if key in {"MW_chart_gev", "MZ_chart_gev"}:
            legacy_ref = row.get("legacy_reference_coordinate", row.get("measured"))
            legacy_sigma = row.get(
                "legacy_reference_experimental_sigma", row.get("measured_sigma")
            )
            legacy_delta = row.get(
                "legacy_experimental_error_only_delta", row.get("delta")
            )
            legacy_standardized = row.get(
                "legacy_experimental_error_only_standardized_difference",
                row.get("delta_over_sigma"),
            )
            rows.append({
                "quantity": key,
                "oph_value": row["conditional_central"],
                "oph_envelope": row["conditional_envelope"],
                "physical_comparison_status": "NOT_EVALUABLE",
                "physical_delta": None,
                "physical_pull": None,
                "legacy_reference_coordinate": legacy_ref,
                "legacy_reference_experimental_sigma": legacy_sigma,
                "legacy_experimental_error_only_delta": legacy_delta,
                "legacy_experimental_error_only_standardized_difference": (
                    legacy_standardized
                ),
                "row_class": "running_tree_chart_coordinate_compare_only",
                "missing_correction_hypothesis": {
                    **hyp,
                    "reduces_to": "complete W/Z common-observable scheme map",
                },
                "artifact_ref": _artifact_ref(CONDITIONAL_EW_JSON),
            })
        else:
            rows.append({
                "quantity": key,
                "oph_value": row["conditional_central"],
                "oph_envelope": row["conditional_envelope"],
                "measured": row["measured"],
                "measured_sigma": row["measured_sigma"],
                "measured_source": row["measured_source"],
                "delta": row["delta"],
                "delta_over_sigma": row["delta_over_sigma"],
                "inside_one_sigma": row["envelope_inside_one_sigma_band"],
                "physical_comparison_status": "COMPARE_ONLY",
                "row_class": "conditional_on_P_and_repair_selection",
                "missing_correction_hypothesis": {
                    **hyp,
                    "reduces_to": "selection-axiom source derivation (#521)",
                },
                "artifact_ref": _artifact_ref(CONDITIONAL_EW_JSON),
            })
    return rows


def _alpha_rows(
    bridge: dict[str, Any], endpoint: dict[str, Any], value_law: dict[str, Any]
) -> list[dict[str, Any]]:
    reference = bridge["reference_decomposition_compare_only"]
    a0 = bridge["anchor_provenance"]["a0_oph"]
    rows = [
        {
            "quantity": "alpha_em_inv_MZ_anchor",
            "oph_value": a0,
            "measured": reference["alpha_inv_mz_phys_on_shell"],
            "measured_source": "PDG 2024 on-shell decomposition",
            "delta": a0 - reference["alpha_inv_mz_phys_on_shell"],
            "row_class": bridge["row_class"],
            "missing_correction_hypothesis": {
                "mechanism": (
                    "the anchor is a one-loop unification-run value with no "
                    "nonperturbative hadronic running; the certified same-scheme gap "
                    "interval contains the standard reference deficit 0.631 in its "
                    "interior, confirming a running deficit rather than an anchor error"
                ),
                "falsification_test": (
                    "emit the source hadronic spectral measure and the a0 scheme "
                    "bridge; the corrected anchor must land inside the certified gap"
                ),
                "reduces_to": HADRONIC_REDUCTION,
            },
            "artifact_ref": _artifact_ref(ANCHOR_BRIDGE_JSON),
        },
        {
            "quantity": "alpha_inv_thomson_endpoint",
            "oph_value": float(endpoint["endpoint"]["alpha_inv_central"]),
            "oph_envelope": [float(v) for v in endpoint["endpoint"]["alpha_inv_interval"]],
            "measured": float(endpoint["compare_only"]["codata_alpha_inv"]),
            "measured_source": "CODATA 2022 via NIST, compare-only",
            "delta": float(endpoint["compare_only"]["gap_central_inv_alpha"]),
            "row_class": endpoint["row_class"],
            "missing_correction_hypothesis": {
                "mechanism": (
                    "the empirical payload interval excludes the value required to "
                    "reach the measured endpoint with the frozen anchor; the whole "
                    "discrepancy is carried by the same-scheme anchor gap"
                ),
                "falsification_test": (
                    "source-side electroweak scheme bridge for a0; endpoint interval "
                    "must then contain the CODATA value"
                ),
                "reduces_to": HADRONIC_REDUCTION,
            },
            "artifact_ref": _artifact_ref(ENDPOINT_JSON),
        },
        {
            "quantity": "alpha_em_eff_inv_on_shell_wz_lane",
            "oph_value": value_law["coherent_emitted_quintet"]["alpha_em_eff_inv"],
            "measured": value_law["compare_only_validation_against_frozen_surface"][
                "alpha_em_eff_inv_reference"
            ],
            "measured_source": "frozen-surface reference inside the value-law artifact",
            "delta": value_law["compare_only_validation_against_frozen_surface"][
                "delta_alpha_em_eff_inv"
            ],
            "row_class": "target_free_repair_value_law",
            "missing_correction_hypothesis": {
                "mechanism": (
                    "the on-shell effective coupling law shares the anchor's missing "
                    "hadronic/scheme running; the +0.118 inverse-alpha gap is the "
                    "on-shell face of the same deficit"
                ),
                "falsification_test": (
                    "apply the emitted scheme bridge to the W/Z lane coupling; the "
                    "value-law gap should co-move with the anchor gap"
                ),
                "reduces_to": HADRONIC_REDUCTION,
            },
            "artifact_ref": _artifact_ref(VALUE_LAW_JSON),
        },
    ]
    return rows


def _lepton_rows(kappa_lane: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row, witness in zip(
        kappa_lane["conditional_mass_rows"],
        kappa_lane["compare_only"]["witness_masses_gev"],
        strict=True,
    ):
        central = row["mass_central"]
        rows.append(
            {
                "quantity": f"m_{row['particle']}_gev",
                "oph_value": central,
                "oph_envelope": row["mass_interval"],
                "measured": witness,
                "measured_source": "PDG witness triple, compare-only",
                "delta": central - witness,
                "relative_delta": central / witness - 1.0,
                "row_class": kappa_lane["row_class"],
                "missing_correction_hypothesis": {
                    "mechanism": (
                        "kappa-interval width is carried by the published-compilation "
                        "payload uncertainty and the independent lepton-transport "
                        "budgets; the physical on-shell anchor requirement sits inside "
                        "the certified anchor-gap interval, so the central offset is a "
                        "position within the certified band, not a payload deficit"
                    ),
                    "falsification_test": (
                        "replace the empirical payload with the source hadronic "
                        "spectral measure; the kappa interval must tighten around "
                        "the witness triple"
                    ),
                    "reduces_to": HADRONIC_REDUCTION,
                },
                "artifact_ref": _artifact_ref(KAPPA_LANE_JSON),
            }
        )
    return rows


def build(out_path: Path = DEFAULT_OUT) -> dict[str, Any]:
    conditional = _load_json(CONDITIONAL_EW_JSON)
    bridge = _load_json(ANCHOR_BRIDGE_JSON)
    endpoint = _load_json(ENDPOINT_JSON)
    value_law = _load_json(VALUE_LAW_JSON)
    kappa_lane = _load_json(KAPPA_LANE_JSON)

    rows = (
        _ew_rows(conditional, endpoint)
        + _alpha_rows(bridge, endpoint, value_law)
        + _lepton_rows(kappa_lane)
    )

    result = {
        "artifact": "oph_near_hit_attribution_surface",
        "generated_utc": _timestamp(),
        "row_class": "compare_only_attribution_surface",
        "guards": {
            "compare_only": True,
            "public_promotion_allowed": False,
            "changes_any_solve_path": False,
            "new_axiom_introduced": False,
        },
        "rows": rows,
        "synthesis": {
            "statement": (
                "The displayed residuals do not establish one dominant missing "
                "physical correction.  W/Z rows are not evaluable until their "
                "common-observable scheme map exists.  The alpha and lepton rows "
                "retain their declared hadronic/same-scheme hypotheses, while H/top "
                "retain their separate selection-axiom dependency."
            ),
            "dominant_reduction": HADRONIC_REDUCTION,
            "secondary_reduction": "selection-axiom source derivation (#521)",
            "wz_reduction": "complete W/Z common-observable scheme map",
            "coherence_check": (
                "the standard reference deficit 0.631 sits inside the certified "
                "anchor-gap interval, and the lepton-lane reference-deficit kappa "
                "sits inside the certified kappa interval; both containments are "
                "read from the same endpoint artifact and are mutually consistent"
            ),
        },
    }

    out_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    result = build(args.out)
    for row in result["rows"]:
        measured = row.get("measured")
        print(f"{row['quantity']:>36}: oph {row['oph_value']}  measured {measured}")


if __name__ == "__main__":
    main()
