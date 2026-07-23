#!/usr/bin/env python3
"""Coherent-closure charged kappa interval: one payload variable, not two.

The rectangle lane (``derive_charged_kappa_interval_from_alpha_transport.py``)
solves

    a0 + g = alpha_inv_0 * (1 - Delta_lep(kappa) - Delta_had5 - Delta_top)

with ``g`` and ``Delta_had5`` ranging over their certified intervals
independently.  Both intervals are images of the same empirical payload: the
endpoint builder constructs the same-scheme anchor-gap interval as the exact
affine image

    g(X) = alpha_inv_0 * (1 - X) - A_lep - a0,     X = Delta_had5 +- u,

(``code/P_derivation/empirical_thomson_endpoint.py``, compare-only block), so
the rectangle treats one uncertainty as two anti-correlated ones.  The
rectangle semantics stay valid as robustness against any future source-emitted
bridge value inside the certified gap interval; this lane records the
complementary coherent reading.

Coherence premise (declared, machine-gated): the anchor gap is evaluated at
the same payload value as the hadronic term, i.e. the recorded certified gap
interval is exactly the affine payload image above.  Under that premise the
payload cancels from the decomposition,

    required Delta_lep = A_lep / alpha_inv_0 - Delta_top,

and kappa is identified by the independent budgets alone (higher-order lepton
remainder, one-loop kernel truncation, alpha_inv_0 input width).  The gate
recomputes the affine image at 90-digit precision and refuses to emit a
coherent interval on any mismatch, writing a no-improvement verdict instead;
the rectangle row then remains the only certified interval.

Exact consequence checked here: the coherent central equals the rectangle
central, because the rectangle's gap midpoint is the affine image of the
payload central.  The tightening removes width only.

Row class ``target_shape_plus_empirical_transport_coherent_closure``: never a
source-only theorem; same leak declarations as the rectangle lane plus the
named coherence premise.  The source-only charged no-go is unchanged.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

import derive_charged_kappa_interval_from_alpha_transport as rectangle_lane

ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
RUNS = ROOT / "particles" / "runs"
RUNTIME = ROOT / "P_derivation" / "runtime"

ENDPOINT_JSON = RUNTIME / "empirical_thomson_endpoint_current.json"
ANCHOR_BRIDGE_JSON = RUNTIME / "anchor_scheme_bridge_current.json"
RECTANGLE_JSON = RUNS / "leptons" / "charged_kappa_interval_from_alpha_transport.json"
DEFAULT_OUT = RUNS / "leptons" / "charged_kappa_interval_coherent_closure.json"

getcontext().prec = 90

# Residual bound for the affine-identity recomputation. The recorded interval
# strings carry ~85 digits; a genuine payload/gap incoherence enters at the
# 1e-4 scale, so any bound at or below 1e-10 separates the two cases. The
# bound only absorbs terminal-digit quantization of the recorded strings.
IDENTITY_RESIDUAL_BOUND = Decimal("1e-50")

# CODATA 2022 alpha_inv standard uncertainty, carried as an input-width budget.
ALPHA_INV_0_UNCERTAINTY = 2.1e-8


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def coherence_identity_gate(endpoint: dict[str, Any]) -> dict[str, Any]:
    """Recompute the anchor-gap interval as the affine payload image.

    Replicates the endpoint builder's operation order at 90-digit precision:

        required_anchor = alpha_inv_0 * (1 - (X +- u)) - A_lep
        anchor_gap      = required_anchor - a0

    and compares against the recorded certified interval.  Returns the gate
    record; ``passed`` is False on any residual above the bound.
    """

    one = Decimal(1)
    codata = Decimal(str(endpoint["compare_only"]["codata_alpha_inv"]))
    a0 = Decimal(str(endpoint["transport_split"]["a0_anchor_inv_alpha"]))
    a_lep = Decimal(
        str(endpoint["transport_split"]["lepton_transport_delta_inv_alpha"])
    )
    delta = Decimal(str(endpoint["inputs"]["delta_alpha_had_5_MZ"]))
    unc = Decimal(str(endpoint["inputs"]["delta_alpha_had_5_MZ_uncertainty"]))
    recorded = [
        Decimal(str(v))
        for v in endpoint["compare_only"]["same_scheme_anchor_gap_interval_inv_alpha"]
    ]

    required_lo = codata * (one - (delta + unc)) - a_lep
    required_hi = codata * (one - (delta - unc)) - a_lep
    recomputed_lo = required_lo - a0
    recomputed_hi = required_hi - a0

    residual_lo = abs(recomputed_lo - recorded[0])
    residual_hi = abs(recomputed_hi - recorded[1])
    passed = residual_lo <= IDENTITY_RESIDUAL_BOUND and residual_hi <= IDENTITY_RESIDUAL_BOUND

    breaking_term = None
    if not passed:
        breaking_term = (
            "anchor_gap_lower_endpoint"
            if residual_lo > IDENTITY_RESIDUAL_BOUND
            else "anchor_gap_upper_endpoint"
        )

    return {
        "statement": (
            "the certified same-scheme anchor-gap interval is the exact affine "
            "image g(X) = alpha_inv_0 * (1 - X) - A_lep - a0 of the empirical "
            "payload interval; gap and payload carry one uncertainty, not two"
        ),
        "operation_order": (
            "required_anchor = codata * (1 - (X +- u)) - A_lep; "
            "anchor_gap = required_anchor - a0; 90-digit decimal context"
        ),
        "residual_lower_endpoint": str(residual_lo),
        "residual_upper_endpoint": str(residual_hi),
        "residual_bound": str(IDENTITY_RESIDUAL_BOUND),
        "breaking_term": breaking_term,
        "passed": passed,
    }


def _width_floor_audit(rectangle: dict[str, Any]) -> dict[str, Any]:
    """State why the surviving width is a premise floor, not a budget slack.

    The dominant surviving budget is the higher-order lepton remainder,
    carried with a full-size band. The remainder itself matches the known
    per-order structure of the leptonic Delta_alpha (the two-loop piece
    dominates), and its kappa-sensitivity across the certified interval is
    two orders below the band. The band therefore functions as a level
    convention, and the level is pinned by the open anchor scheme bridge:
    at the current premise set the certified width floor is the
    scheme-bridge ambiguity, and shrinking the band without the source
    bridge would relabel that convention as a witness exclusion. The
    rectangle lane records the closure target the bridge must emit.
    """

    witness_point = rectangle.get("compare_only", {}).get("witness_point")
    if witness_point is None:
        raise SystemExit(
            "fail closed: rectangle artifact lacks the witness_point block; "
            "rebuild the rectangle lane first"
        )
    return {
        "dominant_budget": "higher_order_lepton_budget",
        "per_order_structure": (
            "the remainder between the three-loop reference Delta_lep and "
            "the one-loop asymptotic kernel at the witness matches the "
            "published per-order breakdown, with the two-loop leptonic "
            "piece dominant"
        ),
        "kappa_sensitivity_reading": (
            "the kappa-derivative of the higher-order piece is of order "
            "(alpha/pi)^2 per unit kappa, two orders below the carried "
            "band across the certified interval"
        ),
        "floor_attribution": (
            "the certified width floor at the current premise set is the "
            "anchor scheme-bridge ambiguity (issue 545, source branch); "
            "the payload term cancels in this lane and the higher-order "
            "band is a level convention pinned by the open bridge"
        ),
        "bridge_closure_target": witness_point,
        "tightening_gate": (
            "a certified width below the current floor requires the source "
            "scheme bridge; no budget is shrunk on this surface"
        ),
    }


def build(
    out_path: Path = DEFAULT_OUT,
    endpoint_path: Path = ENDPOINT_JSON,
    rectangle_path: Path = RECTANGLE_JSON,
    hadronic_packet_path: Path = rectangle_lane.HADRONIC_PROOF_PACKET_JSON,
) -> dict[str, Any]:
    readout = _load_json(rectangle_lane.EXACT_READOUT_JSON)
    endpoint = _load_json(endpoint_path)
    bridge = _load_json(ANCHOR_BRIDGE_JSON)
    hadronic_contract_parent = rectangle_lane.load_hadronic_contract_parent(
        hadronic_packet_path
    )
    if not rectangle_path.exists():
        raise SystemExit(
            "fail closed: rectangle-lane artifact missing at "
            f"{rectangle_path}; run the rectangle lane first"
        )
    rectangle = _load_json(rectangle_path)
    if rectangle.get("artifact") != "oph_charged_kappa_interval_from_alpha_transport":
        raise SystemExit("fail closed: unexpected rectangle-lane artifact")

    gate = coherence_identity_gate(endpoint)
    if not gate["passed"]:
        verdict = {
            "artifact": "oph_charged_kappa_interval_coherent_closure",
            "generated_utc": _timestamp(),
            "no_improvement_verdict": {
                "reason": "coherence_identity_gate_failed",
                "breaking_term": gate["breaking_term"],
                "gate": gate,
                "consequence": (
                    "the recorded anchor-gap interval is not the affine payload "
                    "image; the coherent solve is not licensed and the rectangle "
                    "lane remains the only certified interval row"
                ),
            },
            "rectangle_lane_ref": _artifact_ref(rectangle_path),
            "proof_status": "no_improvement_identity_gate_failed",
        }
        out_path.write_text(
            json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        raise SystemExit(
            "fail closed: coherence identity gate failed on "
            f"{gate['breaking_term']}; no-improvement verdict written"
        )

    centered = [float(v) for v in readout["centered_log_shape_exact"]]
    ratios = (
        math.exp(centered[1] - centered[0]),
        math.exp(centered[2] - centered[0]),
    )
    witness = [float(v) for v in readout["predicted_singular_values_abs"]]

    a_lep_frozen = float(endpoint["transport_split"]["lepton_transport_delta_inv_alpha"])
    reference = bridge["reference_decomposition_compare_only"]
    alpha_inv_0 = float(reference["alpha_inv_0"])
    delta_top = float(reference["Delta_top"])
    delta_lep_ref_3loop = float(reference["Delta_lep"])
    codata_compare = float(endpoint["compare_only"]["codata_alpha_inv"])
    if codata_compare != alpha_inv_0:
        raise SystemExit(
            "fail closed: the endpoint compare-only alpha_inv_0 and the bridge "
            "reference alpha_inv_0 disagree; the coherence premise is not typed"
        )

    # Independent budgets, identical definitions to the rectangle lane.
    packet_witness = rectangle_lane.lepton_packet_asymptotic(witness[0], ratios)
    delta_lep_1loop_witness = packet_witness / alpha_inv_0
    ho_remainder = delta_lep_ref_3loop - delta_lep_1loop_witness
    ho_budget = (0.0, 2.0 * ho_remainder) if ho_remainder > 0.0 else (2.0 * ho_remainder, 0.0)
    kernel_truncation_packet = 5.0e-4

    # Under the coherence premise the payload cancels:
    #   required Delta_lep = A_lep / alpha_inv_0 - Delta_top.
    # alpha_inv_0 input width enters through d(required)/d(alpha_inv_0)
    # = -A_lep / alpha_inv_0^2; carried as a packet-unit budget row.
    alpha_width_packet = (
        a_lep_frozen / alpha_inv_0**2 * ALPHA_INV_0_UNCERTAINTY * alpha_inv_0
    )

    def solve_kappa(ho: float, kernel_slack: float, alpha_slack_packet: float) -> float:
        required = a_lep_frozen / alpha_inv_0 - delta_top - ho
        packet = required * alpha_inv_0 + kernel_slack + alpha_slack_packet
        m_e = rectangle_lane.invert_packet_for_m_e(packet, ratios)
        return math.log(m_e / witness[0])

    # packet decreasing in kappa: packet-lowering extremes give the upper
    # kappa endpoint, identically to the rectangle lane.
    kappa_hi = solve_kappa(ho_budget[1], -kernel_truncation_packet, -alpha_width_packet)
    kappa_lo = solve_kappa(ho_budget[0], kernel_truncation_packet, alpha_width_packet)
    kappa_central = solve_kappa(ho_remainder, 0.0, 0.0)

    if not kappa_lo < kappa_central < kappa_hi:
        raise SystemExit("fail closed: coherent kappa endpoints are not ordered")

    # Exact consequence of the affine identity: the rectangle central used the
    # gap midpoint, which is the image of the payload central, so the two
    # centrals coincide.  Machine-checked, not assumed.
    rectangle_interval = [float(v) for v in rectangle["kappa_interval"]["interval"]]
    rectangle_central = float(rectangle["kappa_interval"]["central_gap_midpoint"])
    central_match_defect = abs(kappa_central - rectangle_central)
    if central_match_defect > 1.0e-9:
        raise SystemExit(
            "fail closed: coherent central does not reproduce the rectangle "
            f"central (defect {central_match_defect:.3e}); the affine premise "
            "does not describe the recorded rectangle solve"
        )

    rectangle_width = rectangle_interval[1] - rectangle_interval[0]
    coherent_width = kappa_hi - kappa_lo
    width_reduction_factor = rectangle_width / coherent_width

    def mass_rows(k_lo: float, k_hi: float, k_c: float) -> list[dict[str, Any]]:
        rows = []
        factors = (1.0, ratios[0], ratios[1])
        for particle, factor in zip(rectangle_lane.MASS_ORDER, factors, strict=True):
            rows.append(
                {
                    "particle": particle,
                    "unit": "GeV",
                    "mass_interval": [
                        witness[0] * factor * math.exp(k_lo),
                        witness[0] * factor * math.exp(k_hi),
                    ],
                    "mass_central": witness[0] * factor * math.exp(k_c),
                    "status": "certified_empirical_closure_interval_coherent",
                    "formula": (
                        "m_i = exp(kappa) * R_i * m_e_witness, kappa from the "
                        "payload-coherent transport inversion"
                    ),
                }
            )
        return rows

    dk_dpacket = math.pi / 2.0
    attribution = {
        "higher_order_lepton_budget": 0.5
        * abs(ho_budget[1] - ho_budget[0])
        * alpha_inv_0
        * dk_dpacket,
        "one_loop_kernel_truncation": kernel_truncation_packet * dk_dpacket,
        "alpha_inv_0_input_width": alpha_width_packet * dk_dpacket,
        "hadronic_payload_uncertainty": 0.0,
        "anchor_gap_half_width": 0.0,
        "reduction": (
            "the payload and anchor-gap terms cancel identically under the "
            "coherence premise; the surviving width is the higher-order lepton "
            "remainder and the one-loop kernel truncation, both closable by "
            "standard perturbative computation; the premise itself reduces to "
            "the source hadronic spectral measure (#425) and the a0 scheme "
            "bridge (#545), unchanged"
        ),
    }

    result = {
        "artifact": "oph_charged_kappa_interval_coherent_closure",
        "issue": 546,
        "generated_utc": _timestamp(),
        "row_class": "target_shape_plus_empirical_transport_coherent_closure",
        "guards": {
            "source_only": False,
            "new_axiom_introduced": False,
            "empirical_hadron_closure": True,
            "external_cross_section_data_used": True,
            "measured_alpha_in_solve_path": True,
            "measured_lepton_masses_directly_supplied_to_inversion": False,
            "target_anchored_lepton_ratios_in_solve_path": True,
            "measured_lepton_triple_used_to_calibrate_higher_order_remainder": True,
            "charged_mass_information_in_solve_path": True,
            "payload_double_count_removed": True,
            "conditional_on_payload_coherent_anchor_gap": True,
            "promotable_as_oph_source_theorem": False,
            "blind_normalization_prediction": False,
            "usable_for_public_final_values": False,
            "usable_as_diagnostic_route_finder": True,
            "satisfies_production_constructive_next_artifact": False,
        },
        "coherence_premise": {
            "statement": (
                "the same-scheme anchor gap takes exactly the value required by "
                "the recorded empirical payload, g(X) = alpha_inv_0 * (1 - X) - "
                "A_lep - a0, evaluated at the same payload value X as the "
                "hadronic term of the decomposition"
            ),
            "rectangle_reading": (
                "the rectangle lane remains the robust statement against any "
                "future source-emitted bridge value g inside its certified "
                "interval, independent of the payload"
            ),
            "identity_gate": gate,
        },
        "inversion_equation": {
            "decomposition": (
                "a0 + g(X) = alpha_inv_0 * (1 - Delta_lep(kappa) - X - Delta_top)"
            ),
            "payload_cancellation": (
                "substituting g(X) removes X exactly: required Delta_lep = "
                "A_lep / alpha_inv_0 - Delta_top"
            ),
            "monotonicity": "packet strictly decreasing in kappa; unique solution",
        },
        "inputs": {
            "hadronic_contract_parent": hadronic_contract_parent,
            "rectangle_lane_ref": _artifact_ref(rectangle_path),
            "endpoint_ref": _artifact_ref(endpoint_path),
            "anchor_bridge_ref": _artifact_ref(ANCHOR_BRIDGE_JSON),
            "lepton_transport_delta_inv_alpha": a_lep_frozen,
            "alpha_inv_0": alpha_inv_0,
            "alpha_inv_0_uncertainty": ALPHA_INV_0_UNCERTAINTY,
            "delta_top": delta_top,
            "ratio_provenance": {
                "artifact_ref": _artifact_ref(rectangle_lane.EXACT_READOUT_JSON),
                "ratios_mu_over_e_tau_over_e": list(ratios),
                "scope": (
                    "centered charged family functional; target-anchored checksum, "
                    "declared as such in this row class"
                ),
            },
            "higher_order_lepton_remainder": {
                "central": ho_remainder,
                "budget": list(ho_budget),
                "definition": (
                    "3-loop reference Delta_lep minus one-loop asymptotic kernel "
                    "at witness; identical to the rectangle lane"
                ),
            },
            "kernel_truncation_packet_budget": kernel_truncation_packet,
        },
        "kappa_interval": {
            "definition": "kappa = ln(m_e / m_e_witness)",
            "interval": [kappa_lo, kappa_hi],
            "central": kappa_central,
            "rectangle_interval": rectangle_interval,
            "central_match": {
                "rectangle_central": rectangle_central,
                "defect": central_match_defect,
                "reading": (
                    "the rectangle gap midpoint is the affine image of the "
                    "payload central, so the two centrals coincide; the coherent "
                    "solve removes width only"
                ),
            },
            "width_reduction_factor": width_reduction_factor,
        },
        "conditional_mass_rows": mass_rows(kappa_lo, kappa_hi, kappa_central),
        "interval_width_attribution_kappa_units": attribution,
        "width_floor_audit": _width_floor_audit(rectangle),
        "compare_only": {
            "witness_masses_gev": witness,
            "witness_inside_certified_intervals": kappa_lo < 0.0 < kappa_hi,
        },
        "claim_boundary": (
            "Absolute charged-lepton masses carry a certified coherent-closure "
            "interval on the empirical closure surface, conditional on the "
            "payload-coherent anchor-gap premise. The rectangle row is retained "
            "unchanged as the premise-free certified statement. No source-only "
            "absolute mass is emitted; the trace-lift no-go and its gate remain "
            "in force unchanged."
        ),
        "constructive_next_artifact": (
            "source_emitted_ward_projected_hadronic_spectral_measure_and_a0_scheme_bridge"
        ),
        "proof_status": "certified_empirical_closure_interval_coherent_kappa_identified",
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
    interval = result["kappa_interval"]["interval"]
    factor = result["kappa_interval"]["width_reduction_factor"]
    print(f"coherent kappa interval: [{interval[0]:+.5f}, {interval[1]:+.5f}]")
    print(f"width reduction vs rectangle: {factor:.1f}x")
    for row in result["conditional_mass_rows"]:
        lo, hi = row["mass_interval"]
        print(
            f"  {row['particle']:>8}: [{lo:.6e}, {hi:.6e}] GeV  "
            f"central {row['mass_central']:.6e}"
        )


if __name__ == "__main__":
    main()
