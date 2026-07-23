#!/usr/bin/env python3
"""Invert the declared Thomson transport for the charged-lepton scale kappa.

The kappa non-identifiability countermodel in ``charged_trace_lift_theorem.json``
rescales ``Y_e -> exp(kappa) Y_e`` while holding every listed antecedent fixed.
That family does not extend to the declared electromagnetic transport lane: the
lepton vacuum-polarization packet obeys ``d(packet)/d(kappa) = -2/pi != 0``, so
holding the emitted Thomson endpoint (equivalently the fixed point
``P = phi + sqrt(pi)/A_Th(P)``) fixed while rescaling the masses is
inconsistent.  On the empirical closure surface the transport therefore
identifies kappa up to certified interval width.

This builder solves, in the on-shell decomposition already recorded compare-only
by ``anchor_scheme_bridge_current.json``,

    a0 + g = alpha_inv_0 * (1 - Delta_lep(kappa) - Delta_had5 - Delta_top)

for kappa, with ``g`` ranging over the certified same-scheme anchor-gap
interval and ``Delta_had5`` over the empirical Ward-projected payload interval.
The mass ratios are consumed from the centered charged family functional, so
``Delta_lep`` is a strictly monotonic function of the single scale kappa and
the inversion is closed-form.

Row class ``target_shape_plus_empirical_transport``: this is never a source-only
theorem, introduces no axiom, and does not satisfy the production
``constructive_next_artifact`` of the source lane.  The solve consumes the
target-anchored charged shape, measured alpha_inv_0, and a higher-order
remainder calibrated at the measured charged triple.  The guard block declares
all three leaks.  The source-only no-go of the trace-lift theorem is unchanged.
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

EXACT_READOUT_JSON = RUNS / "leptons" / "lepton_current_family_exact_readout.json"
TRACE_LIFT_JSON = RUNS / "leptons" / "charged_trace_lift_theorem.json"
ENDPOINT_JSON = RUNTIME / "empirical_thomson_endpoint_current.json"
ANCHOR_BRIDGE_JSON = RUNTIME / "anchor_scheme_bridge_current.json"
HADRONIC_PROOF_PACKET_JSON = (
    RUNS / "hadron" / "ward_projected_spectral_measure_proof_packet.json"
)
DEFAULT_OUT = RUNS / "leptons" / "charged_kappa_interval_from_alpha_transport.json"

# PDG Z pole mass; identical to particles/hadron/ingest_empirical_ee_hadrons.py.
MZ_GEV = 91.1876
MASS_ORDER = ("electron", "muon", "tau")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def lepton_packet_asymptotic(m_e: float, ratios: tuple[float, float]) -> float:
    """Inverse-alpha lepton transport packet, one-loop asymptotic kernel.

    Matches ``paper_math.fermion_transport_kernel_asymptotic`` with Q^2 = 1 and
    multiplicity 1 per charged lepton, evaluated at q = MZ_GEV.
    """

    masses = (m_e, ratios[0] * m_e, ratios[1] * m_e)
    return sum(
        (1.0 / (3.0 * math.pi)) * (math.log(MZ_GEV**2 / m**2) - 5.0 / 3.0)
        for m in masses
    )


def invert_packet_for_m_e(packet: float, ratios: tuple[float, float]) -> float:
    """Closed-form inverse of ``lepton_packet_asymptotic`` in m_e."""

    log_ratio_sum = math.log(ratios[0]) + math.log(ratios[1])
    return math.exp(
        math.log(MZ_GEV) - (3.0 * math.pi * packet + 2.0 * log_ratio_sum + 5.0) / 6.0
    )


def required_delta_lep(
    alpha_inv_0: float, a0: float, gap: float, delta_had5: float, delta_top: float
) -> float:
    """Solve the on-shell decomposition for the leptonic Delta_alpha."""

    return 1.0 - delta_had5 - delta_top - (a0 + gap) / alpha_inv_0


def load_hadronic_contract_parent(
    packet_path: Path = HADRONIC_PROOF_PACKET_JSON,
) -> dict[str, Any]:
    """Gate and cite the issue-#317 spectral-measure proof packet, fail closed.

    The packet certifies the acceptance contract that the production
    ``constructive_next_artifact`` of this lane must satisfy.  It carries no
    numeric payload: ``Delta_had5`` consumed by the solve remains the
    declared-empirical companion of that contract.  If the physical source
    payload ever becomes available, this loader refuses to run so that the
    lane cannot silently keep the empirical payload past that point.
    """

    if not packet_path.exists():
        raise SystemExit(
            "fail closed: issue-317 proof packet missing at "
            f"{packet_path}; run its verifier before this lane"
        )
    packet = _load_json(packet_path)
    if packet.get("artifact") != "oph_ward_projected_spectral_measure_proof_packet":
        raise SystemExit(
            "fail closed: unexpected artifact in issue-317 proof packet: "
            f"{packet.get('artifact')!r}"
        )
    if packet.get("accepted") is not True:
        raise SystemExit(
            "fail closed: issue-317 proof packet is not accepted; the "
            "theorem-and-contract layer must certify before this lane cites it"
        )
    verdicts = packet.get("verdicts", {})
    if verdicts.get("contract_certified", {}).get("value") is not True:
        raise SystemExit(
            "fail closed: issue-317 contract_certified verdict is not true"
        )
    if verdicts.get("physical_source_payload_available", {}).get("value") is True:
        raise SystemExit(
            "fail closed: the physical source payload is reported available; "
            "re-derive this lane against the source-emitted spectral measure "
            "instead of the empirical payload before emitting any interval"
        )
    return {
        "artifact_ref": _artifact_ref(packet_path),
        "artifact": packet["artifact"],
        "issue": 317,
        "contract_certified": verdicts["contract_certified"]["value"],
        "physical_source_payload_available": verdicts[
            "physical_source_payload_available"
        ]["value"],
        "physical_promotion_allowed": verdicts.get(
            "physical_promotion_allowed", {}
        ).get("value"),
        "verifier_command": packet.get("verifier_command"),
        "reading": (
            "the production constructive_next_artifact of this lane carries a "
            "machine-certified acceptance contract (issue #317); the Delta_had5 "
            "payload consumed here remains the declared-empirical companion of "
            "that contract, and the certified reduction target of the interval "
            "width is unchanged (#425, #545)"
        ),
    }


def build(
    out_path: Path = DEFAULT_OUT,
    hadronic_packet_path: Path = HADRONIC_PROOF_PACKET_JSON,
) -> dict[str, Any]:
    readout = _load_json(EXACT_READOUT_JSON)
    endpoint = _load_json(ENDPOINT_JSON)
    bridge = _load_json(ANCHOR_BRIDGE_JSON)
    hadronic_contract_parent = load_hadronic_contract_parent(hadronic_packet_path)

    centered = [float(v) for v in readout["centered_log_shape_exact"]]
    ratios = (
        math.exp(centered[1] - centered[0]),
        math.exp(centered[2] - centered[0]),
    )
    witness = [float(v) for v in readout["predicted_singular_values_abs"]]

    a0 = float(endpoint["transport_split"]["a0_anchor_inv_alpha"])
    frozen_lepton_packet = float(
        endpoint["transport_split"]["lepton_transport_delta_inv_alpha"]
    )
    gap_interval = [
        float(v) for v in endpoint["compare_only"]["same_scheme_anchor_gap_interval_inv_alpha"]
    ]
    delta_had5 = float(endpoint["inputs"]["delta_alpha_had_5_MZ"])
    delta_had5_err = float(endpoint["inputs"]["delta_alpha_had_5_MZ_uncertainty"])

    reference = bridge["reference_decomposition_compare_only"]
    alpha_inv_0 = float(reference["alpha_inv_0"])
    delta_top = float(reference["Delta_top"])
    delta_lep_ref_3loop = float(reference["Delta_lep"])
    reference_gap = float(reference["gap_phys_minus_oph"])

    if delta_had5_err <= 0.0 or gap_interval[0] >= gap_interval[1]:
        raise SystemExit(
            "fail closed: certified hadronic uncertainty and a nondegenerate "
            "anchor-gap interval are required"
        )

    # Higher-order (>= 2 loop) leptonic remainder: the gap between the 3-loop
    # reference Delta_lep and the one-loop asymptotic kernel at the witness
    # triple.  Applied centrally with a [0, 2x] budget, fail-closed wide.
    packet_witness = lepton_packet_asymptotic(witness[0], ratios)
    delta_lep_1loop_witness = packet_witness / alpha_inv_0
    ho_remainder = delta_lep_ref_3loop - delta_lep_1loop_witness
    ho_budget = (0.0, 2.0 * ho_remainder) if ho_remainder > 0.0 else (2.0 * ho_remainder, 0.0)
    # Exact-vs-asymptotic one-loop kernel truncation, in packet units; bounded
    # by the recorded exact/asymptotic split of the frozen internal packet.
    kernel_truncation_packet = 5.0e-4

    def solve_kappa(gap: float, had: float, ho: float, kernel_slack: float) -> float:
        target = required_delta_lep(alpha_inv_0, a0, gap, had, delta_top) - ho
        packet = target * alpha_inv_0 + kernel_slack
        m_e = invert_packet_for_m_e(packet, ratios)
        return math.log(m_e / witness[0])

    # d(packet)/d(kappa) = -2/pi, so kappa is decreasing in the required packet:
    # the upper kappa endpoint takes every packet-lowering extreme and vice versa.
    kappa_hi = solve_kappa(
        gap_interval[1], delta_had5 + delta_had5_err, ho_budget[1], -kernel_truncation_packet
    )
    kappa_lo = solve_kappa(
        gap_interval[0], delta_had5 - delta_had5_err, ho_budget[0], kernel_truncation_packet
    )
    gap_mid = 0.5 * (gap_interval[0] + gap_interval[1])
    kappa_central = solve_kappa(gap_mid, delta_had5, ho_remainder, 0.0)
    kappa_reference_deficit = solve_kappa(reference_gap, delta_had5, ho_remainder, 0.0)
    # The anchor-gap value at which the solve lands exactly on the witness
    # triple (kappa = 0), at the payload central: the closure target a source
    # scheme bridge must emit.
    witness_required_gap = (
        alpha_inv_0 * (1.0 - delta_had5 - delta_top - ho_remainder)
        - packet_witness
        - a0
    )
    witness_closure_defect = solve_kappa(witness_required_gap, delta_had5, ho_remainder, 0.0)
    if abs(witness_closure_defect) > 1.0e-12:
        raise SystemExit(
            "fail closed: the witness-gap inversion does not close "
            f"(defect {witness_closure_defect})"
        )

    if not kappa_lo < kappa_central < kappa_hi:
        raise SystemExit("fail closed: kappa interval endpoints are not ordered")

    def mass_rows(k_lo: float, k_hi: float, k_c: float) -> list[dict[str, Any]]:
        rows = []
        factors = (1.0, ratios[0], ratios[1])
        for particle, factor in zip(MASS_ORDER, factors, strict=True):
            rows.append(
                {
                    "particle": particle,
                    "unit": "GeV",
                    "mass_interval": [
                        witness[0] * factor * math.exp(k_lo),
                        witness[0] * factor * math.exp(k_hi),
                    ],
                    "mass_central": witness[0] * factor * math.exp(k_c),
                    "status": "certified_empirical_closure_interval",
                    "formula": "m_i = exp(kappa) * R_i * m_e_witness, kappa from transport inversion",
                }
            )
        return rows

    # Interval width attribution, in kappa units (d ln m / d packet = -pi/2).
    dk_dpacket = math.pi / 2.0
    attribution = {
        "anchor_gap_half_width": 0.5 * (gap_interval[1] - gap_interval[0]) * dk_dpacket,
        "hadronic_payload_uncertainty": delta_had5_err * alpha_inv_0 * dk_dpacket,
        "higher_order_lepton_budget": 0.5 * abs(ho_budget[1] - ho_budget[0]) * alpha_inv_0 * dk_dpacket,
        "one_loop_kernel_truncation": kernel_truncation_packet * dk_dpacket,
        "reduction": (
            "every dominant term reduces to the source hadronic spectral measure "
            "(#425, hadron backend / lattice lane; its acceptance contract is "
            "certified by the issue-317 proof packet cited in "
            "inputs.hadronic_contract_parent) and the a0 scheme bridge (#545)"
        ),
    }

    # Stage-5 scale consistency: the frozen internal packet against the
    # witness-mass packet; both scale-invariant in their own units.
    stage5_kappa_offset = (packet_witness - frozen_lepton_packet) * dk_dpacket

    result = {
        "artifact": "oph_charged_kappa_interval_from_alpha_transport",
        "issue": 546,
        "generated_utc": _timestamp(),
        "row_class": "target_shape_plus_empirical_transport",
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
            "promotable_as_oph_source_theorem": False,
            "blind_normalization_prediction": False,
            "usable_for_public_final_values": False,
            "usable_as_diagnostic_route_finder": True,
            "satisfies_production_constructive_next_artifact": False,
        },
        "kappa_symmetry_breaking_lemma": {
            "statement": (
                "The one-parameter countermodel Y_e -> exp(kappa) Y_e of the charged "
                "trace-lift theorem does not extend to the declared electromagnetic "
                "transport: the lepton vacuum-polarization packet satisfies "
                "d(packet)/d(kappa) = -2/pi, so the emitted Thomson endpoint and the "
                "fixed point P = phi + sqrt(pi)/A_Th(P) move with kappa. The no-go's "
                "antecedent list omits this closure; source-only the system is the "
                "stiff curve (P(kappa), kappa), and on the empirical closure surface "
                "kappa is identified up to certified interval width."
            ),
            "packet_derivative_in_kappa": -2.0 / math.pi,
            "fixed_point_stiffness_dP_dkappa": (
                math.sqrt(math.pi)
                / float(endpoint["endpoint"]["alpha_inv_central"]) ** 2
            )
            * (2.0 / math.pi),
            "source_only_no_go_status": "unchanged",
        },
        "inversion_equation": {
            "decomposition": "a0 + g = alpha_inv_0 * (1 - Delta_lep(kappa) - Delta_had5 - Delta_top)",
            "packet": "Delta_lep = packet / alpha_inv_0, packet = (1/3pi) * sum_i (ln(MZ^2/m_i^2) - 5/3)",
            "closed_form": "ln m_e = ln MZ - (3pi*packet + 2*(ln R1 + ln R2) + 5)/6",
            "monotonicity": "packet strictly decreasing in kappa; unique solution",
            "mz_gev": MZ_GEV,
            "mz_provenance": "PDG, identical to particles/hadron/ingest_empirical_ee_hadrons.py",
        },
        "inputs": {
            "hadronic_contract_parent": hadronic_contract_parent,
            "a0_anchor_inv_alpha": a0,
            "a0_ref": _artifact_ref(ENDPOINT_JSON),
            "anchor_gap_interval": gap_interval,
            "anchor_gap_ref": _artifact_ref(ANCHOR_BRIDGE_JSON),
            "delta_alpha_had_5_MZ": [delta_had5 - delta_had5_err, delta_had5 + delta_had5_err],
            "delta_top": delta_top,
            "alpha_inv_0": alpha_inv_0,
            "alpha_inv_0_role": "declared empirical closure input inside this solve path (CODATA 2022)",
            "ratio_provenance": {
                "artifact_ref": _artifact_ref(EXACT_READOUT_JSON),
                "ratios_mu_over_e_tau_over_e": list(ratios),
                "scope": (
                    "centered charged family functional; target-anchored checksum per "
                    "charged_trace_lift_theorem.json, declared as such in this row class"
                ),
            },
            "higher_order_lepton_remainder": {
                "central": ho_remainder,
                "budget": list(ho_budget),
                "definition": "3-loop reference Delta_lep minus one-loop asymptotic kernel at witness",
            },
            "kernel_truncation_packet_budget": kernel_truncation_packet,
        },
        "kappa_interval": {
            "definition": "kappa = ln(m_e / m_e_witness)",
            "interval": [kappa_lo, kappa_hi],
            "central_gap_midpoint": kappa_central,
            "reference_deficit_point": {
                "gap": reference_gap,
                "kappa": kappa_reference_deficit,
                "reading": (
                    "with the payload pinned to the published KNT19 compilation the "
                    "physical on-shell anchor requirement sits inside the certified "
                    "interval; the residual kappa at that point is carried by the "
                    "payload width and the anchor-bridge source branch (#425, #545)"
                ),
            },
        },
        "conditional_mass_rows": mass_rows(kappa_lo, kappa_hi, kappa_central),
        "interval_width_attribution_kappa_units": attribution,
        "stage5_scale_consistency": {
            "frozen_internal_packet": frozen_lepton_packet,
            "witness_mass_packet": packet_witness,
            "implied_kappa_offset": stage5_kappa_offset,
            "reading": (
                "the Stage-5 continuation's scale-invariant content already matches the "
                "physical packet at the |kappa| ~ 0.006 level; the quarantined determinant "
                "prescription is consistency-tested, not free"
            ),
        },
        "compare_only": {
            "witness_masses_gev": witness,
            "witness_inside_certified_intervals": kappa_lo < 0.0 < kappa_hi,
            "witness_point": {
                "required_anchor_gap_at_witness_inv_alpha": witness_required_gap,
                "reference_deficit_inv_alpha": reference_gap,
                "scheme_term_difference_inv_alpha": witness_required_gap - reference_gap,
                "inside_certified_interval": (
                    gap_interval[0] <= witness_required_gap <= gap_interval[1]
                ),
                "reading": (
                    "the anchor-gap value at which the certified solve lands "
                    "exactly on the witness triple, at the payload central; a "
                    "source scheme bridge (#545) that emits this value closes "
                    "the lane on the witness, and the difference to the "
                    "standard on-shell reference deficit is the live scheme "
                    "term of the bridge"
                ),
            },
        },
        "claim_boundary": (
            "Absolute charged-lepton masses carry certified intervals on the empirical "
            "closure surface; the continuum kappa freedom is excluded there. No "
            "source-only absolute mass is emitted; the trace-lift no-go and its gate "
            "remain in force unchanged."
        ),
        "constructive_next_artifact": (
            "source_emitted_ward_projected_hadronic_spectral_measure_and_a0_scheme_bridge"
        ),
        "blind_prediction_route": {
            "stage_1_shape": {
                "required": "source_only_charged_mass_ratios_R_mu_e_and_R_tau_e",
                "current_status": "open_target_anchored_shape_only",
                "forbidden_inputs": ["m_e", "m_mu", "m_tau", "PDG", "CODATA"],
            },
            "stage_2_normalization": {
                "required": "source_emitted_Thomson_RG_packet_with_certified_hadronic_and_higher_order_remainders",
                "current_status": "open_empirical_transport_only",
                "mathematical_result": "unique_common_scale_because_d_packet_d_kappa_equals_minus_2_over_pi",
            },
            "stage_3_freeze_then_compare": {
                "required": "no_target_leak_DAG_and_frozen_mass_triple_before_PDG_comparison",
                "current_status": "not_ready",
            },
        },
        "proof_status": "certified_empirical_closure_interval_kappa_identified",
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
    rows = result["conditional_mass_rows"]
    print(f"kappa interval: [{interval[0]:+.4f}, {interval[1]:+.4f}]")
    for row in rows:
        lo, hi = row["mass_interval"]
        print(f"  {row['particle']:>8}: [{lo:.6e}, {hi:.6e}] GeV  central {row['mass_central']:.6e}")


if __name__ == "__main__":
    main()
