#!/usr/bin/env python3
"""Empirical Thomson endpoint on the oph_plus_empirical_hadron_closure surface.

Consumes the empirical e+e- hadronic spectral measure payload
(code/particles/runs/hadron/empirical_ee_hadronic_spectral_measure.json) and
evaluates the OPH fine-structure endpoint map with the frozen source anchor
and lepton transport packet:

    A_L = a0(P*) + Delta_lep(P*)
    denominator convention:  A_Th = A_L / (1 - Delta_alpha_had)
    additive convention:     A_Th = A_L + Delta_alpha_had * A_Th  (solved)
    P = phi + sqrt(pi) / A_Th,   alpha(0) = 1 / A_Th

The certified interval propagates the payload uncertainty through the
self-consistent insertion convention; the first-order truncated insertion
is reported as a diagnostic. The solve path contains no measured alpha value.
The CODATA comparison, the required hadronic shift, and the same-scheme
anchor gap are computed afterwards in an explicitly compare-only block, per
the forbidden-solver-inputs list of the Thomson endpoint contract.

The hadronic input is additionally required to exist as the declared
spectral-measure export (empirical_ward_projected_spectral_measure.json,
the declared-empirical companion of the Ward-projected export contract).
The endpoint refuses to evaluate if the export is absent, fails its own
requadrature consistency gate, or disagrees with the payload integral.
The solve path still reads the payload integral; the export supplies the
explicit spectral object and the spacelike OPH transport packet.

Row class: oph_plus_empirical_hadron_closure. The output is usable for
public final values on the empirical closure surface and is not promotable
as an OPH source theorem. The headline outputs are the certified endpoint
interval and the certified same-scheme anchor-gap interval; the anchor gap
is the number the source-side electroweak scheme bridge has to produce
(structure resolved in runtime/anchor_scheme_bridge_current.json: the gap
is the hadronic plus higher-order running deficit of the one-loop anchor).

Run:
    python3 code/P_derivation/empirical_thomson_endpoint.py
writes code/P_derivation/runtime/empirical_thomson_endpoint_current.json.
"""

from __future__ import annotations

import json
import pathlib
from datetime import datetime, timezone
from decimal import Decimal, getcontext

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
PAYLOAD_PATH = ROOT / "particles" / "runs" / "hadron" / "empirical_ee_hadronic_spectral_measure.json"
MEASURE_PATH = ROOT / "particles" / "runs" / "hadron" / "empirical_ward_projected_spectral_measure.json"
ANCHOR_BRIDGE_PATH = HERE / "runtime" / "anchor_scheme_bridge_current.json"
OUT_PATH = HERE / "runtime" / "empirical_thomson_endpoint_current.json"

getcontext().prec = 90

ONE = Decimal(1)
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")
PHI = (ONE + Decimal(5).sqrt()) / 2
SQRT_PI = PI.sqrt()

# Frozen OPH source-chain packets (identical to the fixed-point demo and the
# Thomson endpoint package; the anchor is a0(P*) = alpha_em^-1(m_Z^2; P*)).
SOURCE_ANCHOR_INV_ALPHA_MZ = Decimal(
    "128.30796547328624820996110874175671618724547618036535646005342169635117784168285644")
LEPTON_TRANSPORT_DELTA_INV_ALPHA = Decimal(
    "4.3093978664522040271317438975344894018487156605576773194711528089665680313257906466129")

# Compare-only reference, never read by the solve path.
CODATA_ALPHA_INV = Decimal("137.035999177")


def endpoint_denominator(a_l: Decimal, delta_had: Decimal) -> Decimal:
    return a_l / (ONE - delta_had)


def endpoint_additive(a_l: Decimal, delta_had: Decimal) -> Decimal:
    """A_Th = A_L + delta_had * A_Th, the additive inverse-alpha insertion."""
    return a_l / (ONE - delta_had)  # identical closed form; kept for clarity


def endpoint_additive_alpha0(a_l: Decimal, delta_had: Decimal) -> Decimal:
    """A_Th = A_L + delta_had / alpha(0) with alpha(0) = 1/A_Th: the strictly
    additive convention, A_Th = A_L + delta_had * A_Th, matches the
    denominator form; the leading-order additive packet uses A_L instead:
    A_Th = A_L + delta_had * A_L."""
    return a_l * (ONE + delta_had)


def pixel_from_endpoint(a_th: Decimal) -> Decimal:
    return PHI + SQRT_PI / a_th


def load_spectral_measure_export(payload: dict) -> dict:
    """Load and gate the declared-empirical Ward-projected spectral-measure
    export. The endpoint arithmetic requires the explicit spectral object,
    not only the integrated payload value."""
    if not MEASURE_PATH.exists():
        raise FileNotFoundError(
            f"missing spectral-measure export {MEASURE_PATH}; run "
            "particles/hadron/derive_empirical_ward_projected_spectral_measure.py")
    with open(MEASURE_PATH, encoding="utf-8") as f:
        measure = json.load(f)
    assert measure["artifact"] == "oph_empirical_ward_projected_hadronic_spectral_measure"
    assert measure["row_class"] == "oph_plus_empirical_hadron_closure"
    guards = measure["guards"]
    assert guards["promotable_as_oph_source_theorem"] is False
    assert guards["surrogate_hadron_artifact"] is False
    assert guards["satisfies_production_constructive_next_artifact"] is False
    consistency = measure["consistency"]
    assert consistency["within_tolerance"] is True, "spectral-measure requadrature gate failed"
    same_release = (measure["provenance"]["data_release"]["release_id"]
                    == payload["data_release"]["release_id"])
    assert same_release, "spectral-measure export and payload releases differ"
    payload_value = Decimal(str(payload["integral"]["value"]))
    export_value = Decimal(str(measure["transport_moments"]["timelike_on_shell_mz"]["value"]))
    tolerance = Decimal(str(consistency["tolerance"]))
    assert abs(export_value - payload_value) <= tolerance, (
        "spectral-measure export disagrees with the payload integral")
    return measure


def _anchor_bridge_block() -> dict:
    """Compare-only pointer to the anchor scheme-bridge analysis (#545). The
    empirical spectral measure supplies the route-A hadronic running; the
    source-only branch stays reduced to the hadron backend."""
    block = {
        "reference_artifact": "P_derivation/runtime/anchor_scheme_bridge_current.json",
        "route": "route_A_empirical_class_supplied_by_this_lane",
        "source_only_branch": "reduces to the OPH hadronic spectral measure, "
                              "blocked on the hadron backend (#425); #545 open "
                              "as that reduction",
    }
    if ANCHOR_BRIDGE_PATH.exists():
        with open(ANCHOR_BRIDGE_PATH, encoding="utf-8") as f:
            bridge = json.load(f)
        block["reference_status"] = bridge.get("verdict", {}).get("issue_545_status")
        block["gap_consistency"] = bridge.get("verdict", {}).get("gap_consistency")
    return block


def evaluate() -> dict:
    with open(PAYLOAD_PATH, encoding="utf-8") as f:
        payload = json.load(f)
    assert payload["artifact"] == "oph_empirical_ee_hadronic_spectral_measure"
    assert payload["guards"]["promotable_as_oph_source_theorem"] is False
    measure = load_spectral_measure_export(payload)

    delta_had = Decimal(str(payload["integral"]["value"]))
    delta_had_unc = Decimal(str(payload["integral"]["uncertainty"]))
    a_l = SOURCE_ANCHOR_INV_ALPHA_MZ + LEPTON_TRANSPORT_DELTA_INV_ALPHA

    # certified interval: payload uncertainty under the self-consistent
    # convention (the strictly additive inverse-alpha insertion solved to all
    # orders coincides with the denominator form; the truncated first-order
    # insertion is a diagnostic, not a convention, and is reported separately)
    a_th_lo = endpoint_denominator(a_l, delta_had - delta_had_unc)
    a_th_hi = endpoint_denominator(a_l, delta_had + delta_had_unc)
    a_th_central = endpoint_denominator(a_l, delta_had)
    truncation_diagnostic = a_th_central - endpoint_additive_alpha0(a_l, delta_had)

    p_lo = pixel_from_endpoint(a_th_hi)   # P decreases as A_Th grows
    p_hi = pixel_from_endpoint(a_th_lo)
    p_central = pixel_from_endpoint(a_th_central)

    # compare-only block (computed after the solve, never fed back)
    gap_central = a_th_central - CODATA_ALPHA_INV
    gap_lo = a_th_lo - CODATA_ALPHA_INV
    gap_hi = a_th_hi - CODATA_ALPHA_INV
    required_delta_had = ONE - a_l / CODATA_ALPHA_INV
    required_anchor_lo = CODATA_ALPHA_INV * (ONE - (delta_had + delta_had_unc)) \
        - LEPTON_TRANSPORT_DELTA_INV_ALPHA
    required_anchor_hi = CODATA_ALPHA_INV * (ONE - (delta_had - delta_had_unc)) \
        - LEPTON_TRANSPORT_DELTA_INV_ALPHA
    anchor_gap_lo = required_anchor_lo - SOURCE_ANCHOR_INV_ALPHA_MZ
    anchor_gap_hi = required_anchor_hi - SOURCE_ANCHOR_INV_ALPHA_MZ

    codata_inside = a_th_lo <= CODATA_ALPHA_INV <= a_th_hi

    return {
        "artifact": "oph_empirical_thomson_endpoint",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "row_class": "oph_plus_empirical_hadron_closure",
        "guards": {
            "source_only": False,
            "empirical_hadron_closure": True,
            "external_cross_section_data_used": True,
            "promotable_as_oph_source_theorem": False,
            "usable_for_public_final_values": True,
            "measured_alpha_in_solve_path": False,
        },
        "inputs": {
            "payload": PAYLOAD_PATH.relative_to(ROOT.parent).as_posix(),
            "payload_release": payload["data_release"]["release_id"],
            "delta_alpha_had_5_MZ": str(delta_had),
            "delta_alpha_had_5_MZ_uncertainty": str(delta_had_unc),
            "source_anchor_inv_alpha_MZ": str(SOURCE_ANCHOR_INV_ALPHA_MZ),
            "lepton_transport_delta_inv_alpha": str(LEPTON_TRANSPORT_DELTA_INV_ALPHA),
            "insertion_convention": "A_Th = A_L / (1 - delta), the "
                "self-consistent form; the solved additive insertion "
                "A_Th = A_L + delta A_Th coincides with it",
            "hadronic_spectral_measure_export": {
                "path": MEASURE_PATH.relative_to(ROOT.parent).as_posix(),
                "artifact": measure["artifact"],
                "profile_id": measure["profile_id"],
                "representation": measure["rho_had_or_measure"]["representation"],
                "positivity_status": measure["rho_had_or_measure"]["positivity_status"],
                "requadrature_abs_difference": str(
                    measure["consistency"]["abs_difference"]),
                "requadrature_tolerance": str(measure["consistency"]["tolerance"]),
            },
        },
        "transport_split": {
            "definition": "Delta_Th split required by the Ward-projected endpoint "
                          "lane: lepton vacuum-polarization transport, hadronic "
                          "spectral object, matching/scheme remainder, and "
                          "certified bounds, on the empirical closure surface",
            "a0_anchor_inv_alpha": str(SOURCE_ANCHOR_INV_ALPHA_MZ),
            "lepton_transport_delta_inv_alpha": str(LEPTON_TRANSPORT_DELTA_INV_ALPHA),
            "hadronic_spectral_object": {
                "artifact": measure["artifact"],
                "path": MEASURE_PATH.relative_to(ROOT.parent).as_posix(),
                "delta_alpha_had_5_MZ_timelike": str(
                    measure["transport_moments"]["timelike_on_shell_mz"]["value"]),
                "delta_alpha_had_5_MZ_spacelike": str(
                    measure["transport_moments"]["spacelike_mz"]["value"]),
                "inverse_alpha_packet_spacelike": str(
                    measure["transport_moments"]["inverse_alpha_packet_spacelike"]["value"]),
                "status": "declared_empirical_not_source_emitted",
            },
            "matching_scheme_remainder": "carried by the same-scheme anchor gap in "
                                         "the compare_only block; structure resolved "
                                         "in runtime/anchor_scheme_bridge_current.json",
            "certified_bounds": "payload uncertainty propagated through the "
                                "self-consistent insertion; see endpoint interval",
        },
        "endpoint": {
            "alpha_inv_central": str(a_th_central),
            "alpha_inv_interval": [str(a_th_lo), str(a_th_hi)],
            "alpha_central": str(ONE / a_th_central),
            "P_central": str(p_central),
            "P_interval": [str(p_lo), str(p_hi)],
            "interval_semantics": "payload uncertainty propagated through "
                                  "the self-consistent insertion convention",
            "first_order_truncation_diagnostic_inv_alpha": str(
                truncation_diagnostic),
        },
        "compare_only": {
            "codata_alpha_inv": str(CODATA_ALPHA_INV),
            "codata_source": "CODATA 2022 via NIST; compare-only, outside the solve path",
            "gap_central_inv_alpha": str(gap_central),
            "gap_interval_inv_alpha": [str(gap_lo), str(gap_hi)],
            "codata_inside_endpoint_interval": bool(codata_inside),
            "required_delta_alpha_had_for_codata": str(required_delta_had),
            "measured_delta_alpha_had_payload": str(delta_had),
            "same_scheme_anchor_gap_interval_inv_alpha": [
                str(anchor_gap_lo), str(anchor_gap_hi)],
            "anchor_gap_reading": "the certified shift the source-side "
                                  "electroweak scheme bridge must emit at the "
                                  "anchor a0(P) for the empirical-closure "
                                  "endpoint to meet the measured value",
            "anchor_bridge": _anchor_bridge_block(),
        },
        "verdict": {
            "status": "empirical_closure_endpoint_certified_gap_localized_at_anchor",
            "statement": "The empirical payload interval for the hadronic "
                         "transport excludes the value required to reach the "
                         "measured endpoint with the frozen source anchor; the "
                         "discrepancy is carried entirely by the same-scheme "
                         "anchor gap. The constructive next artifact is the "
                         "source-side electroweak scheme bridge for a0(P).",
            "constructive_next_artifact": "source_side_electroweak_scheme_bridge_for_a0",
        },
    }


def main() -> int:
    report = evaluate()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
        f.write("\n")
    ep = report["endpoint"]
    cmp_blk = report["compare_only"]
    print(f"alpha^-1 endpoint (empirical closure) = {Decimal(ep['alpha_inv_central']):.9f}")
    print(f"  interval = [{Decimal(ep['alpha_inv_interval'][0]):.6f}, "
          f"{Decimal(ep['alpha_inv_interval'][1]):.6f}]")
    print(f"P = {Decimal(ep['P_central']):.15f}  interval "
          f"[{Decimal(ep['P_interval'][0]):.12f}, {Decimal(ep['P_interval'][1]):.12f}]")
    print(f"compare-only gap to CODATA = {Decimal(cmp_blk['gap_central_inv_alpha']):.6f} "
          f"(CODATA inside interval: {cmp_blk['codata_inside_endpoint_interval']})")
    print(f"same-scheme anchor gap = "
          f"[{Decimal(cmp_blk['same_scheme_anchor_gap_interval_inv_alpha'][0]):.6f}, "
          f"{Decimal(cmp_blk['same_scheme_anchor_gap_interval_inv_alpha'][1]):.6f}]")
    print(f"wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
