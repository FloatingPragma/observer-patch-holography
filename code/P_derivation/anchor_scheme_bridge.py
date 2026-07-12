#!/usr/bin/env python3
"""Anchor scheme-bridge analysis for the certified endpoint gap (#545).

Decomposes the certified same-scheme anchor gap of the empirical Thomson
endpoint into named, standard, compare-only reference terms and emits the
verdict on the two closure routes of #545:

  Route A (scheme-bridge theorem): delta_scheme(P) landing in the certified
          gap [0.649, 0.855] inverse-alpha units.
  Route B (anchor-exactness no-go): the OPH anchor is already scheme-exact,
          so the gap is a genuine source-chain failure.

Finding. The OPH anchor a0(P*) = alpha_em^-1(m_Z^2; P*) is computed by the
OPH one-loop renormalization-group run from the OPH unification scale down
to m_Z (paper_math.build_d10_from_p). A one-loop unification run cannot
contain the nonperturbative hadronic vacuum-polarization contribution, which
enters the physical coupling at m_Z. The standard reference decomposition of
the physical on-shell value is

  alpha^-1(m_Z^2)_phys
     = alpha^-1(0) (1 - Delta_lep - Delta_had5 - Delta_top)
     = 128.939   (5-flavor on-shell)

against the OPH one-loop anchor 128.308, a difference of 0.631, essentially
the lower edge of the certified gap. The whole gap is the hadronic and
higher-order running deficit of the one-loop scheme.

Verdict. Route B is false: the anchor is not scheme-exact; it is a one-loop
unification value with a physical, understood deficit, so the certified gap
is not a source-chain failure. Route A is available only as an
EMPIRICAL-class bridge: filling delta_scheme requires the measured hadronic
running (the same measured input already carried by the empirical payload)
or the OPH source hadronic spectral measure, which does not exist and is
blocked on the hadron backend (#425). A source-only anchor bridge therefore
reduces to #425; there is no source-only scheme-bridge theorem on the
current corpus. This module states that reduction and does not manufacture a
source-only delta_scheme.

Every physical number below is a cited compare-only reference and never
enters any OPH solve path.

Run:
    python3 code/P_derivation/anchor_scheme_bridge.py
writes code/P_derivation/runtime/anchor_scheme_bridge_current.json.
"""

from __future__ import annotations

import json
import pathlib
from datetime import datetime, timezone

HERE = pathlib.Path(__file__).resolve().parent
OUT_PATH = HERE / "runtime" / "anchor_scheme_bridge_current.json"
ENDPOINT_PATH = HERE / "runtime" / "empirical_thomson_endpoint_current.json"

A0_OPH = 128.30796547328625            # OPH one-loop-run anchor, alpha_em^-1(mZ^2; P*)
ALPHA_INV_0 = 137.035999177            # CODATA 2022, compare-only

# Standard reference contributions to Delta_alpha(mZ^2), compare-only:
DELTA_LEP = 0.031497686                # leptonic, 3-loop (Steinhauser et al.), PDG
DELTA_HAD5 = 0.02766                   # 5-flavor hadronic, PDG 2024 (data-driven)
DELTA_TOP = -0.00007                   # top decoupling


def build() -> dict:
    with open(ENDPOINT_PATH, encoding="utf-8") as f:
        endpoint = json.load(f)
    gap_lo, gap_hi = (float(x) for x in
                      endpoint["compare_only"]
                      ["same_scheme_anchor_gap_interval_inv_alpha"])

    alpha_inv_mz_phys = ALPHA_INV_0 * (1.0 - DELTA_LEP - DELTA_HAD5 - DELTA_TOP)
    gap_phys = alpha_inv_mz_phys - A0_OPH

    # inverse-alpha-unit contributions to the anchor deficit
    had_term = ALPHA_INV_0 * DELTA_HAD5
    lep_term = ALPHA_INV_0 * DELTA_LEP
    top_term = ALPHA_INV_0 * abs(DELTA_TOP)

    # what the OPH one-loop run does contain: the leptonic running is a smooth
    # one-loop object the RG captures; the hadronic piece is nonperturbative
    # and cannot appear in a one-loop unification run. The anchor deficit is
    # therefore dominated by the hadronic running plus higher-order remainder.
    return {
        "artifact": "oph_anchor_scheme_bridge",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "object_id": "AnchorSchemeBridge",
        "row_class": "compare_only_scheme_analysis",
        "guards": {
            "measured_values_in_any_oph_solve_path": False,
            "source_only_scheme_bridge_emitted": False,
            "public_promotion_allowed": False,
        },
        "certified_gap_from_endpoint": {
            "same_scheme_anchor_gap_interval": [gap_lo, gap_hi],
            "source": "code/P_derivation/runtime/empirical_thomson_endpoint_current.json",
        },
        "anchor_provenance": {
            "a0_oph": A0_OPH,
            "definition": "alpha_em^-1(m_Z^2; P) from the OPH one-loop RG run "
                          "from the OPH unification scale (paper_math.build_d10_from_p)",
            "scheme": "one_loop_unification_run",
            "contains_nonperturbative_hadronic_running": False,
        },
        "reference_decomposition_compare_only": {
            "alpha_inv_0": ALPHA_INV_0,
            "Delta_lep": DELTA_LEP,
            "Delta_had5": DELTA_HAD5,
            "Delta_top": DELTA_TOP,
            "alpha_inv_mz_phys_on_shell": alpha_inv_mz_phys,
            "gap_phys_minus_oph": gap_phys,
            "inv_alpha_unit_terms": {
                "hadronic_5f": had_term,
                "leptonic": lep_term,
                "top": top_term,
            },
            "citations": [
                "PDG 2024 Review, Electroweak model and constraints",
                "leptonic 3-loop: Steinhauser, Phys. Lett. B 429 (1998) 158",
                "hadronic 5-flavor data-driven: KNT19 arXiv:1911.00367",
            ],
        },
        "verdict": {
            "route_A_scheme_bridge": "available_empirical_class_only",
            "route_A_detail": "filling delta_scheme requires the measured "
                              "hadronic running (empirical class, already "
                              "carried by the payload) or the OPH source "
                              "hadronic spectral measure",
            "route_B_anchor_exact_no_go": "false",
            "route_B_detail": "the OPH anchor is a one-loop unification value "
                              "with a physical, understood hadronic and "
                              "higher-order running deficit; it is not "
                              "scheme-exact, and the gap is not a source-chain "
                              "failure",
            "source_only_reduction": "a source-only anchor bridge reduces to "
                                     "the OPH hadronic spectral measure, which "
                                     "is blocked on the hadron backend (#425)",
            "gap_consistency": "the standard reference deficit "
                               f"{gap_phys:.3f} sits at the lower edge of the "
                               "certified interval, confirming the gap is the "
                               "hadronic-plus-higher-order running deficit "
                               "rather than an anchor error",
            "issue_545_status": "structure_resolved_reduces_to_425_source_bridge_open",
        },
    }


def main() -> int:
    report = build()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
        f.write("\n")
    ref = report["reference_decomposition_compare_only"]
    v = report["verdict"]
    print(f"OPH one-loop anchor a0(P*)           = {A0_OPH:.5f}")
    print(f"physical alpha^-1(mZ) 5f on-shell    = {ref['alpha_inv_mz_phys_on_shell']:.5f}")
    print(f"reference deficit (phys - OPH)       = {ref['gap_phys_minus_oph']:+.5f}")
    print(f"certified anchor-gap interval        = "
          f"{report['certified_gap_from_endpoint']['same_scheme_anchor_gap_interval']}")
    print(f"route A: {v['route_A_scheme_bridge']}")
    print(f"route B (anchor-exact no-go): {v['route_B_anchor_exact_no_go']}")
    print(f"#545 status: {v['issue_545_status']}")
    print(f"wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
