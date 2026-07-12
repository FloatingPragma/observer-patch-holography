#!/usr/bin/env python3
"""Three-scalar interface theorem for the quark kernel program (#377).

Theorem (three-scalar interface). Fix the shared absolute scale g_ch. The
entire quark sector of the candidate chain factors through exactly the
invariant triple

    (r, sigma_u, sigma_d) in (R_>0)^3,

where r is the raw gap ratio of the centered compressed branch generator
and sigma_u, sigma_d are the per-side sector spans. Explicitly:

  (i)   rho_ord = 3/(2 + r) and x2 = (r - 1)/(r + 1), so the ordered ratio
        constant and the mean-law coordinate carry no information beyond r;
  (ii)  the rays v_u, v_d, the mean-law coefficients A_ud, B_ud, the
        sector means g_u, g_d, and the six masses are closed-form
        functions of the triple and g_ch;
  (iii) the map (r, sigma_u, sigma_d) -> (m_u, ..., m_b) is injective,
        with the explicit left inverse

            sigma_u = ln(m_t/m_u)/2,      sigma_d = ln(m_b/m_d)/2,
            r solves  ln(m_c/m_u)/ln(m_t/m_c) = rho/(1 - rho + rho^2/...)

        realized below as: rho = up-sector gap ratio over its own
        complement, r = 3/rho - 2, so the triple is recoverable from the
        sextet in closed form;
  (iv)  no proper sub-tuple suffices: by the kernel admissibility freedom
        theorem the certificate battery leaves the triple free, and the
        axiom-level nonidentifiability theorem shows the axioms emit no
        component of it.

Consequently the K1-K3 program of issue #377 is exactly the derivation of
one point of (R_>0)^3 with clean ancestry; everything downstream is
already a closed computation, and the acceptance harness scores any
candidate point mechanically.

Proof. (i) is algebra from the gap decomposition (g21, g32) =
(s r/(1+r), s/(1+r)). (ii) collects the declared maps of the mean-split
and edge artifacts; each equation is recomputed here against the on-disk
values to tolerance 1e-9 or better. (iii) the left inverse is exhibited
and the round trip triple -> sextet -> triple is executed at the
operating point and at off-band points to 1e-12. (iv) cites the two
obstruction artifacts. Every reference mass appears only inside the
round-trip verification of the left inverse and is labeled compare-only;
the theorem itself emits no numerical triple.

Run:
    python3 code/particles/flavor/derive_quark_kernel_three_scalar_interface_theorem.py
writes code/particles/runs/flavor/quark_kernel_three_scalar_interface_theorem.json.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
from datetime import datetime, timezone

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = ROOT / "particles" / "runs" / "flavor"
DEFAULT_OUT = RUNS / "quark_kernel_three_scalar_interface_theorem.json"
GENERATOR_PATH = RUNS / "generation_bundle_branch_generator.json"
MEAN_SPLIT_PATH = RUNS / "quark_sector_mean_split.json"
WRITEBACK_PATH = RUNS / "charged_shared_absolute_scale_writeback.json"
EDGE_PATH = RUNS / "quark_edge_statistics_spread_candidate.json"

UP_LABELS = ("m_u", "m_c", "m_t")
DOWN_LABELS = ("m_d", "m_s", "m_b")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def rho_of_r(r: float) -> float:
    return 3.0 / (2.0 + r)


def r_of_rho(rho: float) -> float:
    return 3.0 / rho - 2.0


def x2_of_r(r: float) -> float:
    return (r - 1.0) / (r + 1.0)


def rays(rho: float) -> tuple[np.ndarray, np.ndarray]:
    denom = 3.0 * (1.0 + rho)
    v_u = np.asarray([-(2.0 * rho + 1.0), rho - 1.0, rho + 2.0]) / denom
    v_d = np.asarray([-(rho + 2.0), 1.0 - rho, 2.0 * rho + 1.0]) / denom
    return v_u, v_d


def mean_coefficients(r: float) -> tuple[float, float]:
    rho, x2 = rho_of_r(r), x2_of_r(r)
    a_ud = 1.0 / (2.0 * (1.0 + rho - x2 * x2))
    b_ud = 1.0 / (2.0 * (1.0 - x2 * x2 - x2 * x2 / (1.0 + rho)))
    return a_ud, b_ud


def forward_sextet(r: float, sigma_u: float, sigma_d: float,
                   g_ch: float) -> dict[str, float]:
    rho = rho_of_r(r)
    a_ud, b_ud = mean_coefficients(r)
    seed = 0.5 * (sigma_u + sigma_d)
    eta = 0.5 * (sigma_u - sigma_d)
    g_u = g_ch * math.exp(-(a_ud * seed - b_ud * eta))
    g_d = g_ch * math.exp(-(a_ud * seed + b_ud * eta))
    v_u, v_d = rays(rho)
    out = {}
    for labels, g, sigma, ray in ((UP_LABELS, g_u, sigma_u, v_u),
                                  (DOWN_LABELS, g_d, sigma_d, v_d)):
        for label, v in zip(labels, ray, strict=True):
            out[label] = g * math.exp(2.0 * sigma * float(v))
    return out


def recover_triple(sextet: dict[str, float]) -> tuple[float, float, float]:
    """Closed-form left inverse of forward_sextet at fixed g_ch."""
    sigma_u = 0.5 * math.log(sextet["m_t"] / sextet["m_u"])
    sigma_d = 0.5 * math.log(sextet["m_b"] / sextet["m_d"])
    gap_ratio_up = (math.log(sextet["m_c"] / sextet["m_u"])
                    / math.log(sextet["m_t"] / sextet["m_c"]))
    # the up ray has gap ratio rho, so rho = gap_ratio_up and r = 3/rho - 2
    return r_of_rho(gap_ratio_up), sigma_u, sigma_d


def build() -> dict:
    generator = json.loads(GENERATOR_PATH.read_text(encoding="utf-8"))
    mean_split = json.loads(MEAN_SPLIT_PATH.read_text(encoding="utf-8"))
    writeback = json.loads(WRITEBACK_PATH.read_text(encoding="utf-8"))
    edge = json.loads(EDGE_PATH.read_text(encoding="utf-8"))

    spectrum = sorted(
        float(x) for x in
        generator["charged_sector_response_operator_candidate"]["ordered_spectrum"])
    g21, g32 = spectrum[1] - spectrum[0], spectrum[2] - spectrum[1]
    r_repo = g21 / g32
    g_ch = float(writeback["stored_shared_absolute_scale"])

    # (i) and (ii): recomputation checks against the on-disk artifacts
    checks = {}
    checks["rho_of_r_vs_mean_split"] = abs(
        rho_of_r(r_repo) - float(mean_split["rho_ord"]))
    checks["x2_of_r_vs_mean_split"] = abs(
        x2_of_r(r_repo) - float(mean_split["normalized_coordinate_x2"]))
    a_ud, b_ud = mean_coefficients(r_repo)
    checks["A_ud_vs_mean_split"] = abs(a_ud - float(mean_split["A_ud_candidate"]))
    checks["B_ud_vs_mean_split"] = abs(b_ud - float(mean_split["B_ud_candidate"]))
    v_u, v_d = rays(rho_of_r(r_repo))
    checks["v_u_vs_edge_artifact"] = float(np.max(np.abs(
        v_u - np.asarray(edge["profile_rays"]["v_u"]))))
    checks["v_d_vs_edge_artifact"] = float(np.max(np.abs(
        v_d - np.asarray(edge["profile_rays"]["v_d"]))))
    sig_u_ms = float(mean_split["sigma_u_total_log_per_side"])
    sig_d_ms = float(mean_split["sigma_d_total_log_per_side"])
    seed = 0.5 * (sig_u_ms + sig_d_ms)
    eta = 0.5 * (sig_u_ms - sig_d_ms)
    checks["g_u_vs_mean_split"] = abs(
        g_ch * math.exp(-(a_ud * seed - b_ud * eta))
        - float(mean_split["g_u_candidate"]))
    checks["g_d_vs_mean_split"] = abs(
        g_ch * math.exp(-(a_ud * seed + b_ud * eta))
        - float(mean_split["g_d_candidate"]))
    if max(checks.values()) > 1.0e-9:
        raise AssertionError(f"interface determination checks failed: {checks}")

    # (iii): injectivity by explicit left inverse, round trips executed at
    # the operating point and at off-band points
    round_trips = []
    for triple in ((r_repo, sig_u_ms, sig_d_ms),
                   (0.05, 2.0, 7.0),
                   (5.0, 0.3, 0.9)):
        sextet = forward_sextet(*triple, g_ch=g_ch)
        recovered = recover_triple(sextet)
        err = max(abs(a - b) / abs(a) for a, b in zip(triple, recovered,
                                                      strict=True))
        round_trips.append({
            "triple": {"r": triple[0], "sigma_u": triple[1],
                       "sigma_d": triple[2]},
            "max_relative_roundtrip_error": err,
        })
    max_rt = max(row["max_relative_roundtrip_error"] for row in round_trips)
    if max_rt > 1.0e-12:
        raise AssertionError("interface round trip failed")

    return {
        "artifact": "oph_quark_kernel_three_scalar_interface_theorem",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "proof_status": "closed_interface_theorem",
        "claim_tier": "exact_reduction_of_kernel_program_to_three_scalars",
        "row_class": "theorem_grade_reduction",
        "guards": {
            "emits_numerical_triple": False,
            "reference_masses_role": "compare-only round-trip verification "
                                     "of the left inverse; no solve path",
            "public_promotion_allowed": False,
        },
        "theorem_statement": (
            "At fixed shared scale g_ch the quark sector of the candidate "
            "chain factors through exactly the invariant triple "
            "(r, sigma_u, sigma_d) in (R_>0)^3 of the centered compressed "
            "branch generator: rho_ord = 3/(2+r) and x2 = (r-1)/(r+1) "
            "carry no information beyond r; rays, mean-law coefficients, "
            "sector means, and the six masses are closed-form in the "
            "triple; the forward map is injective with an explicit "
            "closed-form left inverse; and by the freedom and "
            "nonidentifiability theorems no proper sub-tuple is "
            "determined by the current corpus. The kernel program is "
            "therefore exactly the source derivation of one point of "
            "(R_>0)^3 with clean ancestry."
        ),
        "interface": {
            "coordinates": ["r (raw gap ratio)", "sigma_u (per-side span)",
                            "sigma_d (per-side span)"],
            "derived_scalars": {
                "rho_ord": "3/(2 + r)",
                "x2": "(r - 1)/(r + 1)",
                "A_ud": "1/(2 (1 + rho - x2^2))",
                "B_ud": "1/(2 (1 - x2^2 - x2^2/(1 + rho)))",
            },
            "left_inverse": {
                "sigma_u": "ln(m_t/m_u)/2",
                "sigma_d": "ln(m_b/m_d)/2",
                "r": "3/rho - 2 with rho = ln(m_c/m_u)/ln(m_t/m_c)",
            },
            "shared_scale": {
                "g_ch": g_ch,
                "provenance": "charged shared absolute scale writeback "
                              "(current-family scope; fourth input with "
                              "its own lane, held fixed by this theorem)",
            },
        },
        "determination_checks_max_abs": checks,
        "round_trips": round_trips,
        "max_roundtrip_relative_error": max_rt,
        "nonredundancy": {
            "freedom_theorem": "family_transport_kernel_admissibility_"
                               "freedom_theorem (certificates leave the "
                               "triple free)",
            "axiom_level_no_go": "quark_axiom_level_yukawa_moduli_"
                                 "nonidentifiability (the axioms emit no "
                                 "component of the triple)",
            "conclusion": "the interface is exactly three-dimensional: "
                          "sufficient by the determination maps, "
                          "recoverable by the left inverse, and free on "
                          "the current corpus in every component",
        },
        "acceptance_binding": (
            "a candidate derivation closes the program by emitting one "
            "clean-ancestry point of the interface; the acceptance "
            "harness gates G2-G6 score exactly that point"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prove the three-scalar interface theorem for the "
                    "quark kernel program.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    report = build()
    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    print("determination checks (max abs deviation vs on-disk artifacts):")
    for name, value in report["determination_checks_max_abs"].items():
        print(f"  {name:28s} {value:.2e}")
    print(f"round-trip max relative error: "
          f"{report['max_roundtrip_relative_error']:.2e}")
    print(f"proof status: {report['proof_status']}")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
