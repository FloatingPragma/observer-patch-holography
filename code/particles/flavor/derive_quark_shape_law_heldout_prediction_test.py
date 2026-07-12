#!/usr/bin/env python3
"""Held-out quark mass tests of the frozen shape law (#377, #379, #380).

The axiom-level nonidentifiability theorem stands: Axioms 1-5 plus fixed P
do not define the quark spread moduli, and this artifact does not claim
otherwise. It runs the strongest intermediate test available on the current
corpus: freeze the candidate structure (profile rays from rho_ord, the
two-scalar affine mean law, the shared absolute scale g_ch), calibrate at
most two quark observables, and score every remaining mass as a held-out
check against the compare-only references.

Declared model. With per-side spreads sigma_q and unit-span trace-zero rays
v_q from the ordered ratio law,

    v_u = (-(2 rho + 1), rho - 1, rho + 2) / (3 (1 + rho)),
    v_d = (-(rho + 2), 1 - rho, 2 rho + 1) / (3 (1 + rho)),
    ln m_{q,i} = ln g_q + 2 sigma_q v_{q,i},

and the two-scalar affine mean law of the mean-split artifact,

    ln g_u = ln g_ch - (A_ud sigma_seed - B_ud eta_ud),
    ln g_d = ln g_ch - (A_ud sigma_seed + B_ud eta_ud),
    sigma_seed = (sigma_u + sigma_d)/2,  eta_ud = (sigma_u - sigma_d)/2.

Four tests, in decreasing strictness of input economy:

  reciprocity (zero parameters): the ray law forces the product of the two
      sector log-gap ratios to equal exactly one; the data product is the
      held-out score.
  version A (four inputs, two held out): sector endpoints m_u, m_t, m_d,
      m_b calibrate (sigma_q, g_q) per sector; the interior masses m_c and
      m_s are predicted.
  version B (two inputs, four held out): only m_t and m_b are consumed;
      the mean law plus g_ch closes the sector means, so (sigma_u, sigma_d)
      solve a linear 2x2 system and m_u, m_c, m_d, m_s are predicted.
  version C (zero quark inputs, six displayed): the edge-statistics
      candidate spreads replace calibration entirely; all six masses are
      emitted and scored.

Status. Everything here is conditional on the template transport kernel
(see the rho_ord sensitivity audit), the edge-coefficient ansatz (version C
spreads), and the current-family shared-scale writeback (g_ch). No version
is source-only, promotion stays blocked, and the nonidentifiability theorem
is unaffected: versions A and B calibrate the free moduli from predeclared
data, version C fills them from a candidate ansatz. Reference values are
consumed only as the predeclared calibration inputs and as compare-only
held-out scores. Reference conventions are the mixed schemes recorded in
particle_reference_values.json; the common-scheme rerun is the declared
follow-up.

Run:
    python3 code/particles/flavor/derive_quark_shape_law_heldout_prediction_test.py
writes code/particles/runs/flavor/quark_shape_law_heldout_prediction_test.json.
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
DEFAULT_OUT = RUNS / "quark_shape_law_heldout_prediction_test.json"
GENERATOR_PATH = RUNS / "generation_bundle_branch_generator.json"
MEAN_SPLIT_PATH = RUNS / "quark_sector_mean_split.json"
WRITEBACK_PATH = RUNS / "charged_shared_absolute_scale_writeback.json"
EDGE_PATH = RUNS / "quark_edge_statistics_spread_candidate.json"
REFERENCE_PATH = ROOT / "particles" / "data" / "particle_reference_values.json"

UP_KEYS = ("up_quark", "charm_quark", "top_quark")
DOWN_KEYS = ("down_quark", "strange_quark", "bottom_quark")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def rays(rho: float) -> tuple[np.ndarray, np.ndarray]:
    denom = 3.0 * (1.0 + rho)
    v_u = np.asarray([-(2.0 * rho + 1.0), rho - 1.0, rho + 2.0]) / denom
    v_d = np.asarray([-(rho + 2.0), 1.0 - rho, 2.0 * rho + 1.0]) / denom
    return v_u, v_d


def forward(g: float, sigma: float, ray: np.ndarray) -> list[float]:
    return [float(g * math.exp(2.0 * sigma * v)) for v in ray]


def comparison_rows(predicted: dict[str, float], refs: dict,
                    held_out: tuple[str, ...],
                    spans: dict[str, float]) -> dict:
    rows = {}
    for key, value in predicted.items():
        entry = refs[key]
        ref = float(entry["value_gev"])
        sigma_ref = 0.5 * (float(entry.get("error_plus_gev") or 0.0)
                           + float(entry.get("error_minus_gev") or 0.0))
        dlog = math.log(value / ref)
        sector_span = spans["up" if key in UP_KEYS else "down"]
        rows[key] = {
            "predicted_gev": value,
            "reference_gev": ref,
            "reference_sigma_gev": sigma_ref,
            "relative_error": value / ref - 1.0,
            "abs_dlog": abs(dlog),
            "dlog_over_sector_log_span": abs(dlog) / sector_span,
            "delta_over_reference_sigma": ((value - ref) / sigma_ref
                                           if sigma_ref > 0.0 else None),
            "role": "held_out" if key in held_out else "calibration_input",
        }
    return rows


def build() -> dict:
    generator = json.loads(GENERATOR_PATH.read_text(encoding="utf-8"))
    mean_split = json.loads(MEAN_SPLIT_PATH.read_text(encoding="utf-8"))
    writeback = json.loads(WRITEBACK_PATH.read_text(encoding="utf-8"))
    edge = json.loads(EDGE_PATH.read_text(encoding="utf-8"))
    refs = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))["entries"]

    # frozen structure with recomputation checks against the artifacts
    matrix = (np.asarray(generator["centered_compressed_branch_generator"]["real"])
              + 1j * np.asarray(generator["centered_compressed_branch_generator"]["imag"]))
    evals = np.linalg.eigvalsh(matrix)
    g21, g32 = float(evals[1] - evals[0]), float(evals[2] - evals[1])
    rho = 3.0 * g32 / (2.0 * g32 + g21)
    x2 = 2.0 * g21 / float(evals[2] - evals[0]) - 1.0
    if abs(rho - float(mean_split["rho_ord"])) > 1.0e-12:
        raise AssertionError("rho_ord disagrees between generator and mean split")
    if abs(x2 - float(mean_split["normalized_coordinate_x2"])) > 1.0e-12:
        raise AssertionError("x2 disagrees between generator and mean split")

    a_ud = 1.0 / (2.0 * (1.0 + rho - x2 * x2))
    b_ud = 1.0 / (2.0 * (1.0 - x2 * x2 - x2 * x2 / (1.0 + rho)))
    if abs(a_ud - float(mean_split["A_ud_candidate"])) > 1.0e-12:
        raise AssertionError("A_ud recomputation disagrees with mean split")
    if abs(b_ud - float(mean_split["B_ud_candidate"])) > 1.0e-12:
        raise AssertionError("B_ud recomputation disagrees with mean split")

    g_ch = float(writeback["stored_shared_absolute_scale"])
    v_u, v_d = rays(rho)
    for mine, recorded in ((v_u, edge["profile_rays"]["v_u"]),
                           (v_d, edge["profile_rays"]["v_d"])):
        if float(np.max(np.abs(mine - np.asarray(recorded)))) > 1.0e-12:
            raise AssertionError("ray formulas disagree with edge artifact")

    def mean_law(sigma_u: float, sigma_d: float) -> tuple[float, float]:
        seed = 0.5 * (sigma_u + sigma_d)
        eta = 0.5 * (sigma_u - sigma_d)
        return (g_ch * math.exp(-(a_ud * seed - b_ud * eta)),
                g_ch * math.exp(-(a_ud * seed + b_ud * eta)))

    # mean-law validation on the recorded diagnostic spreads
    g_u_check, g_d_check = mean_law(
        float(mean_split["sigma_u_total_log_per_side"]),
        float(mean_split["sigma_d_total_log_per_side"]))
    if abs(g_u_check - float(mean_split["g_u_candidate"])) > 1.0e-9:
        raise AssertionError("mean-law recomputation disagrees with g_u_candidate")
    if abs(g_d_check - float(mean_split["g_d_candidate"])) > 1.0e-9:
        raise AssertionError("mean-law recomputation disagrees with g_d_candidate")

    def mass(key: str) -> float:
        return float(refs[key]["value_gev"])

    m = {key: mass(key) for key in UP_KEYS + DOWN_KEYS}
    spans = {
        "up": math.log(m["top_quark"] / m["up_quark"]),
        "down": math.log(m["bottom_quark"] / m["down_quark"]),
    }

    # zero-parameter reciprocity check
    up_ratio = math.log(m["charm_quark"] / m["up_quark"]) / math.log(
        m["top_quark"] / m["charm_quark"])
    down_ratio = math.log(m["strange_quark"] / m["down_quark"]) / math.log(
        m["bottom_quark"] / m["strange_quark"])
    reciprocity = {
        "law": "the ray pair forces gap_ratio_up * gap_ratio_down = 1 exactly, "
               "with no free parameter and no template dependence",
        "data_gap_ratio_up": up_ratio,
        "data_gap_ratio_down": down_ratio,
        "data_product": up_ratio * down_ratio,
        "law_product": 1.0,
        "relative_deviation": up_ratio * down_ratio - 1.0,
    }

    # version A: sector endpoints in, interior masses out
    sigma_u_a = 0.5 * math.log(m["top_quark"] / m["up_quark"])
    sigma_d_a = 0.5 * math.log(m["bottom_quark"] / m["down_quark"])
    g_u_a = m["up_quark"] * math.exp(-2.0 * sigma_u_a * float(v_u[0]))
    g_d_a = m["down_quark"] * math.exp(-2.0 * sigma_d_a * float(v_d[0]))
    pred_a = dict(zip(UP_KEYS, forward(g_u_a, sigma_u_a, v_u), strict=True))
    pred_a.update(zip(DOWN_KEYS, forward(g_d_a, sigma_d_a, v_d), strict=True))
    version_a = {
        "inputs_declared": ["up_quark", "top_quark", "down_quark", "bottom_quark"],
        "held_out": ["charm_quark", "strange_quark"],
        "calibrated": {"sigma_u": sigma_u_a, "sigma_d": sigma_d_a,
                       "g_u": g_u_a, "g_d": g_d_a},
        "uses_mean_law": False,
        "comparison": comparison_rows(
            pred_a, refs, ("charm_quark", "strange_quark"), spans),
    }

    # version B: heavy endpoints in, four light masses out (linear solve)
    half_diff = 0.5 * (a_ud - b_ud)
    half_sum = 0.5 * (a_ud + b_ud)
    system = np.asarray([
        [2.0 * float(v_u[2]) - half_diff, -half_sum],
        [-half_sum, 2.0 * float(v_d[2]) - half_diff],
    ])
    rhs = np.asarray([
        math.log(m["top_quark"] / g_ch),
        math.log(m["bottom_quark"] / g_ch),
    ])
    sigma_u_b, sigma_d_b = (float(x) for x in np.linalg.solve(system, rhs))
    g_u_b, g_d_b = mean_law(sigma_u_b, sigma_d_b)
    pred_b = dict(zip(UP_KEYS, forward(g_u_b, sigma_u_b, v_u), strict=True))
    pred_b.update(zip(DOWN_KEYS, forward(g_d_b, sigma_d_b, v_d), strict=True))
    for key in ("top_quark", "bottom_quark"):
        if abs(pred_b[key] / m[key] - 1.0) > 1.0e-10:
            raise AssertionError("version B fails to reproduce its inputs")
    version_b = {
        "inputs_declared": ["top_quark", "bottom_quark"],
        "held_out": ["up_quark", "charm_quark", "down_quark", "strange_quark"],
        "calibrated": {"sigma_u": sigma_u_b, "sigma_d": sigma_d_b,
                       "g_u": g_u_b, "g_d": g_d_b},
        "uses_mean_law": True,
        "mean_law_conditionality": "g_ch from the current-family shared-scale "
                                   "writeback; A_ud, B_ud from rho_ord and x2",
        "comparison": comparison_rows(
            pred_b, refs,
            ("up_quark", "charm_quark", "down_quark", "strange_quark"), spans),
    }

    # version C: zero quark-mass inputs, candidate spreads from the edge packet
    sigma_u_c = float(edge["candidate_sigmas"]["sigma_u_total_log_per_side"])
    sigma_d_c = float(edge["candidate_sigmas"]["sigma_d_total_log_per_side"])
    g_u_c, g_d_c = mean_law(sigma_u_c, sigma_d_c)
    pred_c = dict(zip(UP_KEYS, forward(g_u_c, sigma_u_c, v_u), strict=True))
    pred_c.update(zip(DOWN_KEYS, forward(g_d_c, sigma_d_c, v_d), strict=True))
    version_c = {
        "inputs_declared": [],
        "held_out": list(UP_KEYS + DOWN_KEYS),
        "spread_source": {
            "artifact": edge.get("artifact"),
            "bridge_status": edge.get("bridge_status"),
            "sigma_u": sigma_u_c,
            "sigma_d": sigma_d_c,
            "note": "edge-coefficient combination is a declared ansatz on "
                    "the cocycle suppression data; see the "
                    "nonidentifiability theorem for why no source rule "
                    "fixes these coefficients on the current corpus",
        },
        "uses_mean_law": True,
        "comparison": comparison_rows(pred_c, refs, UP_KEYS + DOWN_KEYS, spans),
    }
    max_abs_rel_c = max(abs(row["relative_error"])
                        for row in version_c["comparison"].values())
    version_c["max_abs_relative_error"] = max_abs_rel_c
    version_c["max_dlog_over_sector_span"] = max(
        row["dlog_over_sector_log_span"]
        for row in version_c["comparison"].values())

    return {
        "artifact": "oph_quark_shape_law_heldout_prediction_test",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "row_class": "template_and_ansatz_conditional_heldout_test",
        "guards": {
            "source_only_prediction": False,
            "public_promotion_allowed": False,
            "nonidentifiability_theorem_unaffected": True,
            "measured_values_in_solve_paths": "only the predeclared "
                                              "calibration inputs listed per "
                                              "version",
        },
        "frozen_structure": {
            "rho_ord": rho,
            "x2": x2,
            "A_ud": a_ud,
            "B_ud": b_ud,
            "g_ch": g_ch,
            "v_u": [float(v) for v in v_u],
            "v_d": [float(v) for v in v_d],
            "provenance": {
                "rho_ord_x2": "centered compressed branch generator "
                              "(candidate, template kernel ancestry; see "
                              "family_transport_kernel_rho_ord_sensitivity_audit)",
                "g_ch": "charged shared absolute scale writeback "
                        "(current-family scope)",
                "mean_law": "quark sector mean split (A_ud, B_ud closed in "
                            "rho_ord and x2)",
            },
        },
        "reference_policy": {
            "source": "particle_reference_values.json (PDG 2026 API rows)",
            "scheme_note": "mixed conventions as recorded there (light "
                           "quarks MS-bar 2 GeV class, charm and bottom at "
                           "their own scale, top from cross sections); the "
                           "common-scheme rerun is the declared follow-up",
        },
        "reciprocity_zero_parameter_check": reciprocity,
        "version_a_endpoints_in_interior_out": version_a,
        "version_b_two_in_four_out": version_b,
        "version_c_zero_quark_inputs": version_c,
        "notes": [
            "This artifact is the intermediate scientific test between the "
            "closed nonidentifiability theorem and a future source-only "
            "derivation: it measures how much of the quark spectrum the "
            "frozen candidate structure carries once the two free moduli "
            "are supplied by predeclared calibration or by the edge ansatz.",
            "No version emits a source-only quark mass and none is a "
            "public prediction row; the display class is a conditional "
            "candidate audit.",
            "The template dependence of rho_ord is quantified in the "
            "companion sensitivity audit; the issue-377 closure program "
            "(OPH-derived kernel with a normalized spectral scale) is the "
            "declared route to removing every conditional layer at once.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the held-out quark shape-law prediction tests.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    report = build()
    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    rec = report["reciprocity_zero_parameter_check"]
    print(f"reciprocity (0 params): data product = {rec['data_product']:.5f} "
          f"vs law 1.00000  ({rec['relative_deviation']:+.2%})")
    for label, key in (("A (4 in, 2 out)", "version_a_endpoints_in_interior_out"),
                       ("B (2 in, 4 out)", "version_b_two_in_four_out"),
                       ("C (0 quark inputs)", "version_c_zero_quark_inputs")):
        block = report[key]
        print(f"version {label}:")
        for name, row in block["comparison"].items():
            if row["role"] != "held_out":
                continue
            print(f"  {name:14s} pred {row['predicted_gev']:.6g} GeV  "
                  f"ref {row['reference_gev']:.6g}  "
                  f"({row['relative_error']:+.2%}, "
                  f"{row['dlog_over_sector_log_span']:.4f} of sector log span)")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
