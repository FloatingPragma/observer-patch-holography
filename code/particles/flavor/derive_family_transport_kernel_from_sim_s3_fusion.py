#!/usr/bin/env python3
"""Second K1 execution: S3 class-fusion transport kernel (#377, carrier #2).

Carrier candidate #1 (sector population deltas) was falsified by its
conservation null: class-count deltas sum to zero, so that readout has two
physical dimensions. This module executes the named successor readout, the
class-fusion transition structure of the realized gauge field, which has
no sum rule and is the group-algebra object of the same species as the
fusion and statistical-dimension structure behind the electroweak
coefficients.

Predeclaration (fixed before computation):

  carrier      the three S3 conjugacy classes (identity, transposition,
               threecycle) as the same-label sector basis; the
               identification with the generation bundle is a CANDIDATE
               CLAIM, not a theorem.
  transport    T[a,b] = P( class(g_jk g_ij) = b | class(g_ij) = a ) over
               every ordered composable directed-edge pair i->j->k with
               k != i in the FINAL gauge state, with the sim's own
               conventions mirrored exactly: directed label of a reversed
               edge is the group inverse, and two-step path composition is
               S3_MUL[g_ij, g_jk] as in s3_triangle_holonomy. T is the
               realized one-step sector transport of the gauge dynamics.
  descendant   hermitian descendant = T T^t, the schema's T T-dagger rule.
  levels       level 0 = e1_s3_bw_screen_4k, level 1 = e1_s3_bw_screen_64k,
               the same experiment family at two patch counts, fixed
               config seeds, both instrumented with the s3_gauge_state
               export; the branch generator reads the spectrum from the
               latest level.
  intertwiner  identity on the canonical class frame.
  scale        no rescaling; shape invariants (r, rho_ord, x2) are the
               milestone, spans stay in raw transition units until a K2
               clause exists.
  nulls        (a) independent-pair null: the same count assembled from
               global label marginals with node-level correlations
               removed; (b) analytic Haar null: uniform S3 labels, the
               pure class-multiplication-table prediction. Both computed
               at both levels. A physical rho indistinguishable from the
               Haar null means the dynamics add nothing to group
               combinatorics; the physical-vs-null contrast carries the
               evidential weight.
  verdict      level-1 rho_ord against the quoted-convention band of the
               common-scheme scan, plus the acceptance-harness gate
               record (G4 blocking is the expected fail-closed outcome:
               no normalization clause exists).

No quark mass, fitted spread, or flavor template enters any step.

Run (under the sim venv, which carries scipy for the imported module):
    oph-physics-sim/.venv/bin/python \
      code/particles/flavor/derive_family_transport_kernel_from_sim_s3_fusion.py \
      --level0-dir <4k run dir> --level1-dir <64k run dir>
writes code/particles/runs/flavor/family_transport_kernel_sim_s3_fusion_computed.json
and    code/particles/runs/flavor/family_transport_kernel_sim_s3_fusion_verdict.json.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
from datetime import datetime, timezone

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
RER_ROOT = HERE.parents[2]
SIM_ROOT = HERE.parents[3] / "oph-physics-sim"
RUNS = RER_ROOT / "code" / "particles" / "runs" / "flavor"
PAYLOAD_OUT = RUNS / "family_transport_kernel_sim_s3_fusion_computed.json"
VERDICT_OUT = RUNS / "family_transport_kernel_sim_s3_fusion_verdict.json"
SCAN_PATH = RUNS / "quark_common_scheme_shape_law_scan.json"

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(RER_ROOT / "code"))
from particles.artifact_paths import canonicalize_artifact_paths  # noqa: E402

from derive_generation_bundle_branch_generator import (  # noqa: E402
    build_artifact as build_generator,
)
from derive_quark_kernel_normalization_acceptance_harness import (  # noqa: E402
    evaluate as harness_evaluate,
)

_spec = importlib.util.spec_from_file_location(
    "array_s3_holonomy", SIM_ROOT / "oph_fpe" / "defects" / "array_s3_holonomy.py")
_s3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_s3)
S3_MUL, S3_INV, S3_CLASS = _s3.S3_MUL, _s3.S3_INV, _s3.S3_CLASS


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _encode(matrix: np.ndarray) -> dict:
    return {"real": np.real(matrix).tolist(), "imag": np.imag(matrix).tolist()}


def pair_label_counts(left: np.ndarray, right: np.ndarray,
                      gauge: np.ndarray, patch_count: int) -> dict:
    """Exact counts over ordered label pairs of composable directed edges.

    Directed edges: every undirected edge contributes (left->right, g) and
    (right->left, inv g). For node j, every in-edge (i->j) pairs with every
    out-edge (j->k); the single reversal pair (k = i) is excluded exactly.
    """
    src = np.concatenate([left, right])
    dst = np.concatenate([right, left])
    lab = np.concatenate([gauge, S3_INV[gauge]]).astype(np.int64)

    in_hist = np.zeros((patch_count, 6), dtype=np.int64)
    out_hist = np.zeros((patch_count, 6), dtype=np.int64)
    np.add.at(in_hist, (dst, lab), 1)
    np.add.at(out_hist, (src, lab), 1)

    c6 = in_hist.T @ out_hist
    # exact reversal exclusion: the reverse of in-edge (i->j, x) is the
    # out-edge (j->i, inv x), present exactly once per directed edge
    label_totals = np.bincount(lab, minlength=6)
    for x in range(6):
        c6[x, S3_INV[x]] -= label_totals[x]
    if int(c6.min()) < 0:
        raise AssertionError("negative pair count after reversal exclusion")
    return {"c6": c6, "in_hist": in_hist, "out_hist": out_hist,
            "label_totals": label_totals}


def class_fusion_from_c6(c6: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Collapse ordered label-pair counts to the 3x3 class-fusion matrix."""
    counts = np.zeros((3, 3), dtype=np.float64)
    for x in range(6):
        for y in range(6):
            counts[S3_CLASS[x], S3_CLASS[S3_MUL[x, y]]] += c6[x, y]
    rows = counts.sum(axis=1, keepdims=True)
    if np.any(rows <= 0):
        raise AssertionError("empty class row in fusion counts")
    return counts, counts / rows


def independent_null_c6(stats: dict) -> np.ndarray:
    """Node-decorrelated null: global in-label x out-label marginal product,
    scaled to the physical pair total."""
    total_in = stats["in_hist"].sum(axis=0).astype(np.float64)
    total_out = stats["out_hist"].sum(axis=0).astype(np.float64)
    outer = np.outer(total_in, total_out)
    physical_total = float(stats["c6"].sum())
    return outer * (physical_total / outer.sum())


def haar_null_c6() -> np.ndarray:
    return np.ones((6, 6), dtype=np.float64) / 36.0


def readout(descendant: np.ndarray) -> dict:
    centered = descendant - (np.trace(descendant) / 3.0) * np.eye(3)
    evals = np.linalg.eigvalsh(centered)
    g21 = float(evals[1] - evals[0])
    g32 = float(evals[2] - evals[1])
    r = g21 / g32 if g32 > 0 else float("inf")
    return {
        "eigenvalues_centered": [float(x) for x in evals],
        "g21": g21,
        "g32": g32,
        "span": float(evals[2] - evals[0]),
        "raw_gap_ratio_r": r,
        "rho_ord": 3.0 * g32 / (2.0 * g32 + g21) if (2.0 * g32 + g21) > 0
                   else float("nan"),
        "x2": (r - 1.0) / (r + 1.0),
    }


def level_block(run_dir: pathlib.Path) -> dict:
    state = np.load(run_dir / "s3_gauge_state.npz")
    left = np.asarray(state["left"], dtype=np.int64)
    right = np.asarray(state["right"], dtype=np.int64)
    gauge = np.asarray(state["gauge"], dtype=np.int64)
    patch_count = int(np.asarray(state["points"]).shape[0])

    stats = pair_label_counts(left, right, gauge, patch_count)
    rows = {}
    for name, c6 in (("physical", stats["c6"].astype(np.float64)),
                     ("independent_pair_null", independent_null_c6(stats)),
                     ("haar_null", haar_null_c6())):
        counts, t = class_fusion_from_c6(c6)
        d = t @ t.T
        rows[name] = {
            "transition_matrix_row_stochastic": t.tolist(),
            "class_counts": counts.tolist(),
            **readout(d),
        }
        rows[name]["descendant"] = d
    return {
        "run_dir": str(run_dir),
        "patch_count": patch_count,
        "edge_count": int(gauge.shape[0]),
        "pair_total": int(stats["c6"].sum()),
        "rows": rows,
    }


def build(level0_dir: pathlib.Path, level1_dir: pathlib.Path) -> tuple:
    level0 = level_block(level0_dir)
    level1 = level_block(level1_dir)

    d0 = level0["rows"]["physical"].pop("descendant")
    d1 = level1["rows"]["physical"].pop("descendant")
    for block in (level0, level1):
        for row in block["rows"].values():
            row.pop("descendant", None)

    scan = json.loads(SCAN_PATH.read_text(encoding="utf-8"))
    band = scan["findings"]["a2_band_by_convention"]["quoted_mixed"]
    template_rho = float(scan["template_rho_ord"])
    rho1 = level1["rows"]["physical"]["rho_ord"]

    payload = {
        "artifact": "oph_family_transport_kernel",
        "generated_utc": _timestamp(),
        "status": "computed_sim_s3_fusion_candidate",
        "transport_kind": "s3_class_fusion_one_step_edge_transport",
        "proof_status": "computed_candidate_carrier_identification_open",
        "carrier_identification": {
            "claim": "the three S3 conjugacy classes under realized "
                     "one-step fusion transport are the same-label "
                     "generation carrier",
            "status": "candidate_claim_not_theorem",
            "sum_rule_note": "no conservation null: rows of the fusion "
                             "transition are genuine three-dimensional "
                             "distributions",
        },
        "refinements": [
            {"level": 0, "run_dir": level0["run_dir"],
             "patch_count": level0["patch_count"],
             "transport_operator": _encode(np.asarray(
                 level0["rows"]["physical"]["transition_matrix_row_stochastic"],
                 dtype=complex)),
             "hermitian_descendant": _encode(d0.astype(complex))},
            {"level": 1, "run_dir": level1["run_dir"],
             "patch_count": level1["patch_count"],
             "transport_operator": _encode(np.asarray(
                 level1["rows"]["physical"]["transition_matrix_row_stochastic"],
                 dtype=complex)),
             "hermitian_descendant": _encode(d1.astype(complex))},
        ],
        "refinement_intertwiners": [
            {"from_level": 0, "to_level": 1,
             "intertwiner": _encode(np.eye(3, dtype=complex)),
             "justification": "canonical class frame"},
        ],
        "ancestry": {
            "artifacts": [
                f"s3_gauge_state({level0['run_dir'].split('/')[-1]})",
                f"s3_gauge_state({level1['run_dir'].split('/')[-1]})",
            ],
            "attestations": {
                "quark_reference_values_consumed": False,
                "fitted_spreads_consumed": False,
                "numerical_flavor_template_consumed": False,
            },
        },
        "metadata": {
            "note": "Second K1 execution: realized S3 class-fusion "
                    "transport at two patch counts of the same experiment "
                    "family. Carrier identification and the K2 scale "
                    "clause stay open.",
        },
    }

    generator = build_generator(payload)
    gates = harness_evaluate({
        "candidate_id": "sim_s3_class_fusion_kernel",
        "kernel": {"refinements": payload["refinements"]},
        "ancestry": payload["ancestry"],
    })

    verdict = {
        "artifact": "oph_family_transport_kernel_sim_s3_fusion_verdict",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "row_class": "computed_candidate_shape_milestone",
        "guards": {
            "quark_reference_values_consumed": False,
            "fitted_spreads_consumed": False,
            "numerical_flavor_template_consumed": False,
            "public_promotion_allowed": False,
            "post_hoc_construction_changes": "none; predeclared in the "
                                             "module docstring before "
                                             "computation",
        },
        "levels": {"level0": level0, "level1": level1},
        "band_check": {
            "declared_convention_band": band,
            "template_rho_ord": template_rho,
            "computed_rho_ord_level1": rho1,
            "inside_band": bool(band[0] <= rho1 <= band[1]),
            "physical_minus_template": rho1 - template_rho,
        },
        "null_contrast_level1": {
            name: level1["rows"][name]["rho_ord"]
            for name in ("physical", "independent_pair_null", "haar_null")
        },
        "harness_gates": gates,
        "generator_proof_status": generator.get("proof_status"),
        "reading_rules": [
            "the Haar null is the pure class-multiplication prediction; "
            "the physical-vs-Haar contrast measures what the realized "
            "dynamics add beyond group combinatorics",
            "G4 blocking is the expected fail-closed outcome until a K2 "
            "normalization clause exists",
        ],
    }
    return payload, verdict, generator


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute the S3 class-fusion transport kernel from "
                    "instrumented sim runs.")
    parser.add_argument("--level0-dir", required=True)
    parser.add_argument("--level1-dir", required=True)
    parser.add_argument("--payload-out", default=str(PAYLOAD_OUT))
    parser.add_argument("--verdict-out", default=str(VERDICT_OUT))
    args = parser.parse_args()

    payload, verdict, generator = build(
        pathlib.Path(args.level0_dir),
        pathlib.Path(args.level1_dir),
    )
    payload = canonicalize_artifact_paths(payload)
    verdict = canonicalize_artifact_paths(verdict)
    generator = canonicalize_artifact_paths(generator)
    for path, blob in ((pathlib.Path(args.payload_out), payload),
                       (pathlib.Path(args.verdict_out), verdict)):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(blob, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    gen_path = RUNS / "generation_bundle_branch_generator_sim_s3_fusion.json"
    gen_path.write_text(json.dumps(generator, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    for label in ("level0", "level1"):
        block = verdict["levels"][label]
        row = block["rows"]["physical"]
        print(f"{label} ({block['run_dir'].split('/')[-1]}, "
              f"{block['patch_count']} patches, {block['pair_total']} pairs):")
        print(f"  T = {np.round(np.asarray(row['transition_matrix_row_stochastic']), 5).tolist()}")
        print(f"  r = {row['raw_gap_ratio_r']:.6f}  rho_ord = {row['rho_ord']:.6f}  "
              f"x2 = {row['x2']:+.6f}")
    bc = verdict["band_check"]
    print(f"band {bc['declared_convention_band']}  template {bc['template_rho_ord']:.6f}")
    print(f"COMPUTED rho_ord = {bc['computed_rho_ord_level1']:.6f}  "
          f"inside_band = {bc['inside_band']}  "
          f"delta_vs_template = {bc['physical_minus_template']:+.6f}")
    print("null contrast (level 1 rho_ord):")
    for name, value in verdict["null_contrast_level1"].items():
        print(f"  {name:24s} {value:.6f}")
    print(f"harness: first_failed_gate = {verdict['harness_gates']['first_failed_gate']}  "
          f"status = {verdict['harness_gates']['status']}")
    print(f"saved: {args.payload_out}")
    print(f"saved: {args.verdict_out}")
    print(f"saved: {gen_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
