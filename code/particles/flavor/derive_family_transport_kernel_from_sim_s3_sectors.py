#!/usr/bin/env python3
"""First computed family-transport kernel from simulated S3 sector data (#377 K1).

This is the first execution of the K1 computation: a family-transport
kernel whose entries come from realized simulation output instead of the
hand-written template. Everything analysis-defining is predeclared here,
before any spectrum is computed, and the verdict is reported as it comes
out.

Predeclaration (fixed before computation):

  carrier      the three S3 conjugacy-class response channels
               s3_sector_class:class_delta_{0,1,2} of the modular response
               kernel cache. This is the unique exactly-three-channel
               field family in the cache (checkpoint_class has 4,
               stable_flag has 2, record_family has 16), which is what
               makes it the candidate same-label generation carrier. The
               identification of these sector channels with the
               generation bundle is a CANDIDATE CLAIM, not a theorem.
  transport    T = the cache's stored response matrix restricted to the
               three class-delta columns, stacked over (observer, cap,
               time) rows; the stored matrix is the cache's canonical
               signed-robust-zscore object and the only physical response
               array on disk, so the per-column transform is a declared
               cache property shared by all three channels.
  descendant   D = T^t T, the 3x3 second-moment compression, matching the
               kernel schema's T T-dagger hermitian-descendant rule.
  levels       level 0 = runs/oph_universe_4k_viewer_smoke_20260706_v2
               (coarse patch count), level 1 =
               runs/oph_universe_64k_final_audited_20260711 (fine, the
               audited run); the branch-generator convention reads the
               spectrum from the latest level.
  intertwiner  identity on the canonical class frame (identity,
               transposition, threecycle are absolute labels across runs).
  scale        D is not rescaled; the K2 normalization clause does not
               exist yet, so the SHAPE invariants (r, rho_ord, x2) are
               the milestone and the span is recorded in raw cache units.
  nulls        the same construction on every control array of the level-1
               payload (s2 boundary wrong-group, shuffled, shuffled
               response, shuffled observer labels, no modular flow, wrong
               scale); the physical-vs-null contrast is part of the
               verdict, and a physical rho indistinguishable from the
               shuffled ensemble weakens any band coincidence.
  verdict      the emitted rho_ord against the quoted-convention band
               [1.2761, 1.3000] of the common-scheme scan, plus the full
               acceptance-harness gate record (G4 is expected to block:
               no normalization clause exists, which is the correct
               fail-closed outcome for a shape-only milestone).

No quark mass, fitted spread, or flavor template enters any step; the
ancestry is simulation output plus this predeclared compression.

Run:
    python3 code/particles/flavor/derive_family_transport_kernel_from_sim_s3_sectors.py
writes code/particles/runs/flavor/family_transport_kernel_sim_s3_computed.json (payload)
and    code/particles/runs/flavor/family_transport_kernel_sim_s3_shape_verdict.json.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from datetime import datetime, timezone

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
RER_ROOT = HERE.parents[2]
SIM_ROOT = HERE.parents[3] / "oph-physics-sim"
RUNS = RER_ROOT / "code" / "particles" / "runs" / "flavor"
PAYLOAD_OUT = RUNS / "family_transport_kernel_sim_s3_computed.json"
VERDICT_OUT = RUNS / "family_transport_kernel_sim_s3_shape_verdict.json"
SCAN_PATH = RUNS / "quark_common_scheme_shape_law_scan.json"

LEVEL0_RUN = "oph_universe_4k_viewer_smoke_20260706_v2"
LEVEL1_RUN = "oph_universe_64k_final_audited_20260711"
CLASS_FIELDS = tuple(f"s3_sector_class:class_delta_{c}" for c in range(3))

sys.path.insert(0, str(HERE))
from derive_generation_bundle_branch_generator import (  # noqa: E402
    build_artifact as build_generator,
)
from derive_quark_kernel_normalization_acceptance_harness import (  # noqa: E402
    evaluate as harness_evaluate,
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _encode(matrix: np.ndarray) -> dict:
    return {"real": np.real(matrix).tolist(), "imag": np.imag(matrix).tolist()}


def load_run(run_id: str) -> dict:
    run_dir = SIM_ROOT / "runs" / run_id
    cache = json.loads((run_dir / "modular_response_kernel_cache.json")
                       .read_text(encoding="utf-8"))
    payload = np.load(run_dir / "modular_response_kernel_payload.npz")
    columns: dict[tuple[int, int], dict[int, int]] = {}
    for row in cache["feature_rows"]:
        field = row["field"]
        if field in CLASS_FIELDS:
            cls = CLASS_FIELDS.index(field)
            key = (int(row["cap_index"]), int(row["time_index"]))
            columns.setdefault(key, {})[cls] = int(row["feature_index"])
    blocks = sorted(columns)
    for key in blocks:
        if sorted(columns[key]) != [0, 1, 2]:
            raise AssertionError(f"incomplete class triple at {key} in {run_id}")
    return {"run_id": run_id, "cache": cache, "payload": payload,
            "blocks": blocks, "columns": columns}


def sector_descendant(run: dict, array_key: str) -> np.ndarray:
    matrix = np.asarray(run["payload"][array_key], dtype=float)
    stacked = []
    for key in run["blocks"]:
        cols = [run["columns"][key][cls] for cls in range(3)]
        stacked.append(matrix[:, cols])
    t = np.vstack(stacked)
    d = t.T @ t
    if not np.allclose(d, d.T, atol=1.0e-10):
        raise AssertionError("descendant is not symmetric")
    return d


def readout(d: np.ndarray) -> dict:
    centered = d - (np.trace(d) / 3.0) * np.eye(3)
    evals, evecs = np.linalg.eigh(centered)
    g21 = float(evals[1] - evals[0])
    g32 = float(evals[2] - evals[1])
    span = float(evals[2] - evals[0])
    r = g21 / g32 if g32 > 0 else float("inf")
    democratic = np.ones(3) / np.sqrt(3.0)
    overlaps = [float(abs(democratic @ evecs[:, i]) ** 2) for i in range(3)]
    return {
        "eigenvalues_centered": [float(x) for x in evals],
        "g21": g21,
        "g32": g32,
        "span_cache_units": span,
        "raw_gap_ratio_r": r,
        "rho_ord": 3.0 * g32 / (2.0 * g32 + g21) if (2.0 * g32 + g21) > 0
                   else float("nan"),
        "x2": (r - 1.0) / (r + 1.0),
        "democratic_direction_overlap_by_eigenvector": overlaps,
    }


def build(level0_id: str, level1_id: str) -> tuple[dict, dict]:
    level0 = load_run(level0_id)
    level1 = load_run(level1_id)

    d0 = sector_descendant(level0, "matrix")
    d1 = sector_descendant(level1, "matrix")
    read0, read1 = readout(d0), readout(d1)

    controls = {}
    for key in level1["payload"].files:
        if key == "matrix":
            continue
        controls[key] = readout(sector_descendant(level1, key))

    scan = json.loads(SCAN_PATH.read_text(encoding="utf-8"))
    band = scan["findings"]["a2_band_by_convention"]["quoted_mixed"]
    template_rho = float(scan["template_rho_ord"])

    payload = {
        "artifact": "oph_family_transport_kernel",
        "generated_utc": _timestamp(),
        "status": "computed_sim_s3_candidate",
        "transport_kind": "s3_sector_response_second_moment_compression",
        "proof_status": "computed_candidate_carrier_identification_open",
        "carrier_identification": {
            "claim": "the three S3 conjugacy-class response channels are "
                     "the same-label generation carrier",
            "status": "candidate_claim_not_theorem",
            "uniqueness_note": "the unique exactly-three-channel field "
                               "family in the modular response cache",
        },
        "refinements": [
            {"level": 0, "run_id": level0_id,
             "observers": len(level0["cache"]["observer_ids"]),
             "block_count": len(level0["blocks"]),
             "hermitian_descendant": _encode(d0.astype(complex))},
            {"level": 1, "run_id": level1_id,
             "observers": len(level1["cache"]["observer_ids"]),
             "block_count": len(level1["blocks"]),
             "hermitian_descendant": _encode(d1.astype(complex))},
        ],
        "refinement_intertwiners": [
            {"from_level": 0, "to_level": 1,
             "intertwiner": _encode(np.eye(3, dtype=complex)),
             "justification": "canonical class frame: identity, "
                              "transposition, threecycle are absolute "
                              "labels across runs"},
        ],
        "ancestry": {
            "artifacts": [
                f"modular_response_kernel_cache({level0_id})",
                f"modular_response_kernel_cache({level1_id})",
            ],
            "attestations": {
                "quark_reference_values_consumed": False,
                "fitted_spreads_consumed": False,
                "numerical_flavor_template_consumed": False,
            },
        },
        "metadata": {
            "note": "First computed kernel of the issue-377 K1 program: "
                    "second-moment compression of the realized S3 "
                    "sector-response channels at two patch counts. The "
                    "carrier identification and the K2 scale clause stay "
                    "open; shape invariants are the milestone.",
        },
    }

    generator = build_generator(payload)
    candidate = {
        "candidate_id": "sim_s3_sector_second_moment_kernel",
        "kernel": {"refinements": payload["refinements"]},
        "ancestry": payload["ancestry"],
    }
    gates = harness_evaluate(candidate)

    verdict = {
        "artifact": "oph_family_transport_kernel_sim_s3_shape_verdict",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "row_class": "computed_candidate_shape_milestone",
        "guards": {
            "quark_reference_values_consumed": False,
            "fitted_spreads_consumed": False,
            "numerical_flavor_template_consumed": False,
            "public_promotion_allowed": False,
            "post_hoc_construction_changes": "none; the construction was "
                                             "predeclared in the module "
                                             "docstring before computation",
        },
        "levels": {
            "level0": {"run_id": level0_id, **read0},
            "level1": {"run_id": level1_id, **read1},
        },
        "controls_level1": controls,
        "band_check": {
            "declared_convention_band": band,
            "template_rho_ord": template_rho,
            "computed_rho_ord_level1": read1["rho_ord"],
            "inside_band": bool(band[0] <= read1["rho_ord"] <= band[1]),
            "physical_minus_template": read1["rho_ord"] - template_rho,
        },
        "harness_gates": gates,
        "generator_proof_status": generator.get("proof_status"),
        "reading_rules": [
            "inside_band alone is weak evidence: by the sensitivity audit "
            "about a third of generic template-class kernels land inside; "
            "the physical-vs-null contrast and level stability carry the "
            "weight",
            "G4 blocking is the expected fail-closed outcome: no K2 "
            "normalization clause exists, so spans stay unemitted",
        ],
    }
    return payload, verdict, generator


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute the family-transport kernel from simulated "
                    "S3 sector responses.")
    parser.add_argument("--level0-run", default=LEVEL0_RUN)
    parser.add_argument("--level1-run", default=LEVEL1_RUN)
    parser.add_argument("--payload-out", default=str(PAYLOAD_OUT))
    parser.add_argument("--verdict-out", default=str(VERDICT_OUT))
    args = parser.parse_args()

    payload, verdict, generator = build(args.level0_run, args.level1_run)
    for path, blob in ((pathlib.Path(args.payload_out), payload),
                       (pathlib.Path(args.verdict_out), verdict)):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(blob, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    gen_path = RUNS / "generation_bundle_branch_generator_sim_s3_computed.json"
    gen_path.write_text(json.dumps(generator, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    for label in ("level0", "level1"):
        row = verdict["levels"][label]
        print(f"{label} ({row['run_id']}):")
        print(f"  centered spectrum {['%.5g' % v for v in row['eigenvalues_centered']]}")
        print(f"  r = {row['raw_gap_ratio_r']:.6f}  rho_ord = {row['rho_ord']:.6f}  "
              f"x2 = {row['x2']:+.6f}")
        print(f"  democratic-direction overlaps {['%.3f' % v for v in row['democratic_direction_overlap_by_eigenvector']]}")
    bc = verdict["band_check"]
    print(f"band {bc['declared_convention_band']}  template {bc['template_rho_ord']:.6f}")
    print(f"COMPUTED rho_ord = {bc['computed_rho_ord_level1']:.6f}  "
          f"inside_band = {bc['inside_band']}  "
          f"delta_vs_template = {bc['physical_minus_template']:+.6f}")
    print("controls (level 1):")
    for name, row in verdict["controls_level1"].items():
        print(f"  {name:34s} rho_ord = {row['rho_ord']:.6f}  r = {row['raw_gap_ratio_r']:.6f}")
    print(f"harness: first_failed_gate = {verdict['harness_gates']['first_failed_gate']}  "
          f"status = {verdict['harness_gates']['status']}")
    print(f"saved: {args.payload_out}")
    print(f"saved: {args.verdict_out}")
    print(f"saved: {gen_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
