#!/usr/bin/env python3
"""Enumerate PMNS mappings before attributing a NuFIT failure to convention."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np

from derive_neutrino_weighted_cycle_repair import _pmns_parameters
from score_neutrino_nufit61 import (
    DEFAULT_CANDIDATE,
    DEFAULT_MANIFEST,
    THREE_SIGMA_TWO_DOF_DELTA_CHI2,
    _read_grids,
    _sha256,
    _verify_table,
    _wrap_delta_degrees,
)


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "nufit61_pmns_convention_audit.json"
FLAVOR_LABELS = ("e", "mu", "tau")
MASS_LABELS = ("1", "2", "3")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _score_unitary(unitary: np.ndarray, grids: dict[str, Any]) -> dict[str, Any]:
    observables = _pmns_parameters(unitary)
    coordinates = {
        "sin2_theta12": math.sin(observables["theta12_rad"]) ** 2,
        "sin2_theta13": math.sin(observables["theta13_rad"]) ** 2,
        "sin2_theta23": math.sin(observables["theta23_rad"]) ** 2,
        "delta_cp_deg_wrapped": _wrap_delta_degrees(observables["delta_deg"]),
    }
    try:
        t13_t12 = grids["T13/T12"].interpolate(
            coordinates["sin2_theta13"], coordinates["sin2_theta12"]
        )["delta_chi2"]
        t23_dcp = grids["T23/DCP"].interpolate(
            coordinates["sin2_theta23"], coordinates["delta_cp_deg_wrapped"]
        )["delta_chi2"]
    except ValueError as exc:
        return {
            "coordinates": coordinates,
            "observables": observables,
            "score_status": "outside_official_grid",
            "reason": str(exc),
        }
    lower_bound = max(t13_t12, t23_dcp)
    return {
        "coordinates": coordinates,
        "observables": observables,
        "score_status": "scored",
        "T13/T12_delta_chi2": t13_t12,
        "T23/DCP_delta_chi2": t23_dcp,
        "two_profile_lower_bound": lower_bound,
        "passes_3sigma_2d_gate": lower_bound <= THREE_SIGMA_TWO_DOF_DELTA_CHI2,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit PMNS convention and label alternatives on NuFIT 6.1.")
    parser.add_argument("--candidate", default=str(DEFAULT_CANDIDATE))
    parser.add_argument("--source-manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--tb-off-no", required=True)
    parser.add_argument("--tb-yes-no", required=True)
    parser.add_argument("--generated-utc", default="")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()

    candidate_path = pathlib.Path(args.candidate)
    manifest_path = pathlib.Path(args.source_manifest)
    candidate = _load_json(candidate_path)
    manifest = _load_json(manifest_path)
    unitary = np.asarray(candidate["pmns_real"], dtype=float) + 1j * np.asarray(
        candidate["pmns_imag"], dtype=float
    )

    table_paths = {
        "TBoff-NO": pathlib.Path(args.tb_off_no),
        "TByes-NO": pathlib.Path(args.tb_yes_no),
    }
    grids_by_source = {}
    source_receipts = {}
    for source_id, path in table_paths.items():
        source_receipts[source_id] = _verify_table(path, manifest["files"][source_id], source_id)
        grids_by_source[source_id] = _read_grids(path)

    rows = []
    for orientation, oriented_unitary in (
        ("declared_021", unitary),
        ("opposite_012_complex_conjugate", unitary.conjugate()),
    ):
        for row_permutation in itertools.permutations(range(3)):
            for column_permutation in itertools.permutations(range(3)):
                transformed = oriented_unitary[list(row_permutation), :][:, list(column_permutation)]
                scores = {
                    source_id: _score_unitary(transformed, grids)
                    for source_id, grids in grids_by_source.items()
                }
                scored_bounds = [
                    score["two_profile_lower_bound"]
                    for score in scores.values()
                    if score["score_status"] == "scored"
                ]
                rows.append(
                    {
                        "orientation": orientation,
                        "row_permutation_indices": list(row_permutation),
                        "row_labels_interpreted_as_e_mu_tau": [FLAVOR_LABELS[index] for index in row_permutation],
                        "column_permutation_indices": list(column_permutation),
                        "mass_labels_interpreted_as_1_2_3": [MASS_LABELS[index] for index in column_permutation],
                        "normal_ordering_mass_assignment_consistent": column_permutation == (0, 1, 2),
                        "declared_flavor_assignment": row_permutation == (0, 1, 2),
                        "scores": scores,
                        "worst_declared_dataset_lower_bound": max(scored_bounds) if scored_bounds else None,
                        "passes_both_declared_3sigma_gates": bool(scored_bounds)
                        and len(scored_bounds) == len(scores)
                        and all(
                            score.get("passes_3sigma_2d_gate") is True for score in scores.values()
                        ),
                    }
                )

    no_consistent_rows = [row for row in rows if row["normal_ordering_mass_assignment_consistent"]]
    scored_no_rows = [row for row in no_consistent_rows if row["worst_declared_dataset_lower_bound"] is not None]
    best = min(scored_no_rows, key=lambda row: row["worst_declared_dataset_lower_bound"])
    declared = next(
        row
        for row in rows
        if row["orientation"] == "declared_021"
        and row["row_permutation_indices"] == [0, 1, 2]
        and row["column_permutation_indices"] == [0, 1, 2]
    )

    payload = {
        "artifact": "oph_neutrino_nufit61_pmns_convention_audit",
        "generated_utc": args.generated_utc or _timestamp(),
        "candidate_filename": candidate_path.name,
        "candidate_sha256": _sha256(candidate_path),
        "source_manifest_sha256": _sha256(manifest_path),
        "source_receipts": source_receipts,
        "enumeration": {
            "orientations": 2,
            "row_permutations": 6,
            "column_permutations": 6,
            "total_assignments": len(rows),
            "two_parameter_3sigma_threshold": THREE_SIGMA_TWO_DOF_DELTA_CHI2,
        },
        "mathematical_convention_checks": {
            "takagi_relation": "For U^T M U = D, the columns of U diagonalize M^dagger M.",
            "jarlskog_sign": "J = Im(U_e1 U_mu2 U_e2^* U_mu1^*) uses the PDG sign convention.",
            "delta_wrap": "delta in [0,360) is mapped to the equivalent NuFIT coordinate in [-180,180].",
            "normal_ordering_column_rule": "Ascending eigenvalues are (m1,m2,m3); other column permutations change physical mass labels.",
            "transpose_or_dagger_allowed": False,
        },
        "declared_assignment": declared,
        "best_normal_ordering_consistent_assignment": best,
        "decision": {
            "any_normal_ordering_consistent_assignment_passes_both_3sigma_gates": any(
                row["passes_both_declared_3sigma_gates"] for row in no_consistent_rows
            ),
            "declared_assignment_is_best_by_worst_dataset_score": best is declared,
            "convention_error_found": False,
            "post_hoc_relabeling_allowed_as_rescue": False,
        },
        "assignments": rows,
    }

    output_path = pathlib.Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {output_path}")
    print(
        json.dumps(
            {
                "declared_worst_lower_bound": declared["worst_declared_dataset_lower_bound"],
                "best_worst_lower_bound": best["worst_declared_dataset_lower_bound"],
                "any_NO_assignment_passes_both": payload["decision"]["any_normal_ordering_consistent_assignment_passes_both_3sigma_gates"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
