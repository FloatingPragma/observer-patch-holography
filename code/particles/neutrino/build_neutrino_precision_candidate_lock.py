#!/usr/bin/env python3
"""Freeze the weighted-cycle point as a conditional, non-promotable candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import subprocess
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "precision_candidate_lock" / "v1.0.0"

ARTIFACT_INPUTS = {
    "family_transport_kernel": ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json",
    "overlap_edge_line_lift": ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json",
    "overlap_edge_transport_cocycle": ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json",
    "same_label_readback": ROOT / "particles" / "runs" / "neutrino" / "realized_same_label_gap_defect_readback.json",
    "same_label_scalar_certificate": ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json",
    "intrinsic_phase_bundle": ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_mass_eigenstate_bundle_from_scalar_certificate.json",
    "transport_load_selector": ROOT / "particles" / "runs" / "neutrino" / "neutrino_transport_load_segment_selector.json",
    "weighted_cycle_candidate": ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json",
    "correlated_profile_score": ROOT / "particles" / "runs" / "neutrino" / "nufit61_weighted_cycle_retrospective_score.json",
    "pmns_convention_audit": ROOT / "particles" / "runs" / "neutrino" / "nufit61_pmns_convention_audit.json",
    "lane_closure_contract": ROOT / "particles" / "runs" / "neutrino" / "neutrino_lane_closure_contract.json",
}

CODE_INPUTS = {
    "family_kernel_builder": ROOT / "particles" / "flavor" / "derive_family_transport_kernel.py",
    "cocycle_builder": ROOT / "particles" / "flavor" / "derive_overlap_edge_transport_cocycle.py",
    "weighted_cycle_builder": ROOT / "particles" / "neutrino" / "derive_neutrino_weighted_cycle_repair.py",
    "profile_scorer": ROOT / "particles" / "neutrino" / "score_neutrino_nufit61.py",
    "convention_auditor": ROOT / "particles" / "neutrino" / "audit_neutrino_pmns_conventions.py",
    "lock_builder": pathlib.Path(__file__).resolve(),
}

EDGES = [
    ["family_transport_kernel", "overlap_edge_line_lift"],
    ["family_transport_kernel", "overlap_edge_transport_cocycle"],
    ["overlap_edge_line_lift", "overlap_edge_transport_cocycle"],
    ["family_transport_kernel", "same_label_readback"],
    ["overlap_edge_line_lift", "same_label_readback"],
    ["same_label_readback", "same_label_scalar_certificate"],
    ["same_label_scalar_certificate", "intrinsic_phase_bundle"],
    ["overlap_edge_transport_cocycle", "transport_load_selector"],
    ["same_label_scalar_certificate", "weighted_cycle_candidate"],
    ["intrinsic_phase_bundle", "weighted_cycle_candidate"],
    ["transport_load_selector", "weighted_cycle_candidate"],
    ["overlap_edge_transport_cocycle", "weighted_cycle_candidate"],
    ["weighted_cycle_candidate", "correlated_profile_score"],
    ["weighted_cycle_candidate", "pmns_convention_audit"],
    ["correlated_profile_score", "pmns_convention_audit"],
    ["weighted_cycle_candidate", "lane_closure_contract"],
]


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative(path: pathlib.Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def _artifact_node(node_id: str, path: pathlib.Path) -> dict[str, Any]:
    payload = _load_json(path)
    return {
        "id": node_id,
        "kind": "artifact",
        "path": _relative(path),
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
        "artifact": payload.get("artifact"),
        "status": payload.get("status"),
        "proof_status": payload.get("proof_status"),
        "prediction_promotion_allowed": payload.get("prediction_promotion_allowed"),
        "public_surface_candidate_allowed": payload.get("public_surface_candidate_allowed"),
    }


def _code_node(node_id: str, path: pathlib.Path) -> dict[str, Any]:
    return {
        "id": node_id,
        "kind": "code",
        "path": _relative(path),
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _dag_is_acyclic(node_ids: set[str], edges: list[list[str]]) -> bool:
    incoming = {node_id: 0 for node_id in node_ids}
    outgoing = {node_id: [] for node_id in node_ids}
    for source, target in edges:
        if source not in node_ids or target not in node_ids:
            raise ValueError(f"DAG edge refers to an unknown node: {source} -> {target}")
        incoming[target] += 1
        outgoing[source].append(target)
    ready = [node_id for node_id, count in incoming.items() if count == 0]
    visited = 0
    while ready:
        source = ready.pop()
        visited += 1
        for target in outgoing[source]:
            incoming[target] -= 1
            if incoming[target] == 0:
                ready.append(target)
    return visited == len(node_ids)


def _git_head() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _relevant_inputs_clean(paths: list[pathlib.Path]) -> bool:
    result = subprocess.run(
        ["git", "status", "--porcelain", "--", *[str(path) for path in paths]],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return result.returncode == 0 and not result.stdout.strip()


def _write_json(path: pathlib.Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_bundle(output: pathlib.Path, lock_time: str) -> dict[str, Any]:
    if output.exists():
        raise FileExistsError(f"candidate lock is immutable and already exists: {output}")
    output.mkdir(parents=True)

    artifacts = {name: _load_json(path) for name, path in ARTIFACT_INPUTS.items()}
    candidate = artifacts["weighted_cycle_candidate"]
    score = artifacts["correlated_profile_score"]
    convention = artifacts["pmns_convention_audit"]
    kernel = artifacts["family_transport_kernel"]
    line_lift = artifacts["overlap_edge_line_lift"]
    cocycle = artifacts["overlap_edge_transport_cocycle"]

    nodes = [_artifact_node(name, path) for name, path in ARTIFACT_INPUTS.items()]
    nodes.extend(_code_node(name, path) for name, path in CODE_INPUTS.items())
    artifact_node_ids = set(ARTIFACT_INPUTS)
    if not _dag_is_acyclic(artifact_node_ids, EDGES):
        raise ValueError("source DAG contains a cycle")

    pmns = candidate["pmns_observables"]
    ratio_21_32 = float(candidate["dimensionless_ratio_dm21_over_dm32"])
    sin2_theta12 = math.sin(float(pmns["theta12_rad"])) ** 2
    sin2_theta13 = math.sin(float(pmns["theta13_rad"])) ** 2
    sin2_theta23 = math.sin(float(pmns["theta23_rad"])) ** 2
    ratio_21_31 = ratio_21_32 / (1.0 + ratio_21_32)
    cos2_theta12 = 1.0 - sin2_theta12
    ratio_21_ee = ratio_21_32 / (1.0 + cos2_theta12 * ratio_21_32)

    source_failures = []
    if kernel.get("status") == "template":
        source_failures.append("family transport kernel is a hand-written template")
    if line_lift.get("proof_status") != "closed_source_emitted":
        source_failures.append(f"line-lift status is {line_lift.get('proof_status')!r}")
    if cocycle.get("proof_status") != "closed_source_emitted":
        source_failures.append(f"cocycle status is {cocycle.get('proof_status')!r}")
    if candidate.get("historical_target_exposure") is not False:
        source_failures.append("selector-law development used oscillation targets")

    source_dag = {
        "artifact": "oph_neutrino_precision_candidate_source_dag",
        "lock_id": "OPH-NU-CAND-001",
        "version": "1.0.0",
        "lock_time_utc": lock_time,
        "acyclic": True,
        "nodes": nodes,
        "edges": EDGES,
        "frozen_human_choices": {
            "cycle_basis_order": candidate["cycle_basis_order"],
            "holonomy_orientation": candidate["holonomy_orientation"],
            "pmns_row_order": candidate["pmns_row_order_for_pdg"],
            "transport_load_selector": candidate["selected_transport_load_selector"],
            "weight_exponent_formula": candidate["weight_exponent_formula"],
        },
        "target_audit": {
            "historical_target_exposure": True,
            "no_data_after_lock_time_used": True,
            "forbidden_after_lock": [
                "DUNE oscillation likelihood or event spectra released after the lock",
                "Hyper-K oscillation likelihood or event spectra released after the lock",
                "JUNO oscillation likelihood or event spectra released after the lock",
                "post-lock NuFIT or other global-fit updates",
                "any post-lock selector, orientation, row, column, or uncertainty change",
            ],
        },
        "source_closure_gate": {
            "passes": not source_failures,
            "failures": source_failures,
        },
    }

    predictions = {
        "artifact": "oph_neutrino_precision_candidate_predictions",
        "lock_id": "OPH-NU-CAND-001",
        "version": "1.0.0",
        "claim_class": "retrospective_target_informed_template_candidate",
        "prediction_promotion_allowed": False,
        "ordering": "normal",
        "ordering_status": "locked_declared_hypothesis_not_source_derived",
        "independent_coordinates": {
            "sin2_theta12": sin2_theta12,
            "sin2_theta13": sin2_theta13,
            "sin2_theta23": sin2_theta23,
            "delta_cp_rad": float(pmns["delta_rad"]),
            "delta_cp_deg_wrapped": ((float(pmns["delta_deg"]) + 180.0) % 360.0) - 180.0,
            "ratio_delta_m21_sq_over_delta_m32_sq": ratio_21_32,
        },
        "derived_crosschecks": {
            "J": float(pmns["J"]),
            "ratio_delta_m21_sq_over_delta_m31_sq": ratio_21_31,
            "ratio_delta_m21_sq_over_delta_mee_sq": ratio_21_ee,
        },
        "absolute_scale": {
            "status": "one_positive_profiled_nuisance",
            "absolute_mass_prediction_allowed": False,
        },
        "experiment_views": {
            "DUNE": ["ordering", "sin2_theta23", "delta_cp", "sin2_theta13"],
            "Hyper-K": ["ordering", "sin2_theta23", "delta_cp", "sin2_theta13"],
            "JUNO": ["ordering", "sin2_theta12", "sin2_theta13", "ratio_delta_m21_sq_over_delta_mee_sq"],
        },
        "theory_epistemic_distribution": None,
        "numerical_policy": "exact float64 point; no experimental error is re-used as theory uncertainty",
    }

    score_rows = {
        source_id: {
            "T23_DCP_delta_chi2": result["profiles"]["T23/DCP"]["delta_chi2"],
            "lower_bound_delta_chi2": result["joint_fixed_candidate_delta_chi2_lower_bound"],
            "passes_3sigma_2d": result["passes_published_3sigma_2d_compatibility"],
        }
        for source_id, result in score["scores"].items()
    }
    adjudication = {
        "artifact": "oph_neutrino_precision_candidate_adjudication",
        "lock_id": "OPH-NU-CAND-001",
        "version": "1.0.0",
        "baseline": {
            "dataset": "NuFIT 6.1 official marginalized profile tables",
            "retrospective": True,
            "profiles_are_never_summed": True,
            "scores": score_rows,
            "candidate_rejected_by_declared_3sigma_gate": score["decision"]["current_weighted_cycle_candidate_rejected_by_declared_gate"],
            "stored_pmns_relabeling_rescue_found": convention["decision"].get(
                "stored_pmns_relabeling_rescue_found",
                convention["decision"].get("convention_error_found", False),
            ),
            "physical_basis_contract_audited": convention.get("scope", {}).get("weighted_cycle_operator_basis_placement_audited", False),
        },
        "future_no_rescue_rule": {
            "point_may_change_after_lock": False,
            "ordering_may_change_after_lock": False,
            "orientation_or_label_mapping_may_change_after_lock": False,
            "theory_uncertainty_may_be_added_after_lock": False,
        },
        "future_test_statistic": "official joint likelihood ratio against the free three-neutrino model, profiling one common positive oscillation scale and shared experimental nuisances",
        "future_decision_levels": {
            "disfavored": "calibrated p < 0.05",
            "strongly_disfavored": "calibrated p < 0.0027",
            "falsified_candidate": "calibrated p < 2.87e-7 or an official joint 5-sigma exclusion",
        },
        "invalid_lock_conditions": [
            "hash mismatch",
            "post-lock target use",
            "unrecorded convention change",
            "likelihood unavailable or incompatible with the locked observable definition",
        ],
        "oph_core_falsified_by_candidate_rejection": False,
    }

    relevant_paths = list(ARTIFACT_INPUTS.values()) + list(CODE_INPUTS.values())
    repository = {
        "git_head": _git_head(),
        "relevant_inputs_clean_against_head": _relevant_inputs_clean(relevant_paths),
    }
    repository["archival_release_ready"] = repository["relevant_inputs_clean_against_head"]

    _write_json(output / "source_dag.json", source_dag)
    _write_json(output / "predictions.json", predictions)
    _write_json(output / "adjudication.json", adjudication)
    claim = """# Neutrino precision candidate lock v1.0.0

This bundle freezes the weighted-cycle point as an anti-rescue record. Its scientific status is a retrospective, target-informed template candidate. Prediction promotion is prohibited.

The official NuFIT 6.1 correlated profile rejects the frozen point under the declared two-parameter 3σ gate. Exhaustive orientation and label enumeration finds no convention correction. The result rejects this weighted-cycle candidate. It does not falsify the finite OPH core because the flavor kernel is a hand-written template and selector development used oscillation targets.

A promotable neutrino prediction requires a source-emitted flavor kernel, source-derived cycle and label laws, a pre-reference hash lock, and a new likelihood evaluation whose data were unavailable during construction.
"""
    (output / "CLAIM.md").write_text(claim, encoding="utf-8")

    content_files = ["source_dag.json", "predictions.json", "adjudication.json", "CLAIM.md"]
    content_hashes = {name: _sha256(output / name) for name in content_files}
    manifest = {
        "artifact": "oph_neutrino_precision_candidate_lock_manifest",
        "lock_id": "OPH-NU-CAND-001",
        "version": "1.0.0",
        "lock_time_utc": lock_time,
        "immutable": True,
        "prediction_promotion_allowed": False,
        "source_closure_passes": not source_failures,
        "candidate_rejected_by_nufit61": score["decision"]["current_weighted_cycle_candidate_rejected_by_declared_gate"],
        "repository": repository,
        "files": content_hashes,
    }
    _write_json(output / "manifest.json", manifest)
    checksum_lines = [f"{content_hashes[name]}  {name}" for name in sorted(content_hashes)]
    checksum_lines.append(f"{_sha256(output / 'manifest.json')}  manifest.json")
    (output / "SHA256SUMS").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the immutable neutrino precision candidate lock.")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--lock-time", required=True)
    args = parser.parse_args()
    manifest = build_bundle(pathlib.Path(args.output), args.lock_time)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
