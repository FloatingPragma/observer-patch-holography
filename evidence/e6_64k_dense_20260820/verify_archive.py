#!/usr/bin/env python3
"""Fail-closed verifier for the curated e6 65k OPH-FPE run archive."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
CONTROL_FILES = {"README.md", "archive_manifest.json", "verify_archive.py"}


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load_json(name: str) -> Any:
    with (ROOT / name).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_inventory(archive: dict[str, Any]) -> None:
    expected = {row["path"]: row for row in archive["inventory"]}
    actual = {
        path.name
        for path in ROOT.iterdir()
        if path.is_file() and path.name not in CONTROL_FILES
    }
    require(actual == set(expected), f"inventory names differ: expected={set(expected) - actual}, extra={actual - set(expected)}")

    inventory_lines: list[str] = []
    total_bytes = 0
    for name in sorted(expected):
        row = expected[name]
        path = ROOT / name
        size = path.stat().st_size
        digest = sha256_file(path)
        require(size == row["bytes"], f"byte count mismatch: {name}")
        require(digest == row["sha256"], f"SHA-256 mismatch: {name}")
        total_bytes += size
        inventory_lines.append(f"{digest}  {size}  {name}\n")

    curated = archive["curated_archive"]
    require(len(expected) == curated["file_count"] == 45, "curated file count mismatch")
    require(total_bytes == curated["total_bytes"] == 66309382, "curated byte count mismatch")
    tree_digest = hashlib.sha256("".join(inventory_lines).encode("utf-8")).hexdigest()
    require(tree_digest == curated["inventory_sha256"], "curated inventory digest mismatch")


def verify_observer_stream(archive: dict[str, Any]) -> None:
    expected = archive["observer_views"]
    zstd = shutil.which("zstd")
    require(zstd is not None, "zstd executable is required to verify observer_views.jsonl.zst")

    process = subprocess.Popen(
        [zstd, "-q", "-d", "-c", str(ROOT / expected["archived_path"])],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    require(process.stdout is not None and process.stderr is not None, "failed to open zstd pipes")
    digest = hashlib.sha256()
    byte_count = 0
    row_count = 0
    view_types: Counter[str] = Counter()
    try:
        for line in process.stdout:
            digest.update(line)
            byte_count += len(line)
            row_count += 1
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise VerificationError(f"invalid observer JSONL row {row_count}: {exc}") from exc
            view_types[str(row.get("view_type"))] += 1
    finally:
        process.stdout.close()
    stderr = process.stderr.read().decode("utf-8", errors="replace")
    process.stderr.close()
    return_code = process.wait()
    require(return_code == 0, f"zstd decompression failed: {stderr.strip()}")
    require(byte_count == expected["uncompressed_bytes"] == 591815800, "observer stream byte count mismatch")
    require(digest.hexdigest() == expected["uncompressed_sha256"], "observer stream SHA-256 mismatch")
    require(row_count == expected["jsonl_row_count"] == 2052, "observer stream row count mismatch")
    require(view_types == Counter({"patch_observer": 2048, "cap_observer": 4}), f"observer row types mismatch: {view_types}")


def verify_run_claims(archive: dict[str, Any]) -> None:
    expected = archive["run_snapshot"]
    original = archive["source"]["original_run_tree"]
    require(original["file_count"] == 281, "original-run file-count metadata mismatch")
    require(original["total_bytes"] == 1300918262, "original-run byte-count metadata mismatch")
    require(original["inventory_sha256"] == "8c71ec07c6f549ffa5c9cd0a7de732c0097f9e0cb0bfdfe044cec9fd8d9541f1", "original-run informational digest metadata mismatch")
    require(original["verifier_enforced"] is False, "original-run digest must remain explicitly informational")
    require(archive["source"]["original_stdout_log_retained"] is False, "stdout-log custody boundary changed")
    require(archive["curated_archive"]["is_complete_copy_of_original_run"] is False, "curated archive must not claim to be the complete source tree")
    manifest = load_json("manifest.json")
    gauge = load_json("gauge_coupled_dynamics_report.json")
    observers = load_json("observer_consensus_report.json")
    population = load_json("observer_population_report.json")
    replay = load_json("finite_consensus_replay_report.json")
    theorem_core = load_json("theorem_core_receipts.json")
    geometry = load_json("icosahedral_federation_geometry_report.json")
    patch_state = load_json("echosahedral_patch_state_report.json")
    source_observer = load_json("source_dynamics_repair_record_observer_report.json")
    auto_summary = load_json("AUTO_THEOREM_UNIVERSE_SUMMARY.json")
    emergence = load_json("emergence_status_report.json")
    physical = load_json("physical_cmb_promotion_audit_report.json")
    readiness = load_json("large_run_readiness_report.json")

    require(manifest["run_id"] == archive["archive_id"] == "e6_64k_dense_20260820", "run identity mismatch")
    require(manifest["name"] == "e6_axiom_manifest_64k_dense_observers", "configuration name mismatch")
    require(manifest["patch_count"] == expected["patch_count"] == 65536, "patch count mismatch")
    require(manifest["cycles"] == expected["cycles"] == 128, "cycle count mismatch")
    require(manifest["edge_count"] == gauge["edge_count"] == expected["edge_count"] == 390924, "edge count mismatch")
    require(patch_state["ports_per_patch"] == expected["ports_per_patch"] == 12, "ports-per-patch mismatch")
    require(gauge["initial_covariant_mismatch_count"] == expected["initial_covariant_mismatch_count"] == 326047, "initial covariant mismatch mismatch")
    require(gauge["final_covariant_mismatch_count"] == 0, "final covariant mismatch is not zero")
    require(observers["observer_count"] == observers["analyzed_observer_count"] == expected["observer_sample_count"] == 2048, "observer count mismatch")
    require(observers["neighborhood_size"] == expected["observer_neighborhood_size"] == 96, "observer neighborhood mismatch")
    require(population["verbose_jsonl_patch_observer_count"] == 2048, "patch-observer materialization mismatch")
    require(population["verbose_jsonl_cap_observer_count"] == 4, "cap-observer materialization mismatch")

    with (ROOT / "mismatch_trace.csv").open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == 128, "mismatch trace must contain 128 cycle rows")
    require([int(row["cycle"]) for row in rows] == list(range(128)), "mismatch trace cycles are not contiguous 0..127")
    by_cycle = {int(row["cycle"]): row for row in rows}
    require(int(by_cycle[50]["phi"]) == expected["cycle_50_mismatch_count"] == 178280, "cycle-50 mismatch mismatch")
    zero_cycles = [int(row["cycle"]) for row in rows if int(row["phi"]) == 0]
    require(min(zero_cycles) == expected["first_zero_mismatch_cycle"] == 96, "first zero-mismatch cycle mismatch")
    require(int(by_cycle[100]["committed_records"]) == expected["cycle_100_committed_records"] == 28888, "cycle-100 commit mismatch")
    full_commit_cycles = [
        int(row["cycle"]) for row in rows if int(row["committed_records"]) == 65536
    ]
    require(
        min(full_commit_cycles) == expected["first_all_committed_cycle"] == 107,
        "first all-committed cycle mismatch",
    )
    require(int(by_cycle[110]["committed_records"]) == expected["cycle_110_committed_records"] == 65536, "cycle-110 commit mismatch")
    require(all(row["class_conformance_ok"] == "True" for row in rows), "trace contains class-conformance failure")
    require(all(row["sampled_transaction_ok"] == "True" for row in rows), "trace contains sampled-transaction failure")

    evidence = replay["evidence"]
    require(evidence["accepted_theorem_move_count"] == expected["accepted_theorem_move_count"] == 326047, "accepted move count mismatch")
    require(evidence["strict_descent_violation_count"] == expected["strict_descent_violation_count"] == 0, "strict-descent violation mismatch")
    require(replay["local_diamond_checked_pair_count"] == expected["local_diamond_checked_pair_count"] == 512, "local-diamond count mismatch")
    require(replay["shared_node_diamond_checked_pair_count"] == expected["shared_node_diamond_checked_pair_count"] == 256, "shared-node diamond count mismatch")
    require(replay["disjoint_commutation_checked_pair_count"] == expected["disjoint_commutation_checked_pair_count"] == 256, "disjoint-pair count mismatch")
    require(replay["gauge_relabeling_check_count"] == expected["gauge_relabeling_check_count"] == 4, "gauge-relabeling check count mismatch")
    require(replay["gauge_covariance_violation_count"] == expected["gauge_covariance_violation_count"] == 0, "gauge-covariance violation mismatch")
    lyapunov = theorem_core["lyapunov"]
    require(
        lyapunov["cross_cycle_injection_count"]
        == expected["cross_cycle_injection_count"]
        == 69,
        "cross-cycle injection count mismatch",
    )
    require(
        lyapunov["LYAPUNOV_DESCENT_RECEIPT"]
        is expected["global_lyapunov_descent_receipt"]
        is False,
        "global Lyapunov receipt must remain false",
    )

    endpoint = replay["exact_endpoint_branch_check"]
    require(replay["finite_consensus_theorem_receipt"] is False and replay["FINITE_CONSENSUS_THEOREM_RECEIPT"] is False, "finite-consensus receipt must remain false")
    require(endpoint["coverage_complete"] is True, "endpoint structural check is incomplete")
    require(endpoint["structurally_confluent"] is expected["structurally_confluent"] is False, "structural nonconfluence status mismatch")
    require(endpoint["structural_nonconfluence_witness_count"] >= 1 and endpoint["witness"], "structural nonconfluence witness missing")
    require(endpoint["unique_terminal_quotient_hash_count"] == expected["unique_terminal_quotient_hash_count"] == 2, "terminal quotient count mismatch")
    require(endpoint["terminal_orbit_count_semantics"] == expected["terminal_orbit_count_semantics"] == "exact lower bound", "terminal-orbit semantics mismatch")
    require(endpoint["terminal_quotient_hashes"] == [] and expected["exact_endpoint_terminal_hashes_materialized"] == 0, "exact endpoint hash materialization boundary mismatch")
    require(evidence["sampled_unique_terminal_quotient_hash_count"] == expected["sampled_unique_terminal_quotient_hash_count"] == 1, "sampled terminal quotient count mismatch")
    require(len(replay["sampled_terminal_hashes"]) == 1, "sampled terminal hash materialization mismatch")
    witness = endpoint["witness"]
    require(witness["divergent_branches"] == ["repair_left_endpoint", "repair_right_endpoint"], "endpoint witness branches changed")
    require(witness["active_edge_endpoint"] == "left" and witness["node_degree"] == 12, "endpoint witness incidence fields changed")

    require(geometry["geometry_family"] == archive["known_boundaries"]["geometry_family"] == "legacy_fibonacci_knn_control", "legacy geometry label mismatch")
    require(geometry["TRUE_ICOSAHEDRAL_REFINEMENT_TOWER_RECEIPT"] is False, "icosahedral refinement receipt unexpectedly true")
    require(geometry["CARRIER_TO_SUPPORT_REALIZATION_RECEIPT"] is False, "carrier-to-support receipt unexpectedly true")
    require(patch_state["artifact"]["written"] is False and patch_state["artifact"]["path"] is None, "patch-state artifact boundary mismatch")
    require(patch_state["artifact"]["reason"] == "compact_output_profile", "patch-state omission reason mismatch")
    require(patch_state["CARRIER_QUOTIENT_INVARIANCE_RECEIPT"] is False, "carrier quotient receipt unexpectedly true")
    require(patch_state["CARRIER_REFINEMENT_NATURALITY_RECEIPT"] is False, "carrier refinement receipt unexpectedly true")
    require(patch_state["PHYSICAL_STANDARD_MODEL_EMERGENCE_RECEIPT"] is False, "physical Standard Model receipt unexpectedly true")

    boundaries = archive["known_boundaries"]
    dedicated_self_reading = source_observer[
        "OBSERVER_LIKE_SELF_READING_SYSTEM_RECEIPT"
    ]
    aggregate_self_reading = auto_summary["final_receipts"][
        "observer_like_self_reading_system_receipt"
    ]
    require(
        dedicated_self_reading
        is boundaries["dedicated_observer_like_self_reading_system_receipt"]
        is False,
        "dedicated observer-like self-reading receipt must remain false",
    )
    require(
        source_observer["RECORD_READ_AFTER_WRITE_RECEIPT"] is False
        and source_observer["OBSERVER_READBACK_FEEDBACK_CAUSAL_LOOP_RECEIPT"]
        is False,
        "dedicated readback/feedback blockers changed",
    )
    require(
        aggregate_self_reading
        is boundaries[
            "aggregate_auto_summary_observer_like_self_reading_system_receipt"
        ]
        is True,
        "historical aggregate observer-like receipt value changed",
    )
    require(
        boundaries["aggregate_observer_like_receipt_semantic_collision_documented"]
        is True,
        "aggregate receipt disagreement is not documented",
    )

    simulator = manifest["source_provenance"]["simulator"]
    research = manifest["source_provenance"]["research"]
    require(simulator["commit"] == archive["source"]["commit"] == "b52196b296435d704b14d005d1f69caaaa662f97", "simulator commit mismatch")
    require(simulator["dirty"] is False and simulator["untracked_file_count"] == 0, "simulator checkout was not clean")
    require(research["commit"] == archive["source"]["research_checkout"]["commit"], "research commit mismatch")
    require(research["dirty"] is True, "dirty research provenance was lost")
    require(research["worktree_state_sha256"].removeprefix("sha256:") == archive["source"]["research_checkout"]["worktree_state_sha256"], "research worktree digest mismatch")

    require(emergence["status"] == "diagnostic_only_state_derived_controls_failed", "emergence status changed")
    false_emergence_gates = (
        "FINITE_CONSENSUS_THEOREM_RECEIPT",
        "bulk_3d_established",
        "ENDOGENOUS_MODULAR_GENERATOR_RECEIPT",
        "KMS_GEOMETRIC_CLOCK_FIT_RECEIPT",
        "OBJECT_BULK_POPULATION_RECEIPT",
        "EVENT_MANIFOLD_3P1D_RECEIPT",
        "DYNAMIC_DARK_TRANSPORT_RECEIPT",
        "COSMOLOGY_PERTURBATION_RECEIPT",
        "physical_claim",
    )
    require(all(emergence[key] is False for key in false_emergence_gates), "one or more fail-closed emergence gates became true")
    require(physical["physical_cmb_prediction"] is False, "physical CMB prediction gate unexpectedly true")
    require(physical["physical_cmb_promotion_ready"] is False, "physical CMB promotion gate unexpectedly true")
    require(physical["official_likelihood_ready"] is False, "official likelihood gate unexpectedly true")
    require(len(physical["contract_blockers"]) > 0 and len(physical["promotion_blockers"]) > 0, "physical promotion blockers missing")
    require(readiness["cloud_run_safe_for_physical_cmb_prediction"] is False, "physical CMB cloud-readiness gate unexpectedly true")
    require(readiness["cloud_run_safe_for_strict_neutral_bulk_claim"] is False, "strict-neutral cloud-readiness gate unexpectedly true")


def main() -> int:
    try:
        archive = load_json("archive_manifest.json")
        require(archive["schema"] == "oph.curated_run_archive.v1", "unsupported archive schema")
        verify_inventory(archive)
        verify_run_claims(archive)
        verify_observer_stream(archive)
    except (OSError, KeyError, TypeError, ValueError, VerificationError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print(
        "PASS: archive inventory, observer stream, finite diagnostics, provenance, "
        "fail-closed gates, reported witness fields, and documented source/aggregate "
        "receipt disagreement verified (no independent kernel replay)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
