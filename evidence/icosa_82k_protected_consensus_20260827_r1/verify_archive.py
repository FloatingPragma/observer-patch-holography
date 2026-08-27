#!/usr/bin/env python3
"""Fail-closed verifier for the protected-authority 82k run archive.

The exact normal-form replay below is deliberately implemented without
importing oph-physics-sim.  It reconstructs S3 from permutations, derives the
authority-selected terminal edge-slot state from the archived primitive
arrays, and recomputes the source, quotient, authority, and terminal hashes.
"""

from __future__ import annotations

import csv
import hashlib
from itertools import permutations
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
CONTROL_FILES = {"README.md", "archive_manifest.json", "verify_archive.py"}
SIMULATOR_COMMIT = "ce17921eb7504106fef1ba445e1349b5367aa676"
RESEARCH_COMMIT = "ee137dea6fd48487e30863db53bb787d2e5f45e8"
RUN_ID = "icosa_82k_protected_consensus_20260827_r1"
AUTHORITY_SCHEMA = "oph_protected_distinct_node_authority_v1"
PROTECTED_SOURCE_SCHEMA = "oph_protected_authority_coupled_source_hash_v1"
PROTECTED_TERMINAL_SCHEMA = "oph_protected_authority_terminal_quotient_hash_v1"
QUOTIENT_CANONICALIZER = "node_reference_port_frame_v1"


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
    require(actual == set(expected), "archive inventory names differ")
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
    require(len(expected) == curated["file_count"], "file count mismatch")
    require(total_bytes == curated["total_bytes"], "total byte count mismatch")
    digest = hashlib.sha256("".join(inventory_lines).encode("utf-8")).hexdigest()
    require(digest == curated["inventory_sha256"], "inventory digest mismatch")


def s3_tables() -> tuple[np.ndarray, np.ndarray]:
    elements = tuple(permutations((0, 1, 2)))
    index = {element: position for position, element in enumerate(elements)}

    def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(left[right[position]] for position in range(3))

    multiply = np.asarray(
        [[index[compose(left, right)] for right in elements] for left in elements],
        dtype=np.int16,
    )
    inverse = np.asarray(
        [index[tuple(element.index(position) for position in range(3))] for element in elements],
        dtype=np.int16,
    )
    return multiply, inverse


def state_hash(
    left: np.ndarray,
    right: np.ndarray,
    gauge: np.ndarray,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    *,
    schema: str,
) -> str:
    digest = hashlib.sha256()
    digest.update(
        f"oph-covariant-overlap-state-v1\0{schema}\0S3\0{6}\0".encode("ascii")
    )
    for values, dtype in (
        (left, "<i2"),
        (right, "<i2"),
        (gauge, "<i2"),
        (edge_left, "<i8"),
        (edge_right, "<i8"),
    ):
        array = np.ascontiguousarray(np.asarray(values, dtype=dtype))
        digest.update(int(array.size).to_bytes(8, "little", signed=False))
        digest.update(array.tobytes())
    return "sha256:" + digest.hexdigest()


def authority_hash(authority: np.ndarray) -> str:
    values = np.ascontiguousarray(authority, dtype="<i8")
    digest = hashlib.sha256()
    digest.update((AUTHORITY_SCHEMA + "\0").encode("ascii"))
    digest.update(int(values.size).to_bytes(8, "little", signed=False))
    digest.update(values.tobytes())
    return "sha256:" + digest.hexdigest()


def authority_bound_hash(schema: str, value_hash: str, authority: np.ndarray) -> str:
    digest = hashlib.sha256()
    digest.update((schema + "\0").encode("ascii"))
    digest.update(value_hash.encode("ascii"))
    digest.update(b"\0")
    digest.update(authority_hash(authority).encode("ascii"))
    return "sha256:" + digest.hexdigest()


def canonicalize_quotient(
    left: np.ndarray,
    right: np.ndarray,
    gauge: np.ndarray,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    multiply: np.ndarray,
    inverse: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    node_count = int(max(np.max(edge_left), np.max(edge_right))) + 1
    edge_indices = np.arange(left.size, dtype=np.int64)
    references = np.full(node_count, left.size, dtype=np.int64)
    np.minimum.at(references, edge_left, edge_indices)
    np.minimum.at(references, edge_right, edge_indices)
    frames = np.zeros(node_count, dtype=np.int16)
    left_reference = references[edge_left] == edge_indices
    right_reference = references[edge_right] == edge_indices
    frames[edge_left[left_reference]] = inverse[left[left_reference]]
    frames[edge_right[right_reference]] = inverse[right[right_reference]]
    canonical_left = multiply[frames[edge_left], left]
    canonical_right = multiply[frames[edge_right], right]
    canonical_gauge = multiply[
        multiply[frames[edge_left], gauge],
        inverse[frames[edge_right]],
    ]
    return canonical_left, canonical_right, canonical_gauge


def verify_exact_normal_form(archive: dict[str, Any]) -> None:
    replay = load_json("finite_consensus_replay_report.json")
    source_manifest = load_json("finite_consensus_source_manifest.json")
    repair_contract = load_json("repair_kernel_contract.json")
    multiply, inverse = s3_tables()
    with np.load(ROOT / "finite_consensus_source_state.npz", allow_pickle=False) as source:
        left0 = np.asarray(source["initial_port_left"], dtype=np.int16)
        right0 = np.asarray(source["initial_port_right"], dtype=np.int16)
        gauge0 = np.asarray(source["initial_gauge"], dtype=np.int16)
        edge_left = np.asarray(source["edge_left"], dtype=np.int64)
        edge_right = np.asarray(source["edge_right"], dtype=np.int64)
        authority = np.asarray(source["node_repair_authorities"], dtype=np.int64)
    with np.load(ROOT / "protected_repair_authority.npz", allow_pickle=False) as protected:
        protected_authority = np.asarray(protected["node_repair_authorities"], dtype=np.int64)
    with np.load(ROOT / "s3_gauge_state.npz", allow_pickle=False) as terminal:
        terminal_edge_left = np.asarray(terminal["left"], dtype=np.int64)
        terminal_edge_right = np.asarray(terminal["right"], dtype=np.int64)
        final_gauge = np.asarray(terminal["gauge"], dtype=np.int16)
    with np.load(ROOT / "echosahedral_patch_state.npz", allow_pickle=False) as patch:
        final_left = np.asarray(patch["routed_left_state"], dtype=np.int16)
        final_right = np.asarray(patch["routed_right_state"], dtype=np.int16)

    require(left0.shape == right0.shape == gauge0.shape == (122880,), "source edge arrays")
    require(edge_left.shape == edge_right.shape == left0.shape, "source graph arrays")
    require(authority.shape == (81920,), "authority census")
    require(np.array_equal(authority, protected_authority), "authority artifact disagreement")
    require(np.unique(authority).size == authority.size, "authorities are not distinct")
    require(set(authority.tolist()) == set(range(81920)), "authority permutation boundary")
    require(np.all((left0 >= 0) & (left0 < 6)), "left labels outside S3")
    require(np.all((right0 >= 0) & (right0 < 6)), "right labels outside S3")
    require(np.all((gauge0 >= 0) & (gauge0 < 6)), "gauge labels outside S3")
    require(np.all(edge_left != edge_right), "self-loop in source graph")
    require(np.array_equal(terminal_edge_left, edge_left), "terminal left endpoints")
    require(np.array_equal(terminal_edge_right, edge_right), "terminal right endpoints")

    mismatch = left0 != multiply[gauge0, right0]
    require(int(np.count_nonzero(mismatch)) == 102415, "initial mismatch count")
    expected_left = left0.copy()
    expected_right = right0.copy()
    repair_left = authority[edge_right] > authority[edge_left]
    write_left = mismatch & repair_left
    write_right = mismatch & ~repair_left
    expected_left[write_left] = multiply[gauge0[write_left], right0[write_left]]
    expected_right[write_right] = multiply[inverse[gauge0[write_right]], left0[write_right]]
    require(np.array_equal(final_left, expected_left), "terminal left state is not the authority normal form")
    require(np.array_equal(final_right, expected_right), "terminal right state is not the authority normal form")
    require(np.array_equal(final_gauge, gauge0), "protected gauge links changed")
    require(not np.any(final_left != multiply[final_gauge, final_right]), "terminal mismatch remains")

    source_hash = state_hash(
        left0, right0, gauge0, edge_left, edge_right, schema="coupled-representative-v1"
    )
    source_quotient = canonicalize_quotient(
        left0, right0, gauge0, edge_left, edge_right, multiply, inverse
    )
    source_quotient_hash = state_hash(
        *source_quotient,
        edge_left,
        edge_right,
        schema=f"gauge-quotient:{QUOTIENT_CANONICALIZER}",
    )
    terminal_quotient = canonicalize_quotient(
        final_left, final_right, final_gauge, edge_left, edge_right, multiply, inverse
    )
    terminal_quotient_hash = state_hash(
        *terminal_quotient,
        edge_left,
        edge_right,
        schema=f"gauge-quotient:{QUOTIENT_CANONICALIZER}",
    )
    protected_source_hash = authority_bound_hash(
        PROTECTED_SOURCE_SCHEMA, source_hash, authority
    )
    protected_terminal_hash = authority_bound_hash(
        PROTECTED_TERMINAL_SCHEMA, terminal_quotient_hash, authority
    )

    require(authority_hash(authority) == replay["authority_sha256"], "authority hash")
    require(authority_hash(authority) == source_manifest["authority_sha256"], "manifest authority hash")
    require(authority_hash(authority) == repair_contract["authority_sha256"], "contract authority hash")
    require(source_hash == replay["source_state_sha256"], "source state hash")
    require(source_quotient_hash == replay["source_quotient_hash"], "source quotient hash")
    require(protected_source_hash == replay["protected_source_sha256"], "protected source hash")
    require(protected_source_hash == source_manifest["protected_source_sha256"], "manifest protected source hash")
    require(terminal_quotient_hash == replay["terminal_quotient_hash"], "terminal quotient hash")
    require(protected_terminal_hash == replay["terminal_hash"], "protected terminal hash")

    with np.load(ROOT / "echosahedral_patch_state.npz", allow_pickle=False) as patch:
        require(np.array_equal(patch["edge_left"], edge_left), "patch left endpoints")
        require(np.array_equal(patch["edge_right"], edge_right), "patch right endpoints")
        require(bool(np.all(patch["committed"])), "not all patch records are committed")
        routed_slots = {
            (int(node), int(port))
            for node, port in zip(patch["edge_left"], patch["left_port"], strict=True)
        } | {
            (int(node), int(port))
            for node, port in zip(patch["edge_right"], patch["right_port"], strict=True)
        }
        require(len(routed_slots) == 245760, "a routed patch-port slot is shared")

    evidence = replay["evidence"]
    require(replay["receipt"] is True, "finite replay receipt is false")
    require(replay["FINITE_CONSENSUS_THEOREM_RECEIPT"] is True, "finite theorem receipt is false")
    require(evidence["accepted_theorem_move_count"] == 102415, "accepted repair count")
    require(evidence["schedule_replay_count"] == evidence["requested_schedule_replays"] == 16, "schedule replay count")
    zero_fields = (
        "strict_descent_violation_count",
        "accepted_phi_increase_violation_count",
        "disjoint_commutation_violation_count",
        "local_diamond_violation_count",
        "gauge_covariance_violation_count",
        "production_move_contract_violation_count",
        "endpoint_branch_coverage_incomplete_count",
        "endpoint_branch_confluence_violation_count",
        "repair_completeness_violation_count",
        "sector_link_mutation_count",
    )
    require(all(evidence[key] == 0 for key in zero_fields), "replay violation count")
    require(evidence["sector_replay_call_count"] == 102415 * 16, "sector replay call count")
    require(replay["local_diamond_checked_pair_count"] == 1024, "diamond check count")
    require(replay["shared_node_diamond_checked_pair_count"] == 512, "shared-node check count")
    require(replay["disjoint_commutation_checked_pair_count"] == 512, "disjoint check count")
    require(replay["gauge_relabeling_check_count"] == 16, "frame check count")
    exact = replay["exact_normalizer_confluence_check"]
    require(exact["coverage_complete"] is True, "exact normalizer coverage")
    require(exact["structurally_confluent"] is True, "exact normalizer confluence")
    require(exact["endpoint_repair_effective"] is True, "endpoint repair effectiveness")
    require(exact["endpoint_branch_nondeterministic"] is False, "nondeterministic branch remains")

    expected = archive["run_snapshot"]
    require(expected["initial_mismatch_count"] == int(np.count_nonzero(mismatch)), "snapshot mismatch count")
    require(expected["terminal_hash"] == protected_terminal_hash, "snapshot terminal hash")


def verify_reports(archive: dict[str, Any]) -> None:
    manifest = load_json("manifest.json")
    geometry = load_json("icosahedral_federation_geometry_report.json")
    gauge = load_json("gauge_coupled_dynamics_report.json")
    patch = load_json("echosahedral_patch_state_report.json")
    observer = load_json("observer_consensus_report.json")
    population = load_json("observer_population_report.json")
    source_observer = load_json("source_dynamics_repair_record_observer_report.json")
    theorem = load_json("theorem_core_receipts.json")

    require(manifest["run_id"] == archive["archive_id"] == RUN_ID, "run identity")
    require(manifest["patch_count"] == 81920, "patch count")
    require(manifest["edge_count"] == 122880, "edge count")
    require(manifest["cycles"] == 16, "cycle count")
    simulator = manifest["source_provenance"]["simulator"]
    research = manifest["source_provenance"]["research"]
    require(simulator["commit"] == SIMULATOR_COMMIT and simulator["dirty"] is False, "simulator provenance")
    require(simulator["untracked_file_count"] == 0, "simulator untracked files")
    require(research["commit"] == RESEARCH_COMMIT and research["dirty"] is False, "research provenance")
    require(research["untracked_file_count"] == 0, "research untracked files")

    require(geometry["geometry_family"] == "icosahedral_tower", "geometry family")
    require(geometry["patch_basis"] == "cells" and geometry["refinement_level"] == 6, "geometry rung")
    require(geometry["edge_count"] == 122880, "geometry edge count")
    require(geometry["TRUE_ICOSAHEDRAL_REFINEMENT_TOWER_RECEIPT"] is True, "tower receipt")
    require(geometry["CARRIER_TO_SUPPORT_REALIZATION_RECEIPT"] is False, "physical carrier boundary")
    require(gauge["initial_covariant_mismatch_count"] == 102415, "gauge initial mismatch")
    require(gauge["final_covariant_mismatch_count"] == 0, "gauge terminal mismatch")
    require(gauge["gauge_link_changed_count"] == 0, "gauge link mutation")
    require(patch["ports_per_patch"] == 12, "ports per patch")
    require(patch["routed_port_slot_count"] == 245760, "routed slot count")
    require(patch["unrouted_exposed_or_reserved_port_slot_count"] == 737280, "reserved slot count")
    require(patch["artifact"]["written"] is True, "patch-state artifact")
    require(observer["observer_count"] == observer["analyzed_observer_count"] == 2048, "observer count")
    require(observer["neighborhood_size"] == 96, "observer neighborhood")
    require(population["verbose_jsonl_patch_observer_count"] == 2048, "observer rows")
    require(population["verbose_jsonl_cap_observer_count"] == 2, "cap control rows")
    require(source_observer["RECORD_COMMIT_RECEIPT"] is True, "record commit receipt")
    require(source_observer["OBSERVER_LIKE_SELF_READING_SYSTEM_RECEIPT"] is False, "large-run feedback boundary")
    require(source_observer["OBSERVER_READBACK_FEEDBACK_CAUSAL_LOOP_RECEIPT"] is False, "large-run causal feedback boundary")
    require(source_observer["record_observer"]["readback_count"] == 0, "large-run readback count")
    require(source_observer["record_observer"]["feedback_event_count"] == 0, "large-run feedback count")
    require(theorem["finite_consensus_theorem"]["FINITE_CONSENSUS_THEOREM_RECEIPT"] is True, "theorem-core consensus")

    with (ROOT / "mismatch_trace.csv").open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == 16, "trace length")
    require([int(row["cycle"]) for row in rows] == list(range(16)), "trace cycle order")
    require(int(rows[0]["phi_before"]) == 102415 and int(rows[0]["phi"]) == 20495, "cycle-zero trace")
    require(int(rows[1]["phi_before"]) == 20495 and int(rows[1]["phi"]) == 0, "cycle-one trace")
    require(all(int(row["phi"]) == 0 for row in rows[1:]), "terminal mismatch stability")
    require(int(rows[8]["committed_records"]) == 81920, "cycle-eight commits")
    require(all(int(row["committed_records"]) == 81920 for row in rows[8:]), "commit stability")
    require(all(row["class_conformance_ok"] == "True" for row in rows), "scheduler class failure")
    require(all(row["sampled_transaction_ok"] == "True" for row in rows), "sampled transaction failure")


def verify_observer_stream(archive: dict[str, Any]) -> None:
    expected = archive["observer_views"]
    zstd = shutil.which("zstd")
    require(zstd is not None, "zstd executable is required")
    process = subprocess.Popen(
        [zstd, "-q", "-d", "-c", str(ROOT / expected["archived_path"])],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    require(process.stdout is not None and process.stderr is not None, "zstd pipes")
    digest = hashlib.sha256()
    byte_count = 0
    row_count = 0
    row_types: Counter[str] = Counter()
    for line in process.stdout:
        digest.update(line)
        byte_count += len(line)
        row_count += 1
        row_types[str(json.loads(line).get("view_type"))] += 1
    process.stdout.close()
    stderr = process.stderr.read().decode("utf-8", errors="replace")
    process.stderr.close()
    require(process.wait() == 0, f"zstd decompression failed: {stderr.strip()}")
    require(byte_count == expected["uncompressed_bytes"], "observer byte count")
    require(digest.hexdigest() == expected["uncompressed_sha256"], "observer stream hash")
    require(row_count == expected["jsonl_row_count"] == 2050, "observer row count")
    require(row_types == Counter({"patch_observer": 2048, "cap_observer": 2}), "observer row types")


def main() -> int:
    try:
        archive = load_json("archive_manifest.json")
        require(archive["schema"] == "oph.curated_run_archive.v2", "archive schema")
        verify_inventory(archive)
        verify_reports(archive)
        verify_exact_normal_form(archive)
        verify_observer_stream(archive)
    except (OSError, KeyError, TypeError, ValueError, VerificationError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print("PASS: archive bytes, reports, observer rows, and independent authority normal form verify")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
