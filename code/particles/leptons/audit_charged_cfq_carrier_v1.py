#!/usr/bin/env python3
"""Audit the submitted finite CFQ carrier as a digital model witness.

The archive is never extracted or executed by this canonical auditor.  It
verifies archive integrity, the exported graph/register/path data, the central
record dilation, model-relative CFQ receipts, negative controls, and agreement
with the canonical conditional packet theorem.

The crucial boundary is preserved: constructing a finite model whose automaton
contains exactly the desired paths proves existence inside that declared model.
It does not prove that current OPH dynamics or nature selects the model as the
physical charged sector.  The submitted dependency DAG also cannot certify
historical no-target ancestry for a packet devised after the charged spectrum
was known.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import stat
from pathlib import Path, PurePosixPath
from typing import Any
from zipfile import ZipFile, ZipInfo


HERE = Path(__file__).resolve()
CODE_ROOT = HERE.parents[2]
REPO = HERE.parents[3]
WORKSPACE = HERE.parents[4]
DEFAULT_ARCHIVE = WORKSPACE / "correspondence" / "oph_charged_cfq_carrier_v1.zip"
CANONICAL_CFQ = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_source_law_rigidity_conditional.json"
)
DEFAULT_OUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_cfq_carrier_v1_review.json"
)

PREFIX = "oph_charged_cfq_carrier_v1/"
EXPECTED_ARCHIVE_SHA256 = (
    "d4b47d158f80b4bffe318650769812a606e4e065cab1e1a01edd68ff8909291c"
)
EXPECTED_NARRATIVE_SHA256 = (
    "de9fac0fab1b61e5870a9ea67c737b14cbeb5462cc18c2b6874fefd713045d33"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "8dd66843603d25eba2380fef37ce8b6c5eaed2cd12902b0b3b75e8b3231ab3aa"
)
EXPECTED_VERIFICATION_SHA256 = (
    "8612131867c8a0e2203103237444e4e3abcabc4799facae788b5100a0444abb1"
)
EXPECTED_NEGATIVE_CONTROL_SHA256 = (
    "206f24ca93c3da40dc89c543acb3b541f6f6c817890fd35a3cefe15215f6feb4"
)
EXPECTED_FINAL_STATUS_SHA256 = (
    "d3c92f36199e3abdc371606808b05b1a00c37c31e13b69966a354458a69e2331"
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_member(info: ZipInfo) -> bool:
    path = PurePosixPath(info.filename)
    mode = info.external_attr >> 16
    return (
        not path.is_absolute()
        and ".." not in path.parts
        and not stat.S_ISLNK(mode)
        and not (info.flag_bits & 0x1)
    )


def parse_hashes(raw: bytes) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in raw.decode("utf-8").splitlines():
        digest, name = line.split(maxsplit=1)
        rows[name.strip()] = digest
    return rows


def graph_connected(dimension: int, edges: list[list[int]]) -> bool:
    if dimension <= 0:
        return False
    adjacency = [set() for _ in range(dimension)]
    for left, right in edges:
        if not (0 <= left < dimension and 0 <= right < dimension) or left == right:
            return False
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen: set[int] = set()
    pending = [0]
    while pending:
        node = pending.pop()
        if node in seen:
            continue
        seen.add(node)
        pending.extend(adjacency[node] - seen)
    return len(seen) == dimension


def build_review(archive_path: Path, canonical_cfq_path: Path) -> dict[str, Any]:
    archive_raw = archive_path.read_bytes()
    canonical_raw = canonical_cfq_path.read_bytes()
    canonical = json.loads(canonical_raw)

    with ZipFile(archive_path) as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        first_bad_crc = archive.testzip()
        hashes = parse_hashes(archive.read(PREFIX + "SHA256SUMS"))
        internal_hash_checks = {
            name: sha256(archive.read(PREFIX + name)) == digest
            for name, digest in hashes.items()
        }

        def raw(name: str) -> bytes:
            return archive.read(PREFIX + name)

        def packet(name: str) -> dict[str, Any]:
            return json.loads(raw(name))

        def gzip_packet(name: str) -> dict[str, Any]:
            return json.loads(gzip.decompress(raw(name)))

        certificate_raw = raw("cfq_carrier_certificate.json")
        verification_raw = raw("cfq_verification_receipt.json")
        negative_raw = raw("negative_control_receipt.json")
        final_raw = raw("FINAL_STATUS.json")
        certificate = json.loads(certificate_raw)
        verification = json.loads(verification_raw)
        negative = json.loads(negative_raw)
        final_status = json.loads(final_raw)
        carrier = packet("carrier_packet.json")
        dependency = packet("exports/dependency_dag.json")
        patch_tuple = packet("exports/oph_observer_patch_tuple.json")
        path_list = packet("exports/path_list.json")
        automaton = packet("exports/path_automaton.json")
        recovery = packet("exports/recovery_maps.json")
        refinement = packet("exports/refinement_maps.json")
        register_algebras = packet("exports/register_algebras.json")
        graph_packet = gzip_packet("exports/transition_graphs.json.gz")
        basis_labels = gzip_packet("exports/register_basis_labels.json.gz")
        rotations = gzip_packet("exports/rotation_intertwiners.json.gz")
        matrix_units = [
            json.loads(line)
            for line in gzip.decompress(
                raw("exports/transition_matrix_units.jsonl.gz")
            ).decode("utf-8").splitlines()
            if line
        ]

    expected_names = {PREFIX + name for name in hashes} | {PREFIX + "SHA256SUMS"}
    archive_checks = {
        "archive_hash_matches": sha256(archive_raw) == EXPECTED_ARCHIVE_SHA256,
        "safe_non_symlink_unencrypted_members": all(safe_member(info) for info in infos),
        "no_duplicate_member_names": len(names) == len(set(names)),
        "exact_member_set": set(names) == expected_names,
        "crc_check_pass": first_bad_crc is None,
        "forty_five_internal_hashes_declared": len(hashes) == 45,
        "all_internal_hashes_match": all(internal_hash_checks.values()),
        "certificate_hash_matches": sha256(certificate_raw) == EXPECTED_CERTIFICATE_SHA256,
        "verification_hash_matches": sha256(verification_raw) == EXPECTED_VERIFICATION_SHA256,
        "negative_control_hash_matches": sha256(negative_raw) == EXPECTED_NEGATIVE_CONTROL_SHA256,
        "final_status_hash_matches": sha256(final_raw) == EXPECTED_FINAL_STATUS_SHA256,
    }

    expected_dimensions = {
        name: row["dimension"]
        for name, row in canonical["conditional_cfq_packet"]["registers"].items()
    }
    expected_graph_edges = {
        name: row["edges"]
        for name, row in canonical["declared_finite_graph_models"]["graphs"].items()
    }
    graph_checks: dict[str, dict[str, Any]] = {}
    for name, expected_dimension in expected_dimensions.items():
        graph = graph_packet[name]
        algebra = register_algebras["registers"][name]
        labels = basis_labels[name]
        graph_checks[name] = {
            "dimension_matches": graph["dimension"] == expected_dimension,
            "basis_label_count_matches": len(labels) == expected_dimension,
            "edge_count_matches": len(graph["edges"]) == expected_graph_edges[name],
            "connected": graph_connected(graph["dimension"], graph["edges"]),
            "algebra_dimension_matches": algebra["dimension"] == expected_dimension,
            "normalized_trace_matches": (
                algebra["minimal_projector_trace"] == f"1/{expected_dimension}"
            ),
        }

    unit_counts = Counter(row["register"] for row in matrix_units)
    matrix_unit_checks = {
        name: unit_counts[name]
        == expected_dimensions[name] + 2 * expected_graph_edges[name]
        for name in expected_dimensions
    }

    canonical_paths = {
        row["name"]: row for row in canonical["conditional_cfq_packet"]["paths"]
    }
    submitted_paths = {row["id"]: row for row in path_list["paths"]}
    path_checks: dict[str, bool] = {}
    for name, expected in canonical_paths.items():
        submitted = submitted_paths[name]
        coupling = submitted["coupling"]
        monomial = {
            (1, 0): "alpha_U",
            (2, 0): "alpha_U^2",
            (0, 2): "alpha_P^2",
        }[(coupling["alpha_U"], coupling["alpha_P"])]
        if submitted["state_factor"] is not None:
            monomial += f"*{submitted['state_factor']}"
        path_checks[name] = (
            submitted["block"] == expected["block"]
            and submitted["registers"] == expected["registers"]
            and submitted["sign"] == expected["sign"]
            and submitted["signed_trace_weight"] == expected["exact_trace_coefficient"]
            and submitted["checkpoint_ticks"] == 1
            and monomial == expected["monomial"]
        )

    submitted_gate_passes = {
        name: row["pass"] for name, row in verification["gates"].items()
    }
    model_checks = {
        "submitted_verifier_reports_all_eight_gates_pass": (
            len(submitted_gate_passes) == 8 and all(submitted_gate_passes.values())
        ),
        "certificate_reports_all_eight_gates_pass": (
            len(certificate["gates"]) == 8 and all(certificate["gates"].values())
        ),
        "all_register_graph_exports_reproduce": all(
            all(row.values()) for row in graph_checks.values()
        ),
        "all_6467_local_matrix_units_accounted_for": (
            len(matrix_units) == 6467 and all(matrix_unit_checks.values())
        ),
        "all_eight_path_rows_match_conditional_packet": (
            len(path_checks) == 8 and all(path_checks.values())
        ),
        "automaton_declares_exactly_the_eight_packet_paths": (
            set(automaton["declared_primitive_paths"]) == set(canonical_paths)
        ),
        "central_record_dilation_exported_for_ten_occurrences": (
            len(recovery["occurrences"]) == 10
            and all(
                row["central_record_output"]
                and row["cptp"]
                and row["idempotent"]
                for row in recovery["occurrences"]
            )
        ),
        "finite_patch_tuple_declares_central_record_algebra": (
            patch_tuple["record_algebra"]["central"] is True
        ),
        "sixty_rotation_charts_exported": (
            rotations["order"] == 60
            and len(rotations["rotations"]) == 60
            and rotations["face_orbit_size"] == 20
            and rotations["selected_face_stabilizer_size"] == 3
        ),
        "inert_refinement_models_exported": (
            refinement["refinement_kind"] == "inert ancillary tensor refinement"
        ),
        "all_eight_targeted_negative_controls_rejected": (
            negative["all_mutations_rejected"] is True
            and len(negative["mutations"]) == 8
            and all(row["mutation_rejected"] for row in negative["mutations"])
        ),
        "source_law_text_matches_conditional_packet": (
            carrier["source_law"] == verification["derived_source_law"]
        ),
        "submitted_runtime_has_no_external_input_files": (
            carrier["external_runtime_inputs"] == []
            and verification["external_runtime_inputs"] == []
        ),
    }

    historical_ancestry_checks = {
        "submitted_dependency_dag_is_internally_acyclic_and_claims_no_forbidden_node": (
            dependency["forbidden_nodes_present"] is False
            and dependency["path_from_forbidden_node_to_output"] is False
            and dependency["runtime_external_input_files"] == []
        ),
        "dependency_dag_omits_target_informed_parent_history": (
            "historical_formula_discovery" not in {
                row["id"] for row in dependency["nodes"]
            }
            and canonical["historical_charged_target_informed"] is True
        ),
        "strict_historical_no_target_ancestry_pass": False,
        "physical_charged_sector_identified_with_constructed_patch": False,
        "five_oph_axioms_uniquely_select_constructed_patch": False,
    }

    return {
        "artifact": "oph_charged_cfq_carrier_v1_review",
        "schema_version": 1,
        "status": (
            "EXPLICIT_FINITE_DIGITAL_CFQ_MODEL_VERIFIED_"
            "PHYSICAL_CHARGED_SOURCE_SELECTION_AND_HISTORICAL_ANCESTRY_OPEN"
        ),
        "compare_only": True,
        "forbidden_as_candidate_ancestor": True,
        "runtime_charged_reference_consumed": False,
        "historical_charged_target_informed": True,
        "global_source_only": False,
        "branch_tuple_coherent": False,
        "mass_scheme_certified": False,
        "public_prediction_promotion_allowed": False,
        "provenance": {
            "archive_path": (
                str(archive_path.relative_to(WORKSPACE))
                if archive_path.is_relative_to(WORKSPACE)
                else archive_path.name
            ),
            "archive_sha256": sha256(archive_raw),
            "narrative_attachment_sha256": EXPECTED_NARRATIVE_SHA256,
            "canonical_conditional_receipt": str(canonical_cfq_path.relative_to(REPO)),
            "canonical_conditional_receipt_sha256": sha256(canonical_raw),
        },
        "archive": {
            "member_count": len(infos),
            "crc_first_bad_member": first_bad_crc,
            "internal_hash_checks": internal_hash_checks,
            "checks": archive_checks,
            "checks_pass": all(archive_checks.values()),
            "execution_boundary": "Archive programs were not executed by this canonical auditor.",
        },
        "verified_digital_model": {
            "explicit_finite_model_exists": True,
            "model_relative_submitted_cfq_gate_count": sum(submitted_gate_passes.values()),
            "register_graphs": graph_checks,
            "matrix_unit_counts": dict(sorted(unit_counts.items())),
            "matrix_unit_checks": matrix_unit_checks,
            "paths": path_checks,
            "central_record_resolution": (
                "Each noncentral rank-one event is dilated to an accepted/rejected central "
                "D2 outcome in the constructed model."
            ),
            "path_exhaustion_scope": (
                "Exhaustion is exact relative to the automaton whose nodes and edges are "
                "declared by the builder; the automaton is not emitted by prior OPH dynamics."
            ),
            "checks": model_checks,
            "checks_pass": all(model_checks.values()),
        },
        "historical_and_physical_boundary": {
            "checks": historical_ancestry_checks,
            "submitted_no_target_ancestry_claim_accepted": False,
            "reason": (
                "The dependency DAG begins with structural integers and omits the actual "
                "target-informed history of the register dimensions and path catalogue. "
                "Runtime reference separation cannot repair historical ancestry."
            ),
            "digital_existence_consequence": (
                "The broad finite OPH interface admits a deliberately constructed patch "
                "realizing the CFQ schema."
            ),
            "nonconsequence": (
                "Neither current OPH dynamics nor the five broad axioms select this patch "
                "as the physical charged response."
            ),
        },
        "independent_reproduction_observation": {
            "semantic_rebuild_passed_all_submitted_gates": True,
            "negative_controls_rejected_all_eight_mutations": True,
            "byte_reproducible_across_tested_python_numpy_environment": False,
            "submitted_environment": "Python 3.13.5, NumPy 2.3.5, NetworkX 3.6.1",
            "audit_environment": "Python 3.13.7, NumPy 2.3.5, NetworkX 3.6.1",
            "cause": (
                "The build receipt embeds its absolute output path and Python version; "
                "rotation determinant serialization also changed between tested runtimes."
            ),
        },
        "still_open": [
            "derive or uniquely select the CFQ carrier from OPH dynamics before charged comparison",
            "establish honest historical no-target ancestry for any future carrier law",
            "identify the physical charged-generation space with the face-corner bundle",
            "attach digital trace weights to deterministic charged Yukawa response powers",
            "derive rather than define the sqrt(2), 2/9, and 6^-14 baselines",
            "attach the determinant character to the renormalized charged determinant line",
            "freeze one coherent source branch and derive the QED/electroweak pole scheme",
        ],
        "integration_decision": (
            "Retain as a substantial explicit finite digital-model witness and central-record "
            "construction. Do not mark physical CFQ realization, OPH source-law selection, "
            "historical blindness, or charged mass prediction closed."
        ),
        "checks_pass": all(archive_checks.values()) and all(model_checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--canonical-cfq", type=Path, default=CANONICAL_CFQ)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    review = build_review(args.archive.resolve(), args.canonical_cfq.resolve())
    encoded = (json.dumps(review, indent=2, sort_keys=True) + "\n").encode()
    if args.check:
        actual = args.out.read_bytes() if args.out.exists() else None
        ok = actual == encoded
        print(json.dumps({"status": "OK" if ok else "DRIFT"}, indent=2))
        return 0 if ok else 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": review["status"],
                "checks_pass": review["checks_pass"],
                "model_relative_cfq_gates": review["verified_digital_model"][
                    "model_relative_submitted_cfq_gate_count"
                ],
                "submitted_no_target_ancestry_claim_accepted": review[
                    "historical_and_physical_boundary"
                ]["submitted_no_target_ancestry_claim_accepted"],
                "public_prediction_promotion_allowed": review[
                    "public_prediction_promotion_allowed"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
