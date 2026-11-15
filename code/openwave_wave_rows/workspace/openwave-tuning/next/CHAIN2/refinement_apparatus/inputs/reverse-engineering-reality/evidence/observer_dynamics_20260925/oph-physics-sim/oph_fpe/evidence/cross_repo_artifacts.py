"""Hash-pinned imports from the OPH research repository.

Paper and particle artifacts are useful inputs to simulator diagnostics, but an
imported theorem status is not a run receipt.  This module deliberately keeps
those concepts separate: every imported row is informational/diagnostic until
an independent simulator verifier recomputes the hypotheses from run data.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Iterable, Mapping


MANIFEST_SCHEMA = "oph_cross_repo_artifact_manifest_v2"


@dataclass(frozen=True)
class ArtifactSpec:
    key: str
    source_relpath: str
    expected_artifact: str
    epistemic_class: str
    destination_relpath: str | None = None
    optional: bool = False
    identity_path: tuple[str, ...] = ("artifact",)

    @property
    def target_relpath(self) -> str:
        return self.destination_relpath or f"imported_oph_artifacts/{Path(self.source_relpath).name}"


DEFAULT_ARTIFACT_SPECS: tuple[ArtifactSpec, ...] = (
    ArtifactSpec(
        "geometry_cyclic_cap_net",
        "code/geometry/runs/cyclic_cap_net_run_domain.json",
        "oph_cyclic_cap_net_repair_run",
        "paper_side_finite_regression_fixture",
    ),
    ArtifactSpec(
        "geometry_modular_clock",
        "code/geometry/runs/modular_clock_instrumentation_report.json",
        "oph_modular_clock_instrumentation",
        "paper_side_one_particle_regression_fixture",
    ),
    ArtifactSpec(
        "geometry_null_net",
        "code/geometry/runs/null_net_receipt_report.json",
        "oph_null_net_receipt_instrumentation",
        "paper_side_one_particle_regression_fixture",
    ),
    ArtifactSpec(
        "geometry_realized_events",
        "code/geometry/runs/realized_event_receipt_report.json",
        "oph_realized_event_receipts",
        "paper_side_finite_regression_fixture",
    ),
    ArtifactSpec(
        "geometry_bulk_depth",
        "code/geometry/runs/bulk_depth_receipt_report.json",
        "oph_bulk_depth_receipts",
        "paper_side_finite_regression_fixture",
    ),
    ArtifactSpec(
        "einstein_realized_branch",
        "code/geometry/runs/realized_branch_receipt_report.json",
        "einstein_branch_realized_receipt_evaluation",
        "paper_side_finite_witness_status",
        "realized_branch_receipt_report.json",
    ),
    ArtifactSpec(
        "neutrino_lane_closure",
        "code/particles/runs/neutrino/neutrino_lane_closure_contract.json",
        "oph_neutrino_lane_closure_contract",
        "rejected_target_informed_candidate_status",
    ),
    ArtifactSpec(
        "neutrino_nufit61_score",
        "code/particles/runs/neutrino/nufit61_weighted_cycle_retrospective_score.json",
        "oph_neutrino_nufit61_retrospective_profile_score",
        "retrospective_rejection_evidence",
    ),
    ArtifactSpec(
        "conditional_ew_envelope",
        "code/particles/runs/calibration/conditional_ew_predictions_current.json",
        "oph_conditional_ew_predictions",
        "conditional_compare_only_envelope",
    ),
    ArtifactSpec(
        "d10_repair_selection",
        "code/particles/runs/calibration/d10_repair_tuple_selection_theorem.json",
        "oph_d10_repair_tuple_selection_theorem",
        "conditional_selection_theorem",
    ),
    ArtifactSpec(
        "color_amplitude_loop_split",
        "code/particles/runs/calibration/color_amplitude_loop_split.json",
        "oph_color_amplitude_loop_split",
        "diagnostic_color_transport_contract",
    ),
    ArtifactSpec(
        "empirical_hadron_spectral_measure",
        "code/particles/runs/hadron/empirical_ward_projected_spectral_measure.json",
        "oph_empirical_ward_projected_hadronic_spectral_measure",
        "empirical_external_data_closure",
    ),
    ArtifactSpec(
        "empirical_thomson_endpoint",
        "code/P_derivation/runtime/empirical_thomson_endpoint_current.json",
        "oph_empirical_thomson_endpoint",
        "empirical_external_data_closure",
    ),
    ArtifactSpec(
        "anchor_scheme_bridge",
        "code/P_derivation/runtime/anchor_scheme_bridge_current.json",
        "oph_anchor_scheme_bridge",
        "compare_only_scheme_analysis",
    ),
    ArtifactSpec(
        "particle_status_table",
        "code/particles/runs/status/status_table_forward_current.json",
        "oph_status_table_forward_current",
        "particle_frontier_status",
    ),
    ArtifactSpec(
        "particle_gap_ledger",
        "code/particles/runs/status/particle_derivation_gap_ledger.json",
        "oph_particle_derivation_gap_ledger",
        "open_obligation_ledger",
    ),
    ArtifactSpec(
        "charged_trace_no_go",
        "code/particles/runs/leptons/charged_trace_lift_theorem.json",
        "oph_charged_trace_lift_theorem",
        "source_nonidentifiability_no_go",
    ),
    ArtifactSpec(
        "quark_sigma_no_go",
        "code/particles/runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json",
        "oph_quark_sigma_source_nonidentifiability_obstruction",
        "source_nonidentifiability_no_go",
    ),
    ArtifactSpec(
        "quark_scheme_obstruction",
        "code/particles/runs/flavor/quark_running_mass_scheme_convention_obstruction.json",
        "oph_quark_running_mass_scheme_convention_obstruction",
        "scheme_convention_obstruction",
    ),
    ArtifactSpec(
        "lattice_engine_status",
        "code/particles/runs/hadron/lattice_engine_lane_status.json",
        "oph_lattice_engine_lane_status",
        "diagnostic_nonpromoting_engine_status",
        optional=True,
    ),
    ArtifactSpec(
        "lattice_diagnostic_export",
        "code/particles/runs/hadron/lattice_diagnostic_backend_export.json",
        "oph_lattice_diagnostic_backend_export",
        "real_lattice_diagnostic_toy_scale",
        optional=True,
    ),
    ArtifactSpec(
        "a5_screen_sm_closure",
        "code/a5_closure/a5_screen_sm_closure.json",
        "OPH A5 screen/SM closure certificate v1",
        "exact_conditional_a5_coefficient_module_witness",
        identity_path=("schema",),
    ),
    ArtifactSpec(
        "a5_exterior_sm_completion",
        "code/a5_closure/exterior_sm_completion.json",
        "OPH conditional exterior Standard-Model completion certificate v1",
        "exact_conditional_exterior_matter_witness",
        identity_path=("schema",),
    ),
    ArtifactSpec(
        "p_global_uniqueness_extension",
        "code/P_derivation/runtime/p_global_uniqueness_extension_certificate_2026-07-17.json",
        "oph_p_global_uniqueness_extension_certificate",
        "declared_map_global_uniqueness_nonpromoting",
    ),
    ArtifactSpec(
        "reversible_public_checkpoint_control",
        "code/capacity_readback/runtime/reversible_public_checkpoint_packet_certificate.json",
        "PUBLIC_CHECKPOINT_PACKET/v1-reversible",
        "exact_finite_reference_control_not_physical_receipt",
    ),
    ArtifactSpec(
        "scr330_radial_contract_example",
        "code/cosmology/example_radial_receipt.json",
        "scr330-radial-v2",
        "conditional_radial_contract_example_nonpromoting",
        identity_path=("schema_version",),
    ),
    ArtifactSpec(
        "regulator_gluing_fixed_cutoff",
        "code/regulator_gluing/runs/regulator_gluing_evidence_bundle_current.json",
        "oph.regulator_gluing.evidence_bundle/1",
        "exact_fixed_cutoff_regulator_gluing_witness",
        identity_path=("bundle", "schema"),
    ),
)


def import_cross_repo_artifacts(
    source_repo: Path,
    destination: Path,
    *,
    specs: Iterable[ArtifactSpec] = DEFAULT_ARTIFACT_SPECS,
) -> dict[str, Any]:
    source_root = Path(source_repo).resolve()
    target_root = Path(destination).resolve()
    target_root.mkdir(parents=True, exist_ok=True)
    provenance = repository_provenance(source_root)
    rows: list[dict[str, Any]] = []
    missing: list[str] = []

    for spec in specs:
        source = source_root / spec.source_relpath
        if not source.is_file():
            if not spec.optional:
                missing.append(spec.key)
            rows.append(_missing_row(spec))
            continue
        raw = source.read_bytes()
        payload = _json_object(raw, source)
        artifact_id = _nested_identity(payload, spec.identity_path)
        if artifact_id != spec.expected_artifact:
            raise ValueError(
                f"artifact id mismatch for {spec.key}: expected {spec.expected_artifact!r}, got {artifact_id!r}"
            )
        target = target_root / spec.target_relpath
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        digest = _sha256_bytes(raw)
        if _sha256_file(target) != digest:
            raise RuntimeError(f"copy verification failed for {spec.key}")
        rows.append(
            {
                "key": spec.key,
                "artifact": artifact_id,
                "identity_path": list(spec.identity_path),
                "source_relpath": spec.source_relpath,
                "target_relpath": spec.target_relpath,
                "sha256": digest,
                "byte_count": len(raw),
                "source_tracked": _git_path_tracked(source_root, spec.source_relpath),
                "generated_utc": payload.get("generated_utc"),
                "format_version": payload.get("format_version") or payload.get("schema_version"),
                "row_class": payload.get("row_class"),
                "status": _short_status(payload),
                "epistemic_class": spec.epistemic_class,
                "declared_promotion_allowed": _declared_promotion_allowed(payload),
                "simulation_receipt_eligible": False,
                "import_role": "informational_or_diagnostic_only",
                "present": True,
            }
        )

    manifest: dict[str, Any] = {
        "schema": MANIFEST_SCHEMA,
        "source_repository": provenance,
        "paper_release_id": _paper_release_id(source_root / "paper" / "release_info.tex"),
        "artifacts": rows,
        "missing_required_artifacts": missing,
        "all_required_artifacts_present": not missing,
        "all_artifact_hashes_verified": all(row.get("sha256") for row in rows if row.get("present")),
        "source_worktree_clean": provenance["dirty"] is False,
        "run_receipts_promoted_by_import": False,
        "claim_boundary": (
            "Imported paper and particle artifacts carry theorem/frontier provenance only. "
            "They never promote a simulation receipt; local verifiers must recompute every run-specific hypothesis."
        ),
    }
    manifest["manifest_payload_sha256"] = _sha256_json(manifest)
    manifest_path = target_root / "oph_cross_repo_artifact_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def verify_cross_repo_artifact_manifest(root: Path) -> dict[str, Any]:
    base = Path(root)
    path = base / "oph_cross_repo_artifact_manifest.json"
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return _verification(False, ["manifest_missing_or_invalid"], {})
    if not isinstance(manifest, dict):
        return _verification(False, ["manifest_not_object"], {})
    blockers: list[str] = []
    if manifest.get("schema") != MANIFEST_SCHEMA:
        blockers.append("manifest_schema_invalid")
    if manifest.get("all_required_artifacts_present") is not True:
        blockers.append("required_artifacts_missing")
    if manifest.get("source_worktree_clean") is not True:
        blockers.append("source_worktree_not_clean")
    missing_required = manifest.get("missing_required_artifacts")
    if not isinstance(missing_required, list) or missing_required:
        blockers.append("missing_required_artifact_ledger_nonempty_or_invalid")
    expected_payload_hash = manifest.get("manifest_payload_sha256")
    unhashed = dict(manifest)
    unhashed.pop("manifest_payload_sha256", None)
    if expected_payload_hash != _sha256_json(unhashed):
        blockers.append("manifest_payload_hash_mismatch")
    rows = manifest.get("artifacts")
    if not isinstance(rows, list):
        blockers.append("artifact_rows_missing")
        rows = []
    verified_rows: list[dict[str, Any]] = []
    for raw_row in rows:
        if not isinstance(raw_row, dict):
            blockers.append("artifact_row_not_object")
            continue
        row = dict(raw_row)
        if row.get("present") is not True:
            verified_rows.append({**row, "hash_verified": False})
            continue
        relpath = row.get("target_relpath")
        expected = row.get("sha256")
        candidate = base / str(relpath)
        actual = _sha256_file(candidate) if candidate.is_file() else None
        valid = bool(isinstance(expected, str) and expected == actual)
        if not valid:
            blockers.append(f"artifact_hash_mismatch:{row.get('key', 'unknown')}")
        if row.get("simulation_receipt_eligible") is not False:
            blockers.append(f"import_attempted_receipt_promotion:{row.get('key', 'unknown')}")
        verified_rows.append({**row, "hash_verified": valid})
    if manifest.get("run_receipts_promoted_by_import") is not False:
        blockers.append("manifest_attempted_run_receipt_promotion")
    return _verification(not blockers, blockers, {**manifest, "artifacts": verified_rows})


def repository_provenance(repo: Path) -> dict[str, Any]:
    root = Path(repo).resolve()
    commit = _git(root, "rev-parse", "HEAD") or "unknown"
    status = _git(root, "status", "--porcelain=v1", "--untracked-files=all")
    dirty = None if status is None else bool(status)
    patch = _git_bytes(root, "diff", "--binary", "HEAD", "--")
    untracked = _git(root, "ls-files", "--others", "--exclude-standard")
    untracked_rows = sorted(line for line in (untracked or "").splitlines() if line)
    material = bytearray(patch or b"")
    for relpath in untracked_rows:
        material.extend(relpath.encode("utf-8") + b"\0")
        candidate = root / relpath
        if candidate.is_file():
            material.extend(candidate.read_bytes())
    return {
        "checkout_name": root.name,
        "commit": commit,
        "dirty": dirty,
        "worktree_state_sha256": _sha256_bytes(bytes(material)),
        "untracked_file_count": len(untracked_rows),
    }


def _verification(ok: bool, blockers: list[str], manifest: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": "oph_cross_repo_artifact_verification_v1",
        "verified": bool(ok),
        "blockers": blockers,
        "manifest": manifest,
        "simulation_receipts_promoted": False,
    }


def _missing_row(spec: ArtifactSpec) -> dict[str, Any]:
    return {
        "key": spec.key,
        "artifact": spec.expected_artifact,
        "identity_path": list(spec.identity_path),
        "source_relpath": spec.source_relpath,
        "target_relpath": spec.target_relpath,
        "epistemic_class": spec.epistemic_class,
        "optional": spec.optional,
        "present": False,
        "simulation_receipt_eligible": False,
    }


def _json_object(raw: bytes, path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid JSON artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"artifact must be a JSON object: {path}")
    return payload


def _nested_identity(payload: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = payload
    for key in path:
        if not isinstance(value, Mapping):
            return None
        value = value.get(key)
    return value


def _declared_promotion_allowed(payload: dict[str, Any]) -> bool | None:
    candidates: list[Any] = []
    for key in ("public_promotion_allowed", "promotion_allowed", "bridge_prediction_promotion_allowed"):
        if key in payload:
            candidates.append(payload[key])
    guards = payload.get("guards")
    if isinstance(guards, dict):
        for key in (
            "public_promotion_allowed",
            "promotable_as_oph_source_theorem",
            "diagnostic_output_promotable",
        ):
            if key in guards:
                candidates.append(guards[key])
    if not candidates:
        return None
    return all(value is True for value in candidates)


def _short_status(payload: dict[str, Any]) -> Any:
    status = payload.get("status") or payload.get("proof_status") or payload.get("verdict")
    if isinstance(status, dict):
        return status.get("status") or status.get("decision") or status.get("statement")
    if isinstance(status, str) and len(status) > 500:
        return status[:500]
    return status


def _git_path_tracked(repo: Path, relpath: str) -> bool:
    result = subprocess.run(
        ["git", "ls-files", "--error-unmatch", "--", relpath],
        cwd=repo,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def _paper_release_id(path: Path) -> str | None:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    match = re.search(r"\\OPHPaperReleaseID\}\{([^}]+)\}", text)
    return match.group(1) if match else None


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _sha256_bytes(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def _sha256_json(payload: dict[str, Any]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return _sha256_bytes(raw)


def _git(repo: Path, *args: str) -> str | None:
    raw = _git_bytes(repo, *args)
    return raw.decode("utf-8", errors="replace").strip() if raw is not None else None


def _git_bytes(repo: Path, *args: str) -> bytes | None:
    try:
        return subprocess.check_output(["git", *args], cwd=repo, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError):
        return None
