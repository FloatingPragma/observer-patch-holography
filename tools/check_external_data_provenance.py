#!/usr/bin/env python3
"""Validate the issue-553 external-data provenance registry.

The registry pins local artifacts and their loaders. It also distinguishes
raw files whose published bytes are known from values that were manually
transcribed or normalized from a live API. The latter classes must expose a
declared provenance gap instead of inventing an upstream-file hash.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path, PurePosixPath
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = ROOT / "code/audit/external_data_provenance_registry.json"
SCHEMA = "oph.external_data_provenance_registry.v1"
SHA256_RE = re.compile(r"[0-9a-f]{64}")

REQUIRED_DATASET_ARTIFACTS = {
    "pdg-2026-particle-reference-values":
        "code/particles/data/particle_reference_values.json",
    "pdg-legacy-particle-masses-json": "pdg_data/particle_masses.json",
    "pdg-legacy-particle-masses-csv": "pdg_data/particle_masses.csv",
    "knt19-pdg2025-hadronic-spectral-shape":
        "code/particles/runs/hadron/empirical_ee_hadronic_spectral_measure.json",
    "planck2018-table2-lambda-to-N-gaussian-approximation":
        "code/capacity_readback/planck_posterior/planck_lambda_to_N_propagation.json",
    "nufit-6.1-profile-source-manifest":
        "code/particles/neutrino/nufit61_sources.json",
    "flag-2024-light-quark-ratio-fixture":
        "code/particles/data/flag_2024_light_quark_ratio_fixture.json",
    "pdg-2024-vus-kmu2-compare-only-fixture":
        "code/particles/data/pdg_2024_vus_kmu2_fixture.json",
    "codata-2022-inverse-fine-structure-constant":
        "code/P_derivation/codata_2022_alpha_fixture.json",
    "pdg-2026-wz-running-width-target-fixture":
        "code/particles/calibration/wz_pdg_2026_target_fixture.json",
    "bouchard-donagi-threshold-spectrum-literature-packet":
        "code/particles/data/oph_bd_threshold_spectrum_inputs.json",
    "auger-2022-fz12-photon-threshold-diagnostic":
        "code/a5_fingerprint/runtime/fz12_auger_threshold_diagnostic_receipt.json",
}

TOP_LEVEL_FIELDS = {"schema", "issue", "status", "policy", "entries"}
POLICY_FIELDS = {
    "artifact_rule",
    "license_rule",
    "raw_input_rule",
    "determinism_rule",
    "allowed_loader_classifications",
    "allowed_raw_input_classifications",
}
ENTRY_FIELDS = {
    "dataset_id",
    "artifact_path",
    "artifact_sha256",
    "artifact_bytes",
    "artifact_role",
    "provenance_status",
    "source",
    "license",
    "raw_inputs",
    "loader",
}
SOURCE_FIELDS = {"publisher", "version", "citation", "urls"}
LICENSE_FIELDS = {"expression", "note"}
RAW_INPUT_FIELDS = {
    "classification",
    "raw_payloads_vendored",
    "upstream_files",
    "provenance_gap",
}
UPSTREAM_FILE_FIELDS = {"id", "url", "sha256", "bytes"}
LOADER_FIELDS = {
    "path",
    "sha256",
    "bytes",
    "classification",
    "network_required",
    "deterministic_from_declared_inputs",
    "declared_variable_metadata_json_pointers",
    "requirements",
}
NO_RAW_HASH_CLASSES = {
    "hand_transcribed_published_constants_no_raw_payload",
    "live_api_snapshot_raw_responses_not_archived",
}
DETERMINISTIC_LOADER_CLASSES = {
    "deterministic_hand_transcribed_gaussian_approximation",
    "deterministic_transcribed_constants_generator",
    "hash_verifying_external_table_scorer",
}
NETWORK_LOADER_CLASSES = {
    "legacy_optional_live_api_refresh",
    "live_api_refresh_snapshot",
}


class ProvenanceError(ValueError):
    """Raised when the provenance registry fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ProvenanceError(message)


def _require_exact_fields(
    value: Mapping[str, Any],
    expected: set[str],
    path: str,
) -> None:
    actual = set(value)
    _require(
        actual == expected,
        (
            f"{path} fields must be exactly {sorted(expected)}; "
            f"missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        ),
    )


def _require_nonempty_string(value: Any, path: str) -> str:
    _require(isinstance(value, str) and bool(value.strip()), f"{path} must be a non-empty string")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _repository_path(root: Path, raw: Any, field: str) -> Path:
    text = _require_nonempty_string(raw, field)
    pure = PurePosixPath(text)
    _require(not pure.is_absolute(), f"{field} must be repository-relative")
    _require(".." not in pure.parts, f"{field} must not traverse above the repository")
    resolved_root = root.resolve()
    resolved = (root / Path(*pure.parts)).resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise ProvenanceError(f"{field} resolves outside the repository") from exc
    _require(resolved.is_file(), f"{field} does not name a file: {text}")
    return resolved


def _validate_file_pin(
    root: Path,
    record: Mapping[str, Any],
    *,
    path_field: str,
    hash_field: str,
    bytes_field: str,
    label: str,
) -> Path:
    path = _repository_path(root, record.get(path_field), f"{label}.{path_field}")
    expected_hash = record.get(hash_field)
    expected_bytes = record.get(bytes_field)
    _require(
        isinstance(expected_hash, str) and SHA256_RE.fullmatch(expected_hash) is not None,
        f"{label}.{hash_field} must be a lowercase SHA-256",
    )
    _require(
        type(expected_bytes) is int and expected_bytes >= 0,
        f"{label}.{bytes_field} must be a non-negative JSON integer",
    )
    _require(
        path.stat().st_size == expected_bytes,
        f"{label} byte mismatch for {path.relative_to(root)}: "
        f"registry={expected_bytes}, disk={path.stat().st_size}",
    )
    actual_hash = _sha256(path)
    _require(
        actual_hash == expected_hash,
        f"{label} SHA-256 mismatch for {path.relative_to(root)}: "
        f"registry={expected_hash}, disk={actual_hash}",
    )
    return path


def _validate_https_urls(urls: Any, path: str) -> list[str]:
    _require(isinstance(urls, list) and bool(urls), f"{path} must be a non-empty list")
    output: list[str] = []
    for index, value in enumerate(urls):
        url = _require_nonempty_string(value, f"{path}[{index}]")
        _require(url.startswith("https://"), f"{path}[{index}] must use HTTPS")
        output.append(url)
    return output


def _validate_upstream_files(raw: Mapping[str, Any], path: str) -> None:
    classification = raw["classification"]
    upstream = raw["upstream_files"]
    _require(isinstance(upstream, list), f"{path}.upstream_files must be a list")
    if classification in NO_RAW_HASH_CLASSES:
        _require(
            upstream == [],
            (
                f"{path}.upstream_files must be empty for {classification}; "
                "do not assign raw-file hashes to a transcription or an "
                "unarchived live response"
            ),
        )
    if classification == "hash_pinned_external_files_not_vendored":
        _require(bool(upstream), f"{path}.upstream_files must contain the published file pins")

    seen_ids: set[str] = set()
    for index, item in enumerate(upstream):
        item_path = f"{path}.upstream_files[{index}]"
        _require(isinstance(item, Mapping), f"{item_path} must be an object")
        _require_exact_fields(item, UPSTREAM_FILE_FIELDS, item_path)
        source_id = _require_nonempty_string(item["id"], f"{item_path}.id")
        _require(source_id not in seen_ids, f"{item_path}.id is duplicated")
        seen_ids.add(source_id)
        _validate_https_urls([item["url"]], f"{item_path}.url_list")
        _require(
            isinstance(item["sha256"], str)
            and SHA256_RE.fullmatch(item["sha256"]) is not None,
            f"{item_path}.sha256 must be a lowercase SHA-256",
        )
        _require(
            type(item["bytes"]) is int and item["bytes"] > 0,
            f"{item_path}.bytes must be a positive JSON integer",
        )


def _validate_nufit_manifest(entry: Mapping[str, Any], artifact_path: Path) -> None:
    if entry["dataset_id"] != "nufit-6.1-profile-source-manifest":
        return
    manifest = json.loads(artifact_path.read_text(encoding="utf-8"))
    _require(manifest.get("license") == "NOASSERTION", "NuFIT manifest license status drifted")
    expected: dict[str, tuple[str, int, str]] = {}
    files = manifest.get("files")
    _require(isinstance(files, Mapping) and bool(files), "NuFIT manifest files are missing")
    for source_id, record in files.items():
        _require(isinstance(record, Mapping), f"NuFIT source {source_id} must be an object")
        expected[str(source_id)] = (
            str(record.get("sha256")),
            int(record.get("bytes", -1)),
            str(record.get("url")),
        )
    notes = manifest.get("release_notes")
    _require(isinstance(notes, Mapping), "NuFIT release-notes pin is missing")
    expected["release-notes"] = (
        str(notes.get("sha256")),
        int(notes.get("bytes", -1)),
        str(notes.get("url")),
    )
    registered = {
        item["id"]: (item["sha256"], item["bytes"], item["url"])
        for item in entry["raw_inputs"]["upstream_files"]
    }
    _require(
        registered == expected,
        "NuFIT upstream file pins differ between the source manifest and provenance registry",
    )


def _validate_bd_manifest(entry: Mapping[str, Any], artifact_path: Path) -> None:
    if entry["dataset_id"] != "bouchard-donagi-threshold-spectrum-literature-packet":
        return
    packet = json.loads(artifact_path.read_text(encoding="utf-8"))
    expected: dict[str, tuple[str, str]] = {}
    for source in packet.get("external_sources", []):
        digest = source.get("sha256_source_archive", source.get("sha256_pdf"))
        expected[str(source["id"])] = (str(source["url"]), str(digest))
    registered = {
        item["id"]: (item["url"], item["sha256"])
        for item in entry["raw_inputs"]["upstream_files"]
    }
    _require(
        registered == expected,
        "BD upstream URL/hash pins differ between the source packet and provenance registry",
    )


def _validate_content_boundary(entry: Mapping[str, Any], artifact_path: Path) -> None:
    dataset_id = entry["dataset_id"]
    if dataset_id == "pdg-2026-particle-reference-values":
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        _require(
            str(payload.get("source", {}).get("edition")) == "2026",
            "canonical PDG reference artifact is not the declared 2026 edition",
        )
    elif dataset_id == "knt19-pdg2025-hadronic-spectral-shape":
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        _require(
            payload.get("guards", {}).get("external_cross_section_data_used") is True,
            "hadronic artifact no longer declares its empirical-input boundary",
        )
        _require(
            payload.get("source_compilation", {}).get("id") == "knt19_pinned_piecewise_v1",
            "hadronic artifact source-compilation identity drifted",
        )
    elif dataset_id == "planck2018-table2-lambda-to-N-gaussian-approximation":
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        _require(
            set(payload) == {"inputs", "combos"},
            "Planck approximation artifact shape drifted",
        )
        _require(
            set(payload.get("inputs", {}).get("combos", {}))
            == {"TTTEEE_lowE_lensing", "TTTEEE_lowE_lensing_BAO"},
            "Planck approximation likelihood-combination boundary drifted",
        )
    elif dataset_id == "codata-2022-inverse-fine-structure-constant":
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        _require(
            payload.get("claim_status") == "compare_only_empirical_input",
            "CODATA alpha fixture must remain compare-only",
        )
        _require(
            payload.get("source", {}).get("edition") == "CODATA 2022",
            "CODATA alpha fixture edition drifted",
        )
    elif dataset_id == "flag-2024-light-quark-ratio-fixture":
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        _require(
            payload.get("artifact")
            == "oph_flag_2024_light_quark_ratio_fixture"
            and payload.get("status")
            == "COMPARE_ONLY_HAND_TRANSCRIBED_REFERENCE",
            "FLAG light-quark fixture must remain compare-only",
        )
        _require(
            payload.get("derived_quantity", {}).get(
                "input_covariance_available"
            )
            is False,
            "FLAG light-quark fixture must not invent a covariance",
        )
        _require(
            payload.get("claim_boundary", {}).get(
                "oph_theory_uncertainty_supplied"
            )
            is False
            and payload.get("claim_boundary", {}).get(
                "oph_fit_or_selection_input"
            )
            is False
            and payload.get("claim_boundary", {}).get(
                "prediction_preexisted_audit"
            )
            is True
            and payload.get("claim_boundary", {}).get(
                "significance_gate_preregistered"
            )
            is False,
            "FLAG fixture crossed its no-theory-uncertainty/no-fit boundary",
        )
        rows = payload.get("averages", [])
        _require(
            {
                (
                    row.get("nf"),
                    row.get("ms_over_mud", {}).get("published_notation"),
                    row.get("mu_over_md", {}).get("published_notation"),
                )
                for row in rows
            }
            == {
                ("2+1+1", "27.227(81)", "0.465(24)"),
                ("2+1", "27.42(12)", "0.485(19)"),
            },
            "FLAG light-quark transcriptions drifted",
        )
        derived = {
            row["nf"]: row.get("derived_ms_over_md", {})
            for row in rows
        }
        _require(
            derived.get("2+1+1", {}).get("value") == "19.9437775"
            and derived.get("2+1", {}).get("value") == "20.35935",
            "FLAG derived ms/md central values drifted",
        )
    elif dataset_id == "pdg-2024-vus-kmu2-compare-only-fixture":
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        coordinate = payload.get("coordinate", {})
        boundary = payload.get("claim_boundary", {})
        source = payload.get("source", {})
        _require(
            payload.get("artifact") == "oph_pdg_2024_vus_kmu2_fixture"
            and payload.get("status")
            == "COMPARE_ONLY_HAND_TRANSCRIBED_REFERENCE",
            "PDG |V_us| fixture must remain a compare-only transcription",
        )
        _require(
            coordinate.get("determination") == "Kmu2_decay_constant_ratio"
            and coordinate.get("value") == "0.2250"
            and coordinate.get("standard_uncertainty") == "0.0004"
            and coordinate.get("published_notation") == "0.2250 +/- 0.0004",
            "PDG Kmu2 |V_us| coordinate or uncertainty drifted",
        )
        _require(
            "revised April 2024; PDF dated 31 May 2024"
            in source.get("edition", "")
            and source.get("publisher") == "Particle Data Group",
            "PDG Kmu2 |V_us| source edition drifted",
        )
        _require(
            boundary.get("comparison_only") is True
            and boundary.get("used_to_construct_or_select_axes") is False
            and boundary.get("global_ckm_fit_value") is False
            and boundary.get("oph_fit_or_selection_input") is False,
            "PDG Kmu2 |V_us| fixture crossed its compare-only boundary",
        )
    elif dataset_id == "pdg-2026-wz-running-width-target-fixture":
        target = json.loads(artifact_path.read_text(encoding="utf-8")).get(
            "experimental_target", {}
        )
        _require(
            target.get("release") == "PDG-2026-2026-06-01",
            "W/Z target fixture release drifted",
        )
        _require(
            target.get("joint_chi2_licensed") is False
            and "assumptions, not measured facts"
            in target.get("covariance_status", ""),
            "W/Z covariance/license boundary drifted",
        )
    elif dataset_id == "bouchard-donagi-threshold-spectrum-literature-packet":
        packet = json.loads(artifact_path.read_text(encoding="utf-8"))
        _require(
            packet.get("comparison_registry", {}).get("complete_oph_target") is False,
            "BD literature packet unexpectedly claims a complete OPH target",
        )
        _require(
            packet.get("bd_branch", {}).get("selected_moduli_point") is None,
            "BD literature packet unexpectedly claims a selected moduli point",
        )


def validate_registry(payload: Mapping[str, Any], root: Path = ROOT) -> dict[str, Any]:
    """Validate a parsed registry and all repository-local pins."""

    _require(isinstance(payload, Mapping), "registry root must be an object")
    _require_exact_fields(payload, TOP_LEVEL_FIELDS, "$")
    _require(payload["schema"] == SCHEMA, f"registry schema must be {SCHEMA}")
    _require(payload["issue"] == 553, "registry issue must be 553")
    _require_nonempty_string(payload["status"], "$.status")

    policy = payload["policy"]
    _require(isinstance(policy, Mapping), "$.policy must be an object")
    _require_exact_fields(policy, POLICY_FIELDS, "$.policy")
    for field in ("artifact_rule", "license_rule", "raw_input_rule", "determinism_rule"):
        _require_nonempty_string(policy[field], f"$.policy.{field}")
    allowed_loaders = policy["allowed_loader_classifications"]
    allowed_raw = policy["allowed_raw_input_classifications"]
    _require(
        isinstance(allowed_loaders, list)
        and len(allowed_loaders) == len(set(allowed_loaders))
        and all(isinstance(item, str) and item for item in allowed_loaders),
        "$.policy.allowed_loader_classifications must be a unique string list",
    )
    _require(
        isinstance(allowed_raw, list)
        and len(allowed_raw) == len(set(allowed_raw))
        and all(isinstance(item, str) and item for item in allowed_raw),
        "$.policy.allowed_raw_input_classifications must be a unique string list",
    )

    entries = payload["entries"]
    _require(isinstance(entries, list) and bool(entries), "$.entries must be a non-empty list")
    registered_artifacts = {
        entry.get("dataset_id"): entry.get("artifact_path")
        for entry in entries
        if isinstance(entry, Mapping)
    }
    _require(
        registered_artifacts == REQUIRED_DATASET_ARTIFACTS,
        "external-data inventory does not exactly match the mandatory artifact set; "
        f"missing_or_changed={sorted(set(REQUIRED_DATASET_ARTIFACTS.items()) - set(registered_artifacts.items()))}, "
        f"unexpected={sorted(set(registered_artifacts.items()) - set(REQUIRED_DATASET_ARTIFACTS.items()))}",
    )
    seen_ids: set[str] = set()
    seen_artifacts: set[str] = set()
    loader_classes: set[str] = set()
    raw_classes: set[str] = set()

    for index, entry in enumerate(entries):
        entry_path = f"$.entries[{index}]"
        _require(isinstance(entry, Mapping), f"{entry_path} must be an object")
        _require_exact_fields(entry, ENTRY_FIELDS, entry_path)
        dataset_id = _require_nonempty_string(entry["dataset_id"], f"{entry_path}.dataset_id")
        _require(dataset_id not in seen_ids, f"{entry_path}.dataset_id is duplicated")
        seen_ids.add(dataset_id)
        artifact_text = _require_nonempty_string(
            entry["artifact_path"], f"{entry_path}.artifact_path"
        )
        _require(artifact_text not in seen_artifacts, f"{entry_path}.artifact_path is duplicated")
        seen_artifacts.add(artifact_text)
        _require_nonempty_string(entry["artifact_role"], f"{entry_path}.artifact_role")
        _require_nonempty_string(entry["provenance_status"], f"{entry_path}.provenance_status")
        artifact_path = _validate_file_pin(
            root,
            entry,
            path_field="artifact_path",
            hash_field="artifact_sha256",
            bytes_field="artifact_bytes",
            label=entry_path,
        )

        source = entry["source"]
        _require(isinstance(source, Mapping), f"{entry_path}.source must be an object")
        _require_exact_fields(source, SOURCE_FIELDS, f"{entry_path}.source")
        for field in ("publisher", "version", "citation"):
            _require_nonempty_string(source[field], f"{entry_path}.source.{field}")
        _validate_https_urls(source["urls"], f"{entry_path}.source.urls")

        license_record = entry["license"]
        _require(isinstance(license_record, Mapping), f"{entry_path}.license must be an object")
        _require_exact_fields(license_record, LICENSE_FIELDS, f"{entry_path}.license")
        expression = _require_nonempty_string(
            license_record["expression"], f"{entry_path}.license.expression"
        )
        note = _require_nonempty_string(license_record["note"], f"{entry_path}.license.note")
        if expression == "NOASSERTION":
            _require(
                len(note) >= 40,
                f"{entry_path}.license.note must explain why the license is NOASSERTION",
            )

        raw = entry["raw_inputs"]
        _require(isinstance(raw, Mapping), f"{entry_path}.raw_inputs must be an object")
        _require_exact_fields(raw, RAW_INPUT_FIELDS, f"{entry_path}.raw_inputs")
        raw_class = _require_nonempty_string(
            raw["classification"], f"{entry_path}.raw_inputs.classification"
        )
        _require(
            raw_class in allowed_raw,
            f"{entry_path}.raw_inputs.classification is not policy-registered",
        )
        raw_classes.add(raw_class)
        _require(
            type(raw["raw_payloads_vendored"]) is bool,
            f"{entry_path}.raw_inputs.raw_payloads_vendored must be boolean",
        )
        _require(
            raw["raw_payloads_vendored"] is False,
            f"{entry_path} unexpectedly claims vendored raw payloads",
        )
        _require_nonempty_string(
            raw["provenance_gap"],
            f"{entry_path}.raw_inputs.provenance_gap",
        )
        _validate_upstream_files(raw, f"{entry_path}.raw_inputs")

        loader = entry["loader"]
        _require(isinstance(loader, Mapping), f"{entry_path}.loader must be an object")
        _require_exact_fields(loader, LOADER_FIELDS, f"{entry_path}.loader")
        _validate_file_pin(
            root,
            loader,
            path_field="path",
            hash_field="sha256",
            bytes_field="bytes",
            label=f"{entry_path}.loader",
        )
        loader_class = _require_nonempty_string(
            loader["classification"], f"{entry_path}.loader.classification"
        )
        _require(
            loader_class in allowed_loaders,
            f"{entry_path}.loader.classification is not policy-registered",
        )
        loader_classes.add(loader_class)
        _require(
            type(loader["network_required"]) is bool
            and type(loader["deterministic_from_declared_inputs"]) is bool,
            f"{entry_path}.loader network/determinism fields must be booleans",
        )
        _require(
            loader["network_required"] == (loader_class in NETWORK_LOADER_CLASSES),
            f"{entry_path}.loader.network_required conflicts with its classification",
        )
        _require(
            loader["deterministic_from_declared_inputs"]
            == (loader_class in DETERMINISTIC_LOADER_CLASSES),
            f"{entry_path}.loader deterministic flag conflicts with its classification",
        )
        metadata = loader["declared_variable_metadata_json_pointers"]
        requirements = loader["requirements"]
        _require(
            isinstance(metadata, list)
            and all(isinstance(item, str) and item.startswith("/") for item in metadata),
            f"{entry_path}.loader.declared_variable_metadata_json_pointers must contain JSON pointers",
        )
        _require(
            isinstance(requirements, list)
            and bool(requirements)
            and all(isinstance(item, str) and item.strip() for item in requirements),
            f"{entry_path}.loader.requirements must be a non-empty string list",
        )

        _validate_nufit_manifest(entry, artifact_path)
        _validate_bd_manifest(entry, artifact_path)
        _validate_content_boundary(entry, artifact_path)

    _require(
        loader_classes == set(allowed_loaders),
        "policy loader classifications must exactly match the classifications in use",
    )
    _require(
        raw_classes == set(allowed_raw),
        "policy raw-input classifications must exactly match the classifications in use",
    )
    return {
        "pass": True,
        "schema": SCHEMA,
        "issue": 553,
        "entries": len(entries),
        "artifact_pins_checked": len(entries),
        "loader_pins_checked": len(entries),
        "upstream_file_pins_checked": sum(
            len(entry["raw_inputs"]["upstream_files"])
            for entry in entries
        ),
        "license_noassertion_entries": sum(
            entry["license"]["expression"] == "NOASSERTION"
            for entry in entries
        ),
        "loader_classifications": sorted(loader_classes),
        "raw_input_classifications": sorted(raw_classes),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    args = parser.parse_args(argv)
    registry_path = args.registry.resolve()
    try:
        payload = json.loads(registry_path.read_text(encoding="utf-8"))
        summary = validate_registry(payload)
    except (OSError, json.JSONDecodeError, ProvenanceError, TypeError, ValueError) as exc:
        print(f"external data provenance FAILED: {exc}")
        return 1
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
