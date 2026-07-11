#!/usr/bin/env python3
"""Seal complete IBM Runtime harvest journals into the frozen analysis packet.

This adapter is intentionally strict.  It requires the public manifest catalog
to equal the disjoint union of locked dynamic rows and locked diagnostic
calibration circuits, requires every submitted job and PUB to have one complete
harvest, converts seven-bit Qiskit counts through the analysis module's single
canonical converter, and binds every output row to journal, job, logical-
circuit, compiled-QPY, backend, and layout receipts.  Logical identity is the
preregistered normalized-OpenQASM digest; compiled QPY is retained only as an
execution artifact.  It never connects to IBM and never writes reveal material
into the resulting blinded packet.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import cayley_blind_likelihood_analysis as analysis
import record_gated_cayley_runtime as runtime


PACKET_ADAPTER_SCHEMA = "oph.runtime-analysis-packet-adapter.v1"
SUBMISSION_ROLES = ("development", "heldout")
DYNAMIC_JOINED_LAYOUTS = {
    "cayley": "final[3]|decision[1]|heated[3]",
    "mh": "after[3]|accepted[1]|before[3]",
}


class PacketAdapterError(runtime.RuntimeSafetyError):
    """A runtime receipt cannot be transformed without guessing or omission."""


def _read_bytes(path: Path, label: str) -> bytes:
    try:
        return path.read_bytes()
    except OSError as exc:
        raise PacketAdapterError(f"could not read {label}") from exc


def _submission_journal_digest(paths: Mapping[str, Path]) -> str:
    if set(paths) != set(SUBMISSION_ROLES):
        raise PacketAdapterError("submission journals must cover development and heldout")
    return runtime.sha256_json(
        {
            role: runtime.sha256_bytes(_read_bytes(paths[role], f"{role} submission journal"))
            for role in SUBMISSION_ROLES
        }
    )


def _parse_time(value: Any, label: str) -> datetime:
    if not isinstance(value, str) or not value:
        raise PacketAdapterError(f"{label} timestamp is unavailable")
    normalized = value[:-1] + "+00:00" if value.endswith("Z") else value
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError as exc:
        raise PacketAdapterError(f"{label} timestamp is invalid") from exc
    if parsed.tzinfo is None:
        raise PacketAdapterError(f"{label} timestamp must include a timezone")
    return parsed.astimezone(timezone.utc)


def _calibration_age_seconds(harvest: Mapping[str, Any]) -> float:
    metrics = harvest.get("metrics")
    if not isinstance(metrics, Mapping):
        raise PacketAdapterError("job metrics are unavailable for calibration age")
    timestamps = metrics.get("timestamps")
    if not isinstance(timestamps, Mapping):
        raise PacketAdapterError("job metrics omit timestamps")
    running = _parse_time(timestamps.get("running"), "job running")

    calibration = harvest.get("calibration_at_execution")
    properties = calibration.get("properties") if isinstance(calibration, Mapping) else None
    last_update = properties.get("last_update_date") if isinstance(properties, Mapping) else None
    if last_update is None:
        # This is conservative: the committed selection snapshot is no newer
        # than the execution calibration.  It is never replaced by wall-clock
        # time or an analyst-supplied value.
        last_update = harvest.get("properties_last_update")
    calibrated = _parse_time(last_update, "backend calibration")
    age = (running - calibrated).total_seconds()
    if not math.isfinite(age) or age < 0:
        raise PacketAdapterError("backend calibration timestamp is later than job execution")
    return age


def _validate_lock_catalog(bundle: runtime.VerifiedBundle) -> tuple[dict[str, Any], set[str]]:
    lock = bundle.analysis_lock_document
    runtime.validate_hardened_analysis_lock(bundle)
    expected_rows_raw = lock.get("expected_rows")
    calibrations = lock.get("calibrations")
    if not isinstance(expected_rows_raw, list) or not isinstance(calibrations, Mapping):
        raise PacketAdapterError("analysis lock lacks expected rows or calibrations")
    expected_rows: dict[str, Any] = {}
    for row in expected_rows_raw:
        if not isinstance(row, Mapping):
            raise PacketAdapterError("analysis lock contains a malformed dynamic row")
        opaque_id = row.get("opaque_id")
        if row.get("row_id") != opaque_id or not isinstance(opaque_id, str):
            raise PacketAdapterError("locked dynamic rows must be keyed one-to-one by opaque_id")
        if opaque_id in expected_rows:
            raise PacketAdapterError("analysis lock repeats a dynamic opaque_id")
        expected_rows[opaque_id] = row

    diagnostic_ids: list[str] = []
    for calibration in calibrations.values():
        if not isinstance(calibration, Mapping):
            raise PacketAdapterError("analysis lock contains a malformed calibration")
        control = calibration.get("control_rule")
        ids = control.get("diagnostic_opaque_ids") if isinstance(control, Mapping) else None
        if not isinstance(ids, list):
            raise PacketAdapterError("calibration control lacks diagnostic opaque IDs")
        diagnostic_ids.extend(ids)
    if len(diagnostic_ids) != len(set(diagnostic_ids)):
        raise PacketAdapterError("diagnostic opaque IDs are not globally unique")
    diagnostic_set = set(diagnostic_ids)
    dynamic_set = set(expected_rows)
    if dynamic_set & diagnostic_set:
        raise PacketAdapterError("dynamic and diagnostic catalogs overlap")
    manifest_set = {row["opaque_id"] for row in bundle.manifest["circuits"]}
    if dynamic_set | diagnostic_set != manifest_set:
        raise PacketAdapterError("manifest catalog is not exactly classified by the analysis lock")
    coverage = lock.get("catalog_coverage")
    expected_coverage = {
        "dynamic_analysis_opaque_ids": sorted(dynamic_set),
        "diagnostic_calibration_opaque_ids": sorted(diagnostic_set),
        "all_catalog_circuits_classified": True,
    }
    if coverage != expected_coverage:
        raise PacketAdapterError("analysis lock catalog_coverage differs from the manifest union")
    for opaque_id in dynamic_set:
        if bundle.circuit_families.get(opaque_id) not in {"cayley", "mh"}:
            raise PacketAdapterError("locked dynamic row maps to a non-dynamic circuit")
    for opaque_id in diagnostic_set:
        if bundle.circuit_families.get(opaque_id) != "readout_calibration":
            raise PacketAdapterError("locked diagnostic ID maps to a non-calibration circuit")
    return expected_rows, diagnostic_set


def _load_submissions(
    bundle: runtime.VerifiedBundle,
    paths: Mapping[str, Path],
) -> tuple[dict[str, dict[str, Any]], dict[str, tuple[dict[str, Any], dict[str, Any]]]]:
    registrations: dict[str, dict[str, Any]] = {}
    by_opaque: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {}
    manifest_rows = {row["opaque_id"]: row for row in bundle.manifest["circuits"]}
    for role in SUBMISSION_ROLES:
        events = runtime.read_event_journal(paths[role], bundle.manifest_sha256)
        if any(
            event.get("event_type") == "submission_failed_without_registered_job"
            for event in events
        ):
            raise PacketAdapterError("submission journal contains an unresolved failed submission")
        role_registrations = [
            event for event in events if event.get("event_type") == "job_registered"
        ]
        if not role_registrations:
            raise PacketAdapterError(f"{role} submission journal has no registered jobs")
        started = {
            event.get("group_id")
            for event in events
            if event.get("event_type") == "submission_started"
        }
        registered_groups = {event.get("group_id") for event in role_registrations}
        if started != registered_groups:
            raise PacketAdapterError("started and registered submission groups differ")
        for registration in role_registrations:
            job_id = registration.get("job_id")
            if not isinstance(job_id, str) or not job_id or job_id in registrations:
                raise PacketAdapterError("provider job IDs are missing or duplicated")
            runtime._require_hex64(
                registration.get("operator_source_sha256"),
                "submission operator_source_sha256",
            )
            runtime._require_hex64(
                registration.get("plan_sha256"), "submission plan_sha256"
            )
            if registration.get("backend_role") != role:
                raise PacketAdapterError("submission journal contains the wrong semantic role")
            opaque_ids = registration.get("opaque_ids")
            bindings = registration.get("circuit_bindings")
            if (
                not isinstance(opaque_ids, list)
                or not isinstance(bindings, list)
                or len(opaque_ids) != len(bindings)
                or len(set(opaque_ids)) != len(opaque_ids)
            ):
                raise PacketAdapterError("submission PUB bindings are malformed")
            registrations[job_id] = registration
            for opaque_id, binding in zip(opaque_ids, bindings):
                if opaque_id in by_opaque or opaque_id not in manifest_rows:
                    raise PacketAdapterError("submitted opaque IDs are duplicated or unregistered")
                if not isinstance(binding, Mapping) or binding.get("opaque_id") != opaque_id:
                    raise PacketAdapterError("submission circuit binding order is inconsistent")
                public = manifest_rows[opaque_id]
                if (
                    binding.get("logical_circuit_sha256")
                    != public["logical_circuit_sha256"]
                ):
                    raise PacketAdapterError(
                        "submitted logical-circuit digest differs from the manifest"
                    )
                if registration.get("shots") != public["shots"]:
                    raise PacketAdapterError("submitted shots differ from the manifest")
                by_opaque[opaque_id] = (registration, dict(binding))
    if set(by_opaque) != set(manifest_rows):
        raise PacketAdapterError("submission journals omit one or more manifest circuits")
    return registrations, by_opaque


def _load_harvests(
    bundle: runtime.VerifiedBundle,
    harvest_journal: Path,
    registrations: Mapping[str, Mapping[str, Any]],
) -> dict[str, tuple[dict[str, Any], dict[str, Any]]]:
    events = runtime.read_event_journal(harvest_journal, bundle.manifest_sha256)
    if any(event.get("event_type") == "job_lookup_failed" for event in events):
        raise PacketAdapterError("harvest journal contains a failed provider-job lookup")
    by_job: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for event in events:
        if event.get("event_type") == "job_harvest":
            job_id = event.get("job_id")
            if job_id not in registrations:
                raise PacketAdapterError("harvest journal contains an unregistered provider job")
            by_job[str(job_id)].append(event)
    if set(by_job) != set(registrations):
        raise PacketAdapterError("harvest journal omits one or more submitted jobs")

    results: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {}
    for job_id, registration in registrations.items():
        observations = by_job[job_id]
        completed = [event for event in observations if event.get("status") == "DONE"]
        if not completed:
            raise PacketAdapterError("a submitted provider job has no completed harvest")
        # Reharvesting is allowed, but every completed receipt must contain the
        # same raw result hashes.  Use the last append-only receipt only after
        # this equality check.
        completed_hash_sets = []
        for event in completed:
            event_results = event.get("results")
            if not isinstance(event_results, list):
                raise PacketAdapterError("completed harvest has no PUB results")
            completed_hash_sets.append(
                tuple(
                    (row.get("opaque_id"), row.get("raw_joined_counts_sha256"))
                    for row in event_results
                    if isinstance(row, Mapping)
                )
            )
        if len(set(completed_hash_sets)) != 1:
            raise PacketAdapterError("repeated completed harvests contain different counts")
        harvest = completed[-1]
        for field in (
            "operator_source_sha256",
            "plan_sha256",
            "group_id",
            "backend_name",
            "backend_role",
            "backend_role_opaque_id",
            "layout_id",
            "physical_layout",
            "properties_last_update",
            "shots",
            "opaque_ids",
            "circuit_bindings",
        ):
            if harvest.get(field) != registration.get(field):
                raise PacketAdapterError(f"harvest differs from submission field {field}")
        if harvest.get("submission_event_sha256") != registration.get("event_sha256"):
            raise PacketAdapterError("harvest is not bound to its submission event")
        event_results = harvest["results"]
        if [row.get("opaque_id") for row in event_results] != registration["opaque_ids"]:
            raise PacketAdapterError("harvest PUB order differs from submission order")
        bindings = {row["opaque_id"]: row for row in registration["circuit_bindings"]}
        for result in event_results:
            opaque_id = result["opaque_id"]
            binding = bindings[opaque_id]
            if (
                result.get("logical_circuit_sha256")
                != binding["logical_circuit_sha256"]
                or result.get("compiled_qpy_sha256") != binding["compiled_qpy_sha256"]
            ):
                raise PacketAdapterError(
                    "harvest logical/compiled circuit binding differs from submission"
                )
            joined_counts = result.get("joint_counts")
            if not isinstance(joined_counts, Mapping) or not joined_counts:
                raise PacketAdapterError("harvest PUB has no raw joined counts")
            if runtime.sha256_json(joined_counts) != result.get("raw_joined_counts_sha256"):
                raise PacketAdapterError("harvest raw joined-count digest does not verify")
            if opaque_id in results:
                raise PacketAdapterError("harvest repeats an opaque circuit result")
            results[opaque_id] = (harvest, result)
    return results


def _dynamic_rows(
    expected_rows: Mapping[str, Mapping[str, Any]],
    executions: Mapping[str, tuple[Mapping[str, Any], Mapping[str, Any]]],
) -> list[dict[str, Any]]:
    rows = []
    for expected in expected_rows.values():
        opaque_id = expected["opaque_id"]
        harvest, result = executions[opaque_id]
        family = expected.get("family")
        expected_layout = DYNAMIC_JOINED_LAYOUTS.get(family)
        if expected_layout is None:
            raise PacketAdapterError("locked dynamic row has an unsupported family")
        if harvest.get("family") != family:
            raise PacketAdapterError("runtime dynamic family differs from the lock")
        raw_counts = result["joint_counts"]
        converted = analysis.qiskit_joined_counts_to_analysis_counts(
            raw_counts,
            expected["valid_codes"],
            expected_shots=expected["shots"],
        )
        if result.get("joined_count_layout") != expected_layout:
            raise PacketAdapterError("dynamic result lacks the frozen seven-bit joined layout")
        if result.get("joined_count_width") != 7:
            raise PacketAdapterError("dynamic result joined-count width is not seven")
        identity = {
            "opaque_id": opaque_id,
            "logical_circuit_sha256": result["logical_circuit_sha256"],
            "backend_role": harvest["backend_role"],
            "backend_name": harvest["backend_name"],
            "layout_id": harvest["layout_id"],
            "physical_layout": harvest["physical_layout"],
        }
        for key, value in identity.items():
            if expected.get(key) != value:
                raise PacketAdapterError(f"runtime dynamic identity differs from lock field {key}")
        retrieved = sum(int(value) for value in raw_counts.values())
        rows.append(
            {
                "row_id": expected["row_id"],
                **identity,
                "job_id": harvest["job_id"],
                "group_id": harvest["group_id"],
                "submission_event_sha256": harvest["submission_event_sha256"],
                "harvest_event_sha256": harvest["event_sha256"],
                "raw_joined_counts_sha256": result["raw_joined_counts_sha256"],
                "calibration_age_seconds": _calibration_age_seconds(harvest),
                "declared_shots": int(expected["shots"]),
                "submitted_shots": int(harvest["shots"]),
                "retrieved_shots": retrieved,
                "excluded_shots": 0,
                "postselected": False,
                "all_outcomes_included": True,
                "counts": converted,
            }
        )
    return rows


def _calibration_results(
    lock: Mapping[str, Any],
    bundle: runtime.VerifiedBundle,
    executions: Mapping[str, tuple[Mapping[str, Any], Mapping[str, Any]]],
) -> dict[str, dict[str, Any]]:
    results = {}
    for calibration_id, calibration in lock["calibrations"].items():
        control = calibration["control_rule"]
        diagnostic_ids = list(control["diagnostic_opaque_ids"])
        counts_by_id: dict[str, Mapping[str, Any]] = {}
        provider_job_ids = []
        receipt_rows = []
        for opaque_id in diagnostic_ids:
            harvest, result = executions[opaque_id]
            descriptor = bundle.circuit_descriptors[opaque_id]
            parameters = descriptor["parameters"]
            if (
                descriptor["family"] != "readout_calibration"
                or parameters.get("evidentiary_role") != "diagnostic_only"
                or parameters.get("basis_code")
                != control["expected_basis_code_by_opaque_id"][opaque_id]
            ):
                raise PacketAdapterError("diagnostic circuit differs from its frozen basis code")
            if result.get("joined_count_width") != 4:
                raise PacketAdapterError("diagnostic result is not an exact four-bit count table")
            counts = result["joint_counts"]
            if sum(int(value) for value in counts.values()) != int(descriptor["shots"]):
                raise PacketAdapterError("diagnostic result shot count is incomplete")
            counts_by_id[opaque_id] = counts
            provider_job_ids.append(harvest["job_id"])
            receipt_rows.append(
                {
                    "opaque_id": opaque_id,
                    "job_id": harvest["job_id"],
                    "group_id": harvest["group_id"],
                    "submission_event_sha256": harvest["submission_event_sha256"],
                    "harvest_event_sha256": harvest["event_sha256"],
                    "raw_joined_counts_sha256": result["raw_joined_counts_sha256"],
                }
            )
        provider_job_ids = list(dict.fromkeys(provider_job_ids))
        receipt_sha256 = runtime.sha256_json(
            {
                "schema_version": PACKET_ADAPTER_SCHEMA,
                "calibration_id": calibration_id,
                "diagnostic_rows": receipt_rows,
            }
        )
        try:
            result = analysis.basis_calibration_control_result(
                diagnostic_counts_by_opaque_id=counts_by_id,
                control_rule=control,
                provider_job_ids=provider_job_ids,
                calibration_receipt_sha256=receipt_sha256,
            )
        except Exception as exc:
            raise PacketAdapterError("basis calibration control result could not be sealed") from exc
        results[calibration_id] = result
    return results


def seal_runtime_analysis_packet(
    *,
    bundle: runtime.VerifiedBundle,
    submission_journals: Mapping[str, Path],
    harvest_journal: Path,
    created_utc: str | None = None,
) -> dict[str, Any]:
    """Build and fully validate a blinded hardware data packet."""

    expected_rows, _ = _validate_lock_catalog(bundle)
    registrations, submitted = _load_submissions(bundle, submission_journals)
    executions = _load_harvests(bundle, harvest_journal, registrations)
    if set(executions) != set(submitted) or set(executions) != {
        row["opaque_id"] for row in bundle.manifest["circuits"]
    }:
        raise PacketAdapterError("submitted, harvested, and manifest circuit sets differ")
    rows = _dynamic_rows(expected_rows, executions)
    calibration_results = _calibration_results(
        bundle.analysis_lock_document, bundle, executions
    )
    submission_sha = _submission_journal_digest(submission_journals)
    harvest_sha = runtime.sha256_bytes(_read_bytes(harvest_journal, "harvest journal"))
    try:
        packet = analysis.seal_data_packet(
            analysis_lock_sha256=bundle.analysis_lock_sha256,
            rows=rows,
            manifest_sha256=bundle.manifest_sha256,
            submission_journal_sha256=submission_sha,
            harvest_journal_sha256=harvest_sha,
            source_kind="ibm_qpu_hardware",
            calibration_results=calibration_results,
            created_utc=created_utc,
        )
        analysis.validate_data_packet(bundle.analysis_lock_document, packet)
    except Exception as exc:
        raise PacketAdapterError("sealed Runtime data packet does not validate") from exc
    return packet


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Seal complete IBM Runtime journals into the blinded analysis packet."
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reveal", type=Path, required=True)
    parser.add_argument("--analysis-lock", type=Path, required=True)
    parser.add_argument("--development-submission-journal", type=Path, required=True)
    parser.add_argument("--heldout-submission-journal", type=Path, required=True)
    parser.add_argument("--harvest-journal", type=Path, required=True)
    parser.add_argument(
        "--created-utc",
        help=(
            "reuse the original packet creation timestamp for an exact "
            "content-hash replay; defaults to the current UTC time"
        ),
    )
    parser.add_argument("--json-out", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        # Preserve the operator's pinned Qiskit/runtime version gate before
        # rebuilding and hashing logical circuits, even though this offline
        # adapter never invokes the Runtime primitive itself.
        runtime.load_runtime_bindings()
        bundle = runtime.verify_operator_bundle(
            args.manifest,
            args.reveal,
            args.analysis_lock,
        )
        packet = seal_runtime_analysis_packet(
            bundle=bundle,
            submission_journals={
                "development": args.development_submission_journal,
                "heldout": args.heldout_submission_journal,
            },
            harvest_journal=args.harvest_journal,
            created_utc=args.created_utc,
        )
        runtime.write_new_json(args.json_out, packet)
        print(
            json.dumps(
                {
                    "data_packet_sha256": packet["data_packet_sha256"],
                    "json_out": str(args.json_out),
                    "row_count": len(packet["rows"]),
                    "source_kind": packet["source_kind"],
                },
                sort_keys=True,
            )
        )
        return 0
    except runtime.RuntimeSafetyError as exc:
        print(f"packet safety refusal: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"packet internal refusal: {type(exc).__name__}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
