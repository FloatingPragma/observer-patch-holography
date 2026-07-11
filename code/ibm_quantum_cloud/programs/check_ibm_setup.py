#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import qiskit
import qiskit_aer
import qiskit_ibm_runtime
from qiskit_ibm_runtime import QiskitRuntimeService


def load_api_key(credentials_file: Path) -> str:
    text = credentials_file.read_text().strip()
    match = re.fullmatch(r"IBM cloud API key:\s*(\S+)", text)
    api_key = match.group(1) if match else text
    if not api_key or any(character.isspace() for character in api_key):
        raise ValueError(
            f"IBM cloud API key not found in {credentials_file}; expected one raw token "
            "or 'IBM cloud API key: <token>'"
        )
    return api_key


def instance_name_for_active_instance(
    active_instance_id: str | None, instances: list[dict]
) -> str | None:
    if not active_instance_id:
        return None
    for instance in instances:
        if instance.get("crn") == active_instance_id:
            return instance.get("name")
    return None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify IBM Quantum API access without submitting a job."
    )
    parser.add_argument(
        "--credentials-file",
        type=Path,
        default=Path("IBM_cloud.txt"),
        help="Path to the local credential note that contains the IBM Cloud API key.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=None,
        help="Optional path to write the verification result JSON.",
    )
    parser.add_argument(
        "--min-qubits",
        type=int,
        default=3,
        help="Minimum qubit count for the least-busy backend lookup.",
    )
    args = parser.parse_args()

    api_key = load_api_key(args.credentials_file)
    service = QiskitRuntimeService(
        channel="ibm_quantum_platform",
        token=api_key,
    )

    instances = service.instances()
    active_instance_id = service.active_instance()
    active_instance_name = instance_name_for_active_instance(active_instance_id, instances)
    usage = service.usage()
    backends = service.backends(simulator=False, operational=True)
    backend_rows = []
    for backend in sorted(backends, key=lambda item: (item.status().pending_jobs, item.name)):
        status = backend.status()
        backend_rows.append(
            {
                "name": backend.name,
                "num_qubits": int(getattr(backend.configuration(), "n_qubits", 0)),
                "pending_jobs": int(status.pending_jobs),
                "operational": bool(status.operational),
            }
        )

    least_busy = service.least_busy(
        operational=True,
        simulator=False,
        min_num_qubits=args.min_qubits,
    )

    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "verification_type": "auth_instance_usage_backend_discovery_only",
        "job_submission_test_performed": False,
        "credentials_file": str(args.credentials_file),
        "service_channel": service.channel,
        "active_instance_name": active_instance_name,
        "active_instance_plan": next(
            (item.get("plan") for item in instances if item.get("name") == active_instance_name),
            None,
        ),
        "usage_consumed_seconds": int(usage["usage_consumed_seconds"]),
        "usage_remaining_seconds": int(usage["usage_remaining_seconds"]),
        "usage_limit_seconds": int(usage["usage_limit_seconds"]),
        "usage_limit_reached": bool(usage["usage_limit_reached"]),
        "usage_period_start": usage["usage_period"]["start_time"],
        "usage_period_end": usage["usage_period"]["end_time"],
        "least_busy_backend_min_qubits": args.min_qubits,
        "least_busy_backend": least_busy.name,
        "available_operational_backends": backend_rows,
        "local_environment": {
            "python": ".".join(str(part) for part in __import__("sys").version_info[:3]),
            "qiskit": qiskit.__version__,
            "qiskit_ibm_runtime": qiskit_ibm_runtime.__version__,
            "qiskit_aer": qiskit_aer.__version__,
        },
    }

    rendered = json.dumps(result, indent=2, sort_keys=True)
    print(rendered)
    if args.json_out is not None:
        args.json_out.write_text(rendered + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
