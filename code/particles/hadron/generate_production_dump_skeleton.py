#!/usr/bin/env python3
"""Generate the backend correlator-dump skeleton for hadron production."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.hadron.backend_export_bundle import load_backend_input_artifact
from particles.hadron.production_execution_support import (
    build_backend_manifest,
    build_production_dump,
)


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def build_skeleton(receipt: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    task_by_id = {str(t["ensemble_id"]): t for t in (payload.get("support_realization_schedule") or {}).get("ensemble_schedule", [])}
    out: dict[str, Any] = {
        "artifact": "backend_correlator_dump.production",
        "generated_utc": receipt.get("generated_utc"),
        "production_execution": True,
        "dry_run": False,
        "surrogate_execution": False,
        "tiny_geometry_pilot": False,
        "receipt_artifact": receipt.get("artifact"),
        "manifest_artifact": "oph_hadron_production_backend_manifest",
        "ensembles": {},
        "notes": [
            "This is a shape/provenance skeleton only.",
            "Replace every null array payload with real backend-emitted values before ingestion.",
        ],
    }
    for sched in (receipt.get("execution_contract") or {}).get("ensemble_schedule", []):
        ensemble_id = str(sched["ensemble_id"])
        task = task_by_id[ensemble_id]
        t_extent = int(task["T"])
        l_extent = int(task["L"])
        sources = [
            {"src_id": "src0", "coord": [0, 0, 0, 0]},
            {"src_id": "src1", "coord": [l_extent // 2, l_extent // 2, l_extent // 2, t_extent // 2]},
        ]
        ens = {
            "ensemble_id": ensemble_id,
            "family_index": task.get("family_index"),
            "beta": task.get("beta"),
            "L": l_extent,
            "T": t_extent,
            "aLambda_msbar3": task.get("aLambda_msbar3"),
            "am_l": task.get("am_l"),
            "am_s": task.get("am_s"),
            "cfgs": {},
        }
        for cfg_id in sched.get("cfg_ids", []):
            cfg_id = str(cfg_id)
            cfg_entry = {
                "trajectory_stop": (sched.get("trajectory_stop_by_cfg") or {}).get(cfg_id),
                "sources": {},
            }
            for src in sources:
                cfg_entry["sources"][src["src_id"]] = {
                    "coord": src["coord"],
                    "expected_t_extent": t_extent,
                    "pi_iso": None,
                    "N_iso_direct": None,
                    "N_iso_exchange": None,
                }
            ens["cfgs"][cfg_id] = cfg_entry
        out["ensembles"][ensemble_id] = ens
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate a production backend correlator dump skeleton.")
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--payload", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--backend-input",
        default=None,
        help=(
            "Optional backend correlator export to normalize into the production dump schema. "
            "Accepts either the existing inline JSON format or a raw backend bundle manifest/directory."
        ),
    )
    parser.add_argument(
        "--manifest-output",
        default=None,
        help="Optional path for the normalized oph_hadron_production_backend_manifest artifact.",
    )
    args = parser.parse_args()
    receipt = _load_json(args.receipt)
    payload = _load_json(args.payload)
    if args.backend_input:
        backend_input = load_backend_input_artifact(args.backend_input)
        skeleton = build_production_dump(receipt, payload, backend_input)
        if args.manifest_output:
            manifest = build_backend_manifest(
                receipt,
                payload,
                backend_input,
                backend_input_path=args.backend_input,
            )
            Path(args.manifest_output).write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
    else:
        skeleton = build_skeleton(receipt, payload)
    Path(args.output).write_text(json.dumps(skeleton, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
