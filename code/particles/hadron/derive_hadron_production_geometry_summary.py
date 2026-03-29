#!/usr/bin/env python3
"""Summarize the frozen hadron production geometry and storage boundary.

Chain role: translate the emitted hadron execution contract into concrete
ensemble geometry, correlator-dump shapes, and naive storage estimates.

Mathematics: deterministic bookkeeping from the fixed support schedule; no
physics values are inferred here.

OPH-derived inputs: the stable-channel cfg/source payload and the runtime
receipt that fix the seeded 2+1 ensemble family and execution schema.

Output: a machine-readable production geometry summary used to state the honest
hadron runtime boundary without pretending that production execution already
exists.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"
DEFAULT_RECEIPT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "production_geometry_summary.json"
RAW_GAUGE_BYTES_PER_SITE = 4 * 9 * 16


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the hadron production geometry summary.")
    parser.add_argument("--payload", default=str(DEFAULT_PAYLOAD))
    parser.add_argument("--receipt", default=str(DEFAULT_RECEIPT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.payload).read_text(encoding="utf-8"))
    receipt = json.loads(Path(args.receipt).read_text(encoding="utf-8"))
    schedule = list((payload.get("support_realization_schedule") or {}).get("ensemble_schedule", []))
    ensembles = []
    total_sites_all_cfg = 0
    total_raw_gauge_bytes_all_cfg = 0
    total_correlator_bytes = 0

    for entry in schedule:
        l_extent = int(entry["L"])
        t_extent = int(entry["T"])
        n_cfg = int(entry["n_cfg"])
        n_src_per_cfg = int(entry["n_src_per_cfg"])
        sites_per_cfg = l_extent**3 * t_extent
        total_sites = sites_per_cfg * n_cfg
        correlator_scalar_count = 3 * n_cfg * n_src_per_cfg * t_extent
        correlator_bytes = correlator_scalar_count * 8
        raw_bytes_per_cfg = sites_per_cfg * RAW_GAUGE_BYTES_PER_SITE
        raw_bytes_all_cfg = raw_bytes_per_cfg * n_cfg
        source_coordinates = [
            [0, 0, 0, 0],
            [l_extent // 2, l_extent // 2, l_extent // 2, t_extent // 2],
        ]
        ensembles.append(
            {
                "ensemble_id": entry["ensemble_id"],
                "family_index": entry["family_index"],
                "beta": entry["beta"],
                "L": l_extent,
                "T": t_extent,
                "aLambda_msbar3": entry["aLambda_msbar3"],
                "am_l": entry["am_l"],
                "am_s": entry["am_s"],
                "n_cfg": n_cfg,
                "n_src_per_cfg": n_src_per_cfg,
                "sites_per_cfg": sites_per_cfg,
                "total_sites_all_cfg": total_sites,
                "source_coordinates": source_coordinates,
                "correlator_arrays_required_from_backend": [
                    "pi_iso",
                    "N_iso_direct",
                    "N_iso_exchange",
                ],
                "correlator_scalar_count_backend_dump": correlator_scalar_count,
                "correlator_bytes_float64_backend_dump": correlator_bytes,
                "raw_gauge_bytes_per_cfg_naive": raw_bytes_per_cfg,
                "raw_gauge_bytes_all_cfg_naive": raw_bytes_all_cfg,
                "notes": [
                    "Correlator dump size is tiny compared with gauge-field generation/storage.",
                    "Raw gauge estimate is naive and backend-independent; actual storage may differ with compression/checkpointing.",
                ],
            }
        )
        total_sites_all_cfg += total_sites
        total_raw_gauge_bytes_all_cfg += raw_bytes_all_cfg
        total_correlator_bytes += correlator_bytes

    summary = {
        "artifact": "oph_hadron_production_geometry_summary",
        "generated_utc": _timestamp(),
        "receipt_artifact": receipt.get("artifact"),
        "manifest_artifact": "oph_hadron_production_backend_manifest",
        "storage_formula": {
            "raw_gauge_bytes_per_site_naive": RAW_GAUGE_BYTES_PER_SITE,
            "formula": "4 links/site * 9 complex/link * 16 bytes/complex = 576 bytes/site (naive full double-complex storage, no compression)",
        },
        "ensembles": ensembles,
        "totals": {
            "n_ensembles": len(ensembles),
            "total_cfg": sum(int(entry["n_cfg"]) for entry in schedule),
            "total_sites_all_cfg": total_sites_all_cfg,
            "total_raw_gauge_bytes_all_cfg_naive": total_raw_gauge_bytes_all_cfg,
            "total_correlator_bytes_float64_backend_dump": total_correlator_bytes,
        },
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
