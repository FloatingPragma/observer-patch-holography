#!/usr/bin/env python3
"""Stage-3 comparison for the 2026-07-14 construction run (A7 shape).

This file is the only place in the construction run where the reference
capacities appear. It runs after the candidate certificates are emitted and
reads them as frozen inputs; the candidate evaluation cones never import this
file. References: the Lambda-located SL-4 basin 3.31e122 and the CL-3
electroweak-bridge value 3.53e122 (both in nats of capacity).
"""

from __future__ import annotations

import json
from pathlib import Path

from mpmath import mp, mpf

PRECISION = 40

REFERENCES = {
    "sl4_lambda_located": "3.31e122",
    "cl3_ew_bridge": "3.53e122",
}

FAMILY_FILES = {
    "capL": "F_candidate_capL_certificates.json",
    "capP": "F_candidate_capP_certificates.json",
    "capK": "F_candidate_capK_certificates.json",
}


def main() -> int:
    mp.dps = PRECISION
    runtime = Path(__file__).resolve().parent / "runtime"
    refs = {k: mpf(v) for k, v in REFERENCES.items()}

    rows = []
    for family, filename in FAMILY_FILES.items():
        data = json.loads((runtime / filename).read_text(encoding="utf-8"))
        for r in data["rows"]:
            entry = {"branch": r["branch"], "status": r["status"]}
            if r["status"] == "fixed_point_certified":
                fp = mpf(r["fixed_point_nats"])
                entry["fixed_point_nats"] = mp.nstr(fp, 25)
                for name, ref in refs.items():
                    entry[f"log10_ratio_to_{name}"] = mp.nstr(mp.log10(fp / ref), 8)
                entry["lands_in_either_basin"] = bool(
                    any(abs(mp.log10(fp / ref)) < 1 for ref in refs.values())
                )
            else:
                entry["lands_in_either_basin"] = False
            rows.append(entry)

    landed = [r for r in rows if r["lands_in_either_basin"]]
    certified = [r for r in rows if r["status"] == "fixed_point_certified"]
    fps = [mpf(r["fixed_point_nats"]) for r in certified]
    report = {
        "artifact": "oph_capacity_readback_construction_comparison_2026-07-14",
        "protocol": "basin-then-contract stage 3; comparison executed after all branch fixed points were recorded",
        "references_nats": REFERENCES,
        "landing_criterion": "|log10(N_fp / N_ref)| < 1 for either reference",
        "total_rows": len(rows),
        "certified_rows": len(certified),
        "landed_rows": len(landed),
        "certified_fixed_point_range_nats": {
            "lo": mp.nstr(min(fps), 20),
            "hi": mp.nstr(max(fps), 20),
        },
        "max_log10_shortfall_to_sl4_basin": mp.nstr(
            mp.log10(refs["sl4_lambda_located"] / max(fps)), 8
        ),
        "verdict": (
            "no branch lands: every certified fixed point sits at O(1)-O(10^3) nats "
            "against reference capacities of O(10^122); CL-7 stays open"
            if not landed else "landing recorded; see rows"
        ),
        "rows": rows,
    }
    out = runtime / "F_construction_comparison_2026-07-14.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {k: v for k, v in report.items() if k != "rows"}
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
