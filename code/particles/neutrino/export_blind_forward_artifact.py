#!/usr/bin/env python3
"""Export the blind neutrino forward artifact before any compare surface."""

from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser(description="Export the blind neutrino forward artifact.")
    ap.add_argument("--scale-anchor", default="particles/runs/neutrino/neutrino_scale_anchor.json")
    ap.add_argument("--family", default="particles/runs/neutrino/family_response_tensor.json")
    ap.add_argument("--lift", default="particles/runs/neutrino/majorana_holonomy_lift.json")
    ap.add_argument("--pullback-metric", default="particles/runs/neutrino/majorana_phase_pullback_metric.json")
    ap.add_argument("--envelope", default="particles/runs/neutrino/majorana_phase_envelope.json")
    ap.add_argument("--majorana", default="particles/runs/neutrino/forward_majorana_matrix.json")
    ap.add_argument("--splittings", default="particles/runs/neutrino/forward_splittings.json")
    ap.add_argument("--pmns", default="particles/runs/neutrino/pmns_from_shared_basis.json")
    ap.add_argument("--out", default="particles/runs/neutrino/blind_forward_artifact.json")
    args = ap.parse_args()

    paths = {}
    arg_names = {
        "scale_anchor": "scale_anchor",
        "family": "family",
        "lift": "lift",
        "pullback_metric": "pullback_metric",
        "envelope": "envelope",
        "majorana": "majorana",
        "splittings": "splittings",
        "pmns": "pmns",
    }
    for key, arg_name in arg_names.items():
        value = pathlib.Path(getattr(args, arg_name))
        paths[key] = value if value.is_absolute() else ROOT / value
    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    majorana = load_json(paths["majorana"])
    splittings = load_json(paths["splittings"])
    payload = {
        "scale_anchor": load_json(paths["scale_anchor"]),
        "family": load_json(paths["family"]),
        "lift": load_json(paths["lift"]) if paths["lift"].exists() else None,
        "pullback_metric": load_json(paths["pullback_metric"]) if paths["pullback_metric"].exists() else None,
        "envelope": load_json(paths["envelope"]) if paths["envelope"].exists() else None,
        "majorana": majorana,
        "splittings": splittings,
        "pmns": load_json(paths["pmns"]) if paths["pmns"].exists() else None,
        "phase_mode": majorana.get("phase_mode"),
        "selector_used": majorana.get("selector_used"),
        "envelope_used": majorana.get("envelope_used"),
        "certification_status": majorana.get("certification_status"),
        "phase_certificate_source": splittings.get("phase_certificate_source"),
    }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
