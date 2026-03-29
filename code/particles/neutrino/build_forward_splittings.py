#!/usr/bin/env python3
"""Convert the forward Majorana matrix into masses, splittings, and ordering.

Chain role: expose the numerical neutrino outputs that sit directly on top of
the blind Majorana matrix and the phase certificates.

Mathematics: sorted mass extraction, `Delta m^2` computation, splitting ratios,
and ordering classification from collective-mode overlap.

OPH-derived inputs: the local forward Majorana matrix plus the phase-law or
phase-envelope certificate emitted upstream in `/particles`.

Output: the forward neutrino splitting surface consumed by the closure bundle
and table-side gap reporting.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser(description="Build blind neutrino splittings/order from the Majorana artifact.")
    ap.add_argument("--majorana", default="particles/runs/neutrino/forward_majorana_matrix.json")
    ap.add_argument("--envelope", default="particles/runs/neutrino/majorana_phase_envelope.json")
    ap.add_argument("--pullback-metric", default="particles/runs/neutrino/majorana_phase_pullback_metric.json")
    ap.add_argument("--out", default="particles/runs/neutrino/forward_splittings.json")
    args = ap.parse_args()

    majorana_path = pathlib.Path(args.majorana)
    if not majorana_path.is_absolute():
        majorana_path = ROOT / args.majorana
    envelope_path = pathlib.Path(args.envelope)
    if not envelope_path.is_absolute():
        envelope_path = ROOT / args.envelope
    pullback_metric_path = pathlib.Path(args.pullback_metric)
    if not pullback_metric_path.is_absolute():
        pullback_metric_path = ROOT / args.pullback_metric
    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    majorana = load_json(majorana_path)
    envelope = load_json(envelope_path) if envelope_path.exists() else None
    pullback_metric = load_json(pullback_metric_path) if pullback_metric_path.exists() else None
    masses = [float(value) for value in majorana.get("masses_sorted_gev", [])]
    if len(masses) != 3:
        raise ValueError("majorana artifact must provide three sorted masses")
    dm21 = (masses[1] ** 2) - (masses[0] ** 2)
    dm31 = (masses[2] ** 2) - (masses[0] ** 2)
    r = None if abs(dm31) <= 1.0e-30 else dm21 / dm31
    overlaps = [float(value) for value in majorana.get("collective_mode_overlap_by_eigenvector", [])]
    dominant_index = None if len(overlaps) != 3 else int(max(range(3), key=lambda idx: overlaps[idx]))
    if dominant_index == 2:
        ordering = "normal_like_collective_dominance"
    elif dominant_index == 0:
        ordering = "inverted_like_collective_dominance"
    else:
        ordering = "undetermined"
    selector_point_certified = bool(majorana.get("selector_point_certified", False)) or str(
        majorana.get("certification_status", "")
    ).startswith("selector_closed_")
    selector_law_certified = bool(majorana.get("selector_law_certified", False))
    if selector_law_certified:
        phase_certificate_source = str(pullback_metric_path) if pullback_metric is not None else majorana.get("inputs", {}).get("pullback_metric_artifact")
    elif selector_point_certified:
        phase_certificate_source = majorana.get("inputs", {}).get("lift_artifact")
    else:
        phase_certificate_source = str(envelope_path) if envelope is not None else None
    projector_certificate_source = (
        majorana.get("inputs", {}).get("lift_artifact") if selector_point_certified else str(envelope_path) if envelope is not None else None
    )
    if selector_law_certified or selector_point_certified:
        ordering_phase_certified = ordering
    else:
        ordering_phase_certified = ordering if envelope is not None and envelope.get("ordering_phase_stable") else None
    payload = {
        "status": "blind_forward_splittings",
        "majorana_artifact": str(majorana_path),
        "phase_certificate_source": phase_certificate_source,
        "projector_certificate_source": projector_certificate_source,
        "masses_gev_sorted": masses,
        "delta_m21_sq_gev2": dm21,
        "delta_m31_sq_gev2": dm31,
        "splitting_ratio_r": r,
        "ordering_real_seed": ordering,
        "ordering_phase_certified": ordering_phase_certified,
        "selector_point_certified": selector_point_certified,
        "selector_law_certified": selector_law_certified,
        "collective_mode_dominance": overlaps,
        "ordering_theorem_status": (
            "selector_law_certified"
            if selector_law_certified
            else "selector_point_certified"
            if selector_point_certified
            else "phase_certified"
            if ordering_phase_certified
            else "blind_real_seed_output"
        ),
        "phase_mode": majorana.get("phase_mode"),
        "certification_status": majorana.get("certification_status"),
        "phase_law_certificate_source": str(pullback_metric_path) if selector_law_certified and pullback_metric is not None else None,
        "mass_bounds_gev": None if envelope is None else envelope.get("mass_bounds_gev"),
        "delta_m21_sq_bounds_gev2": None if envelope is None else envelope.get("delta_m21_sq_bounds_gev2"),
        "delta_m31_sq_bounds_gev2": None if envelope is None else envelope.get("delta_m31_sq_bounds_gev2"),
        "collective_projector_phase_stable": None if envelope is None else envelope.get("collective_projector_phase_stable"),
        "gap_vs_radius_certificate": None if envelope is None else envelope.get("gap_vs_radius_certificate"),
        "notes": [
            "Ordering_real_seed is a report-only spectral output until the phase-envelope or selector certificate closes.",
            "The complex Majorana phase theorem remains separate from this real-first splitting surface.",
        ],
    }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
