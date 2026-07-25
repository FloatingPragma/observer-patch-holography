#!/usr/bin/env python3
"""Convert the forward Majorana matrix into ascending singular values and gaps.

Chain role: expose the numerical neutrino outputs that sit directly on top of
the blind Majorana matrix and the phase certificates.

Mathematics: sorted singular-value extraction and ascending squared-gap
reporting. Collective-mode overlap is retained as a spectral
diagnostic; it does not select the physical solar pair or atmospheric sign.

Declared pipeline inputs: the local forward Majorana matrix plus the phase-law
or phase-envelope certificate emitted upstream in `/particles`. This spectral
step does not promote source-open inputs.

Output: the forward neutrino splitting surface consumed by the closure bundle
and table-side gap reporting.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Build ascending neutrino singular-value gaps from the Majorana artifact; physical ordering remains open."
    )
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
    ascending_gaps = {
        "s1_minus_s0": (masses[1] ** 2) - (masses[0] ** 2),
        "s2_minus_s0": (masses[2] ** 2) - (masses[0] ** 2),
        "s2_minus_s1": (masses[2] ** 2) - (masses[1] ** 2),
    }
    overlaps = [float(value) for value in majorana.get("collective_mode_overlap_by_eigenvector", [])]
    dominant_index = None if len(overlaps) != 3 else int(max(range(3), key=lambda idx: overlaps[idx]))
    collective_mode_location = (
        "s2" if dominant_index == 2 else "s0" if dominant_index == 0 else "s1" if dominant_index == 1 else None
    )
    physical_assignment = None
    dm21 = dm31 = dm32 = None
    r = None
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
    # Phase or selector certification can stabilize a spectral subspace, but
    # it cannot by itself attach the physical nu_i labels.
    ordering_phase_certified = None
    source_closure_status = dict(majorana.get("source_closure_status") or {"closed": False})
    source_closed = source_closure_status.get("closed") is True
    payload = {
        "artifact": "oph_neutrino_forward_splittings",
        "status": "blind_forward_splittings",
        "proof_scope": "ascending_spectral_algebra_conditional_on_declared_majorana_matrix",
        "source_only_physical_input_eligible": source_closed,
        "public_surface_candidate_allowed": False,
        "source_closure_status": source_closure_status,
        "majorana_artifact": str(majorana_path),
        "phase_certificate_source": phase_certificate_source,
        "projector_certificate_source": projector_certificate_source,
        "masses_gev_sorted": masses,
        "ascending_state_labels": ["s0", "s1", "s2"],
        "ascending_mass_sq_gaps_gev2": ascending_gaps,
        "ascending_gap_ratio_s10_over_s20": (
            None
            if abs(ascending_gaps["s2_minus_s0"]) <= 1.0e-30
            else ascending_gaps["s1_minus_s0"] / ascending_gaps["s2_minus_s0"]
        ),
        "mass_eigenstate_label_status": "ascending_singular_states_only",
        "missing_mass_label_object": "source_derived_solar_pair_and_atmospheric_sign_rule",
        "physical_mass_label_assignment": physical_assignment,
        "delta_m21_sq_gev2": dm21,
        "delta_m31_sq_gev2": dm31,
        "delta_m32_sq_gev2": dm32,
        "splitting_ratio_r": r,
        "ordering_real_seed": "unresolved_without_mass_eigenstate_label_rule",
        "collective_mode_location": collective_mode_location,
        "ordering_phase_certified": ordering_phase_certified,
        "selector_point_certified": selector_point_certified,
        "selector_law_certified": selector_law_certified,
        "collective_mode_dominance": overlaps,
        "ordering_theorem_status": "not_established_mass_label_rule_absent",
        "phase_mode": majorana.get("phase_mode"),
        "certification_status": majorana.get("certification_status"),
        "phase_law_certificate_source": str(pullback_metric_path) if selector_law_certified and pullback_metric is not None else None,
        "mass_bounds_gev": None if envelope is None else envelope.get("mass_bounds_gev"),
        "delta_m21_sq_bounds_gev2": None,
        "delta_m31_sq_bounds_gev2": None,
        "ascending_mass_sq_gap_bounds_gev2": None
        if envelope is None
        else {
            "s1_minus_s0": envelope.get("ascending_gap_s10_sq_bounds_gev2"),
            "s2_minus_s0": envelope.get("ascending_gap_s20_sq_bounds_gev2"),
        },
        "collective_projector_phase_stable": None if envelope is None else envelope.get("collective_projector_phase_stable"),
        "gap_vs_radius_certificate": None if envelope is None else envelope.get("gap_vs_radius_certificate"),
        "notes": [
            "Ascending singular-value sorting does not select normal or inverted physical mass labels.",
            "Collective-mode dominance is a spectral diagnostic and is not an ordering theorem.",
            "Selector and phase certificates stabilize the declared spectral construction but do not supply the missing solar-pair and atmospheric-sign label rule.",
            "The complex Majorana phase theorem remains separate from this real-first splitting surface.",
        ],
    }
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
