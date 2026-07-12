#!/usr/bin/env python3
"""Compare the frozen Stage-5 candidate with charged-lepton references.

This audit is deliberately downstream of the frozen candidate and may never be
used as an ancestor of the candidate builder.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CANDIDATE = ROOT / "particles" / "runs" / "leptons" / "charged_stage5_frozen_candidate.json"
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_stage5_frozen_candidate_audit.json"
CURRENT_PIXEL_BRANCH = (
    ROOT / "particles" / "runs" / "calibration" / "d10_ew_common_transport_gap_report.json"
)
ORDER = ("electron", "muon", "tau")


def _comparison_rows(
    candidate_values: dict[str, Decimal], references: dict[str, Any]
) -> tuple[dict[str, Any], list[Decimal]]:
    rows: dict[str, Any] = {}
    errors: list[Decimal] = []
    for name in ORDER:
        predicted = candidate_values[name]
        observed = Decimal(str(references["entries"][name]["value_gev"]))
        relative = predicted / observed - Decimal(1)
        errors.append(abs(relative))
        rows[name] = {
            "candidate_gev": str(predicted),
            "reference_gev": str(observed),
            "relative_error": str(relative),
            "relative_error_ppm": str(relative * Decimal(1_000_000)),
        }
    return rows, errors


def build_audit(
    candidate: dict[str, Any], references: dict[str, Any], current_pixel_branch: dict[str, Any]
) -> dict[str, Any]:
    candidate_values = {
        name: Decimal(candidate["candidate_mass_rows"][name]["value_gev"])
        for name in ORDER
    }
    rows, errors = _comparison_rows(candidate_values, references)
    frozen_v = Decimal(candidate["inputs"]["v_from_source_transmutation_gev"])
    current_v = Decimal(str(current_pixel_branch["running_family"]["v_report_gev"]))
    current_values = {
        name: value * current_v / frozen_v for name, value in candidate_values.items()
    }
    current_rows, current_errors = _comparison_rows(current_values, references)
    ratios: dict[str, Any] = {}
    for label, numerator, denominator in (
        ("muon_over_electron", "muon", "electron"),
        ("tau_over_electron", "tau", "electron"),
        ("tau_over_muon", "tau", "muon"),
    ):
        candidate_ratio = candidate_values[numerator] / candidate_values[denominator]
        reference_ratio = (
            Decimal(str(references["entries"][numerator]["value_gev"]))
            / Decimal(str(references["entries"][denominator]["value_gev"]))
        )
        ratios[label] = {
            "candidate": str(candidate_ratio),
            "reference": str(reference_ratio),
            "relative_error_ppm": str(
                (candidate_ratio / reference_ratio - Decimal(1)) * Decimal(1_000_000)
            ),
        }
    candidate_sum = sum(candidate_values.values(), Decimal(0))
    candidate_root_sum = sum((value.sqrt() for value in candidate_values.values()), Decimal(0))
    reference_values = {
        name: Decimal(str(references["entries"][name]["value_gev"])) for name in ORDER
    }
    reference_sum = sum(reference_values.values(), Decimal(0))
    reference_root_sum = sum((value.sqrt() for value in reference_values.values()), Decimal(0))
    candidate_koide = candidate_sum / candidate_root_sum**2
    reference_koide = reference_sum / reference_root_sum**2
    return {
        "artifact": "oph_charged_stage5_frozen_candidate_audit",
        "status": "COMPARE_ONLY_RETRODICTIVE_ACCURACY_AUDIT",
        "candidate_artifact": candidate["artifact"],
        "candidate_mutation_allowed": False,
        "public_prediction_promotion_allowed": False,
        "reference_packet": {
            name: references["entries"][name]["source"] for name in ORDER
        },
        "mass_rows": rows,
        "ratios": ratios,
        "koide_invariant": {
            "candidate": str(candidate_koide),
            "reference": str(reference_koide),
            "reference_relative_to_two_thirds_ppm": str(
                (reference_koide / (Decimal(2) / Decimal(3)) - Decimal(1))
                * Decimal(1_000_000)
            ),
            "interpretation": (
                "Q=2/3 is exact on the imposed balanced carrier; it is not independently derived."
            ),
        },
        "max_absolute_relative_error": str(max(errors)),
        "max_absolute_relative_error_percent": str(max(errors) * Decimal(100)),
        "pixel_branch_sensitivity": {
            "frozen_legacy_branch": {
                "P": candidate["inputs"]["frozen_pixel_P"],
                "v_gev": str(frozen_v),
                "max_absolute_relative_error_percent": str(max(errors) * Decimal(100)),
            },
            "current_public_pixel_probe": {
                "P": str(current_pixel_branch["p"]),
                "v_gev": str(current_v),
                "mass_rows": current_rows,
                "max_absolute_relative_error_percent": str(
                    max(current_errors) * Decimal(100)
                ),
                "source_artifact": (
                    "code/particles/runs/calibration/d10_ew_common_transport_gap_report.json"
                ),
            },
            "interpretation": (
                "The 0.02284 percent headline belongs to the frozen legacy P=1.63094 branch. "
                "On the current public-pixel probe the same formula remains sub-permille but is less accurate."
            ),
        },
        "mass_scheme_audit": {
            "candidate_scheme_declared": False,
            "charged_rg_threshold_map_present": False,
            "direct_pole_mass_comparison_theorem_validated": False,
            "conclusion": (
                "The numerical comparison is descriptive only until OPH specifies the candidate's "
                "renormalization coordinate and derives its QED/electroweak conversion to pole masses."
            ),
        },
        "interpretation": (
            "The frozen formula is numerically accurate, but this comparison cannot convert "
            "a historically post-hoc continuation into a prospective prediction."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, default=CANDIDATE)
    parser.add_argument("--references", type=Path, default=REFERENCE_JSON)
    parser.add_argument("--current-pixel-branch", type=Path, default=CURRENT_PIXEL_BRANCH)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_audit(
        json.loads(args.candidate.read_text(encoding="utf-8")),
        json.loads(args.references.read_text(encoding="utf-8")),
        json.loads(args.current_pixel_branch.read_text(encoding="utf-8")),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(artifact, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
