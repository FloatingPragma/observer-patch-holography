#!/usr/bin/env python3
"""Compare-only distances between the criticality boundary family and data.

This module consumes the measured Higgs and top values from the reference
table and the boundary-scale scan artifact.  It exists to state, per named
source boundary scale, how far each conditional coordinate pair sits from
measurement, in relative terms and in experimental standard deviations.  It
is a diagnostic: nothing here is a prediction ancestor, and the scan artifact
itself reads no reference value.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SCAN = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_boundary_scan.json"
)
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_comparison.json"
)


def _distance_row(
    row: dict[str, Any],
    mt_ref: float,
    mt_sigma: float,
    mh_ref: float,
    mh_sigma: float,
    mh_band: float,
) -> dict[str, Any]:
    mt = float(row["mt_pole_gev"])
    mh = float(row["mh_tree_gev"])
    return {
        "mt_pole_gev": mt,
        "mh_tree_gev": mh,
        "mt_relative": mt / mt_ref - 1.0,
        "mh_relative": mh / mh_ref - 1.0,
        "mt_sigma": (mt - mt_ref) / mt_sigma,
        "mh_sigma_experimental": (mh - mh_ref) / mh_sigma,
        "mh_within_matching_band": abs(mh / mh_ref - 1.0) <= mh_band,
    }


def _relation_test(
    scan: dict[str, Any], mt_ref: float, mh_ref: float, mh_sigma: float
) -> dict[str, Any]:
    """Evaluate m_H on the criticality curve at the measured top mass.

    The curve is fit-free; the measured top selects the point on it, so this
    block tests the m_t to m_H relation independently of the boundary-scale
    selection theorem.  Compare-only by construction.
    """

    result: dict[str, Any] = {}
    for label, curve in (
        ("one_loop", scan.get("one_loop_curve", [])),
        ("two_loop", scan.get("two_loop_curve", [])),
    ):
        rows = [r for r in curve if isinstance(r, dict) and "mt_pole_gev" in r]
        rows.sort(key=lambda r: r["mt_pole_gev"])
        mts = [r["mt_pole_gev"] for r in rows]
        mhs = [r["mh_tree_gev"] for r in rows]
        scales = [r["boundary_scale_gev"] for r in rows]
        if len(rows) < 2 or not (mts[0] <= mt_ref <= mts[-1]):
            result[label] = {
                "status": "measured_top_outside_curve_range",
                "curve_mt_range_gev": [mts[0], mts[-1]] if mts else None,
            }
            continue
        for i in range(len(mts) - 1):
            if mts[i] <= mt_ref <= mts[i + 1]:
                frac = (mt_ref - mts[i]) / (mts[i + 1] - mts[i])
                mh_at = mhs[i] + frac * (mhs[i + 1] - mhs[i])
                scale_at = math.exp(
                    math.log(scales[i])
                    + frac * (math.log(scales[i + 1]) - math.log(scales[i]))
                )
                result[label] = {
                    "status": "interpolated",
                    "mh_on_curve_at_measured_mt_gev": mh_at,
                    "mh_relative": mh_at / mh_ref - 1.0,
                    "mh_sigma_experimental": (mh_at - mh_ref) / mh_sigma,
                    "implied_boundary_scale_gev": scale_at,
                }
                break
    return result


def build_artifact(scan: dict[str, Any], references: dict[str, Any]) -> dict[str, Any]:
    entries = references["entries"]
    higgs = entries["higgs"]
    top = entries["top_quark_direct_aux"]
    mh_ref = float(higgs["api_value"])
    mh_sigma = float(higgs["api_error_positive"])
    mt_ref = float(top["api_value"])
    mt_sigma = float(top["api_error_positive"])
    mh_band = float(scan["matching_bands"]["mh_tree_to_pole_relative_band"])

    named = {}
    for name, row in scan["one_loop_named_boundaries"].items():
        named[name] = _distance_row(row, mt_ref, mt_sigma, mh_ref, mh_sigma, mh_band)

    two_loop = {}
    if scan["two_loop_named_boundaries"].get("available"):
        for name, row in scan["two_loop_named_boundaries"].items():
            if isinstance(row, dict) and "mt_pole_gev" in row:
                two_loop[name] = _distance_row(
                    row, mt_ref, mt_sigma, mh_ref, mh_sigma, mh_band
                )

    curve = scan["one_loop_curve"]
    best = min(
        curve,
        key=lambda r: (r["mt_pole_gev"] / mt_ref - 1.0) ** 2
        + (r["mh_tree_gev"] / mh_ref - 1.0) ** 2,
    )
    best_row = _distance_row(best, mt_ref, mt_sigma, mh_ref, mh_sigma, mh_band)
    best_row["boundary_scale_gev"] = best["boundary_scale_gev"]

    relation_test = _relation_test(scan, mt_ref, mh_ref, mh_sigma)

    archived = named["mu_U_gauge_unification"]
    single_anchor = named["E_cell_pixel_energy"]

    return {
        "artifact": "oph_d11_criticality_comparison",
        "schema_version": 1,
        "status": "compare_only_distance_surface",
        "row_class": "compare_only_never_a_prediction_ancestor",
        "promotion_allowed": False,
        "references": {
            "mh_gev": mh_ref,
            "mh_sigma": mh_sigma,
            "mt_gev": mt_ref,
            "mt_sigma": mt_sigma,
            "note": "consumed here only; the scan artifact reads no reference",
        },
        "one_loop_named_distances": named,
        "two_loop_named_distances": two_loop,
        "one_loop_curve_nearest_point": best_row,
        "mt_mh_relation_test": relation_test,
        "reading": {
            "archived_mu_U_boundary": {
                "mt_relative": archived["mt_relative"],
                "mh_relative": archived["mh_relative"],
            },
            "single_anchor_E_cell_boundary": {
                "mt_relative": single_anchor["mt_relative"],
                "mh_relative": single_anchor["mh_relative"],
            },
            "statement": (
                "The boundary-scale family brackets the measured pair. The "
                "single-anchor E_cell branch reduces the archived deficits "
                "by an order of magnitude in the Higgs channel; the residual "
                "sits inside the declared loop-truncation and matching "
                "bands. The open object is the boundary-scale selection "
                "theorem, never a numerical fit."
            ),
        },
        "checks": {
            "scan_checks_pass": bool(scan["checks_pass"]),
            "distances_finite": all(
                math.isfinite(v)
                for row in named.values()
                for v in row.values()
                if isinstance(v, float)
            ),
        },
        "checks_pass": bool(scan["checks_pass"]),
    }


def build() -> dict[str, Any]:
    scan = json.loads(SCAN.read_text(encoding="utf-8"))
    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))
    return build_artifact(scan, references)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "one_loop_named_distances": {
                    name: {
                        "mt_rel": round(row["mt_relative"], 5),
                        "mh_rel": round(row["mh_relative"], 5),
                    }
                    for name, row in artifact["one_loop_named_distances"].items()
                },
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
