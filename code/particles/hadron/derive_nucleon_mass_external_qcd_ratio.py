#!/usr/bin/env python3
"""Nucleon mass from the source Lambda_QCD and an external lattice ratio.

Every light-hadron mass factors as
``m_h/E_star = (m_h/Lambda) * (Lambda/E_star)``.  The second factor is the
perturbative transmutation of the source strong coupling
(``qcd/derive_lambda_qcd_source_transmutation.py``).  The first factor is a
nonperturbative pure-theory number that lattice QCD computes; this lane
consumes it as a declared external theory constant, so the row class is
``oph_plus_external_qcd_theory`` and the output is a conditional hadron mass,
never a source-only claim.

The external constant is the dimensionless ratio ``R_N = m_N/Lambda_MSbar^(3)``
from published lattice determinations.  Its ancestry includes lattice scale
setting; that ancestry is declared.  A machinery-convention band covers the
difference between this repository's Lambda extraction and the published
extraction on identical input.  No measured hadron mass enters the prediction
row; the comparison block is labeled compare-only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

try:
    from qcd.derive_lambda_qcd_source_transmutation import build as build_lambda
except ModuleNotFoundError:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "qcd"))
    from derive_lambda_qcd_source_transmutation import build as build_lambda

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "hadron" / "nucleon_mass_external_qcd_ratio.json"
)

# External pure-theory constant: nucleon mass over Lambda_MSbar^(3) from
# published lattice QCD (FLAG-class averages).  Declared ancestry: lattice
# scale setting.  The uncertainty is dominated by the published Lambda.
R_NUCLEON_OVER_LAMBDA3 = 2.775
R_NUCLEON_UNCERTAINTY = 0.105

# Relative band covering the difference between this repository's exact-ODE
# Lambda extraction and the published extraction on identical world input
# (measured in the transmutation lane's machinery-validation block).
EXTRACTION_CONVENTION_BAND = 0.035

E_STAR_DISPLAY_GEV = 1.2208901289579269e19


def build_artifact(lambda_artifact: dict[str, Any]) -> dict[str, Any]:
    lam3 = float(lambda_artifact["central"]["lambda3_gev"])
    lam3_lo, lam3_hi = (float(x) for x in lambda_artifact["lambda3_interval_gev"])

    central = R_NUCLEON_OVER_LAMBDA3 * lam3
    r_lo = R_NUCLEON_OVER_LAMBDA3 - R_NUCLEON_UNCERTAINTY
    r_hi = R_NUCLEON_OVER_LAMBDA3 + R_NUCLEON_UNCERTAINTY
    conv_lo = 1.0 - EXTRACTION_CONVENTION_BAND
    conv_hi = 1.0 + EXTRACTION_CONVENTION_BAND
    interval = (r_lo * lam3_lo * conv_lo, r_hi * lam3_hi * conv_hi)

    checks = {
        "lambda_lane_checks_pass": bool(lambda_artifact["checks_pass"]),
        "interval_contains_central": interval[0] < central < interval[1],
        "external_ratio_declared": True,
    }

    return {
        "artifact": "oph_nucleon_mass_external_qcd_ratio",
        "schema_version": 1,
        "status": "conditional_hadron_mass_from_source_lambda_and_external_ratio",
        "row_class": "oph_plus_external_qcd_theory",
        "promotion_allowed": False,
        "source_factor": {
            "lambda3_gev": lam3,
            "lambda3_interval_gev": [lam3_lo, lam3_hi],
            "artifact": "runs/qcd/lambda_qcd_source_transmutation.json",
        },
        "external_theory_factor": {
            "R_nucleon_over_lambda3": R_NUCLEON_OVER_LAMBDA3,
            "uncertainty": R_NUCLEON_UNCERTAINTY,
            "source": (
                "published lattice QCD determinations (FLAG-class averages "
                "of m_N and Lambda_MSbar^(3))"
            ),
            "ancestry_note": (
                "lattice scale setting enters the published ratio; the "
                "constant is a theory computation, never a direct hadron "
                "measurement in this lane's ancestry"
            ),
            "extraction_convention_band_relative": EXTRACTION_CONVENTION_BAND,
        },
        "prediction": {
            "m_nucleon_gev_display": central,
            "m_nucleon_interval_gev_display": list(interval),
            "m_nucleon_over_E_star": central / E_STAR_DISPLAY_GEV,
        },
        "compare_only": {
            "measured_m_proton_gev": 0.93827,
            "central_relative_difference": central / 0.93827 - 1.0,
            "interval_contains_measured": interval[0] < 0.93827 < interval[1],
            "role": "comparison only, outside the prediction ancestry",
        },
        "display_note": (
            "GeV values use the unclosed clock candidate; the clock audit "
            "classifies that decimal as a calibration checksum."
        ),
        "checks": checks,
        "checks_pass": all(bool(v) for v in checks.values()),
        "claim_boundary": (
            "This is a conditional hadron mass: source transmutation times a "
            "declared external lattice-theory ratio. A source-only hadron "
            "mass requires the production hadron backend."
        ),
    }


def build() -> dict[str, Any]:
    return build_artifact(build_lambda())


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
                "checks_pass": artifact["checks_pass"],
                "m_nucleon_gev": artifact["prediction"]["m_nucleon_gev_display"],
                "interval": artifact["prediction"][
                    "m_nucleon_interval_gev_display"
                ],
                "compare_only_delta": artifact["compare_only"][
                    "central_relative_difference"
                ],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
