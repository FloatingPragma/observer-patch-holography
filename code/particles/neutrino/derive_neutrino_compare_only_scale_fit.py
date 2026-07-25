#!/usr/bin/env python3
"""Fit compare-only absolute neutrino scales against the repaired branch.

Chain role: take the scale-free repaired weighted-cycle neutrino branch and
evaluate whether one positive normalization scalar can match representative
oscillation central values.

Mathematics: the repaired branch fixes the dimensionless PMNS/hierarchy
pattern and therefore fixes `Delta_hat_21 / Delta_hat_32`. A single positive
scale `lambda_nu` can rescale all absolute masses and splittings, but it
cannot change that ratio. This script records the resulting compare-only fit
surface and the central-value mismatch that remains after the repaired branch.

Declared input: the rejected target-informed weighted-cycle candidate. The
oscillation central values below are explicit comparison targets.

Output: a compare-only fit artifact for diagnostics and CLI/status reporting.
It must not feed back into theorem state or public promoted neutrino masses.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps
DEFAULT_REPAIR_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_compare_only_scale_fit.json"

# PDG 2025 neutrino review Table 14.7, Ref. [193], normal ordering,
# SK-ATM and IC24 representative central values used for compare-only fitting.
PDG_2025_NO_CENTRAL = {
    "source": "PDG 2025 neutrino review Table 14.7, Ref. [193] with SK-ATM and IC24, normal ordering, representative central values",
    "delta_m21_sq_eV2": 7.49e-5,
    "delta_m21_sq_sigma_eV2": 0.195e-5,
    "delta_m32_sq_eV2": 2.438e-3,
    "delta_m32_sq_sigma_eV2": 0.020e-3,
}


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _fit_one_scale(
    *,
    label: str,
    z: float,
    mhat: list[float],
    dh21: float,
    dh31: float,
    dh32: float,
    ref21: float,
    ref32: float,
    sig21: float,
    sig32: float,
) -> dict[str, Any]:
    lam = math.sqrt(z)
    d21 = dh21 * z
    d31 = dh31 * z
    d32 = dh32 * z
    return {
        "fit_kind": label,
        "lambda_nu": lam,
        "masses_eV": [lam * x for x in mhat],
        "delta_m_sq_eV2": {
            "21": d21,
            "31": d31,
            "32": d32,
        },
        "ratio_21_over_32": d21 / d32,
        "residual_sigma": {
            "21": (d21 - ref21) / sig21,
            "32": (d32 - ref32) / sig32,
        },
        "chi2_21_32": ((d21 - ref21) / sig21) ** 2 + ((d32 - ref32) / sig32) ** 2,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Fit compare-only absolute neutrino scales against the repaired weighted-cycle branch.")
    parser.add_argument("--repair", default=str(DEFAULT_REPAIR_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    repair = _load_json(Path(args.repair))
    if repair.get("artifact") != "oph_neutrino_weighted_cycle_repair":
        raise SystemExit("repair artifact mismatch")

    mhat = [float(x) for x in repair["scale_free_mass_normal_form"]["masses"]]
    dh = {key: float(val) for key, val in repair["scale_free_dm2_normal_form"]["dm2"].items()}
    ref21 = float(PDG_2025_NO_CENTRAL["delta_m21_sq_eV2"])
    ref32 = float(PDG_2025_NO_CENTRAL["delta_m32_sq_eV2"])
    sig21 = float(PDG_2025_NO_CENTRAL["delta_m21_sq_sigma_eV2"])
    sig32 = float(PDG_2025_NO_CENTRAL["delta_m32_sq_sigma_eV2"])

    z_solar = ref21 / dh["21"]
    z_atmos = ref32 / dh["32"]
    z_weighted = (
        dh["21"] * ref21 / (sig21**2) + dh["32"] * ref32 / (sig32**2)
    ) / (
        dh["21"] ** 2 / (sig21**2) + dh["32"] ** 2 / (sig32**2)
    )

    ratio_pred = dh["21"] / dh["32"]
    ratio_ref = ref21 / ref32

    payload = {
        "artifact": "oph_neutrino_compare_only_scale_fit",
        "generated_utc": _timestamp(),
        "status": "compare_only",
        "source_artifact": str(Path(args.repair)),
        "theorem_boundary": {
            "status": "scale_parameter_open_compare_only_fit_surface",
            "statement": (
                "The repaired weighted-cycle branch fixes the PMNS pattern and the dimensionless splitting ratio. "
                "A single positive lambda_nu rescales all masses and all Delta m^2 values together, but it cannot change the fixed ratio Delta_hat_21 / Delta_hat_32."
            ),
            "forbidden_feedback": "compare_only_fit_must_not_feed_back_into_theorem_state_or_lambda_nu_emission",
        },
        "reference_central_values": PDG_2025_NO_CENTRAL,
        "scale_free_branch": {
            "m_hat_eV": mhat,
            "delta_hat_m_sq_eV2": dh,
            "ratio_21_over_32": ratio_pred,
            "pmns_observables": repair["pmns_observables"],
        },
        "exact_central_match_possible_with_single_lambda_nu": False,
        "central_ratio_mismatch": {
            "predicted_ratio_21_over_32": ratio_pred,
            "reference_ratio_21_over_32": ratio_ref,
            "relative_difference": (ratio_pred - ratio_ref) / ratio_ref,
            "statement": (
                "One lambda_nu cannot hit both reference central splittings exactly because the repaired branch fixes a ratio that differs from the reference central ratio."
            ),
        },
        "fits": {
            "solar_only": _fit_one_scale(
                label="solar_only",
                z=z_solar,
                mhat=mhat,
                dh21=dh["21"],
                dh31=dh["31"],
                dh32=dh["32"],
                ref21=ref21,
                ref32=ref32,
                sig21=sig21,
                sig32=sig32,
            ),
            "atmospheric_only": _fit_one_scale(
                label="atmospheric_only",
                z=z_atmos,
                mhat=mhat,
                dh21=dh["21"],
                dh31=dh["31"],
                dh32=dh["32"],
                ref21=ref21,
                ref32=ref32,
                sig21=sig21,
                sig32=sig32,
            ),
            "weighted_least_squares": _fit_one_scale(
                label="weighted_least_squares",
                z=z_weighted,
                mhat=mhat,
                dh21=dh["21"],
                dh31=dh["31"],
                dh32=dh["32"],
                ref21=ref21,
                ref32=ref32,
                sig21=sig21,
                sig32=sig32,
            ),
        },
        "existing_compare_only_adapter": repair["compare_only_atmospheric_anchor"],
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
