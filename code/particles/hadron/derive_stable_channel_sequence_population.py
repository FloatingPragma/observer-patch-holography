#!/usr/bin/env python3
"""Define the per-ensemble stable-channel sequence families for hadrons.

Chain role: attach the `pi_iso` and `N_iso` correlator/effective-mass schema to
each seeded unquenched ensemble.

Mathematics: time-support enumeration and channel-specific observable layout in
`a*Lambda_MSbar^(3)` units.

OPH-derived inputs: the seeded unquenched ensemble family and the current
stable-channel operator/contraction choices.

Output: the stable-channel sequence population shell that awaits realized
cfg/source measurements.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FULL_UNQUENCHED = ROOT / "particles" / "runs" / "hadron" / "full_unquenched_correlator.json"
DEFAULT_GROUNDSTATE = ROOT / "particles" / "runs" / "hadron" / "stable_channel_groundstate_readout.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_population.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(full_unquenched: dict, groundstate_readout: dict | None = None) -> dict:
    groundstate_readout = groundstate_readout or {}
    family = []
    for ensemble in full_unquenched.get("ensemble_family", []):
        t_extent = int(ensemble["T"])
        a_lambda = float(ensemble["aLambda_msbar3"])
        family.append(
            {
                "ensemble_id": ensemble["ensemble_id"],
                "family_index": ensemble["family_index"],
                "beta": ensemble["beta"],
                "L": ensemble["L"],
                "T": ensemble["T"],
                "aLambda_msbar3": ensemble["aLambda_msbar3"],
                "am_l": ensemble["am_l"],
                "am_s": ensemble["am_s"],
                "t_support": list(range(t_extent)),
                "t_lambda_msbar3": [float(t_idx) * a_lambda for t_idx in range(t_extent)],
                "n_cfg": None,
                "n_src_per_cfg": None,
                "cfg_average_formula": "(1/(n_cfg*n_src_per_cfg)) * sum_{cfg,src} c_X^{cfg,src}(t)",
                "resampling": "jackknife",
                "pi_iso": {
                    "connected_only": True,
                    "valence_rule": "u=d=l",
                    "corr_formula": "C_pi^(n)(t) = sum_x <P^a(x,t) P^{a†}(0)>_{mu_n}",
                    "corr_estimator_formula": "(1/(n_cfg*n_src_per_cfg)) * sum_{cfg,src} c_pi^{cfg,src}(t)",
                    "corr_t": [],
                    "corr_t_stderr": [],
                    "am_eff_formula": "am_eff_pi^(n)(t) = log(C_pi^(n)(t) / C_pi^(n)(t+1))",
                    "am_eff_t": [],
                    "am_eff_t_stderr": [],
                    "am_ground": None,
                    "ratio_to_lambda_msbar3": None,
                    "mass_gev": None,
                    "population_status": "awaiting_measure_evaluation",
                },
                "N_iso": {
                    "contraction_plan_artifact": "N_iso.full_baryon_contractions",
                    "valence_rule": "u=d=l",
                    "corr_direct_formula": "C_N,dir^(n)(t)",
                    "corr_exchange_formula": "C_N,ex^(n)(t)",
                    "corr_formula": "C_N^(n)(t) = C_N,dir^(n)(t) - C_N,ex^(n)(t)",
                    "corr_direct_t": [],
                    "corr_direct_t_stderr": [],
                    "corr_exchange_t": [],
                    "corr_exchange_t_stderr": [],
                    "corr_t": [],
                    "corr_t_stderr": [],
                    "am_eff_formula": "am_eff_N^(n)(t) = log(|C_N^(n)(t)| / |C_N^(n)(t+1)|)",
                    "am_eff_t": [],
                    "am_eff_t_stderr": [],
                    "am_ground": None,
                    "ratio_to_lambda_msbar3": None,
                    "mass_gev": None,
                    "population_status": "awaiting_measure_evaluation",
                },
            }
        )
    return {
        "artifact": "oph_hadron_stable_channel_sequence_population",
        "generated_utc": _timestamp(),
        "sequence_emission_law_id": "oph_qcd_2p1_stable_channel_sequence_emission",
        "status": "law_closed_waiting_measure_evaluation",
        "proof_status": "law_closed_waiting_measure_evaluation",
        "predictive_promotion_allowed": False,
        "source_artifact": full_unquenched.get("artifact"),
        "population_law_id": "oph_qcd_2p1_lambda_ratio_ensemble_population",
        "operator_fixing_rule": groundstate_readout.get(
            "operator_fixing_rule",
            "lowest_local_color_singlet_by_channel_quantum_numbers",
        ),
        "readout_fixing_rule": groundstate_readout.get(
            "readout_fixing_rule",
            "full_correlator_sequence_plus_t_to_infinity_log_ratio_limit",
        ),
        "lambda_msbar_3_gev": (full_unquenched.get("qcd_inputs") or {}).get("Lambda_MSbar_3_gev"),
        "rho_l": (full_unquenched.get("dimensionless_mass_ratios") or {}).get("rho_l"),
        "rho_s": (full_unquenched.get("dimensionless_mass_ratios") or {}).get("rho_s"),
        "ensemble_sequences": family,
        "promote_after": "StableChannelForwardWindowConvergence",
        "notes": [
            "This artifact sits strictly between seeded ensemble production and stable-channel convergence.",
            "The sequence-emission law is fixed; what remains is measure evaluation of the per-ensemble correlator arrays.",
            "Ground-state masses stay unset until StableChannelForwardWindowConvergence closes.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the stable-channel sequence-population artifact.")
    parser.add_argument("--full-unquenched", default=str(DEFAULT_FULL_UNQUENCHED))
    parser.add_argument("--groundstate-readout", default=str(DEFAULT_GROUNDSTATE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    full_unquenched = json.loads(Path(args.full_unquenched).read_text(encoding="utf-8"))
    groundstate_path = Path(args.groundstate_readout)
    groundstate_readout = json.loads(groundstate_path.read_text(encoding="utf-8")) if groundstate_path.exists() else None
    artifact = build_artifact(full_unquenched, groundstate_readout)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
