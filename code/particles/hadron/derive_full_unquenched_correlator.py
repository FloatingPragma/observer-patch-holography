#!/usr/bin/env python3
"""Seed the active 2+1 hadron ensemble family from local particle inputs.

Chain role: turn the local quark masses and QCD scale bridge into the seeded
unquenched ensemble family that every hadron artifact depends on.

Mathematics: simple lattice-unit seeding, asymptotic-scaling beta proxy, and
ensemble-family population in `Lambda_MSbar` units.

OPH-derived inputs: the local `/particles` quark masses and the QCD
`Lambda_MSbar^(3)` descendant.

Output: the seeded unquenched correlator/ensemble producer for the stable
channel and resonance lanes.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FORWARD_YUKAWAS = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
DEFAULT_CONTRACTION_PLAN = ROOT / "particles" / "runs" / "hadron" / "proton_contraction_plan.json"
DEFAULT_RHO_BASIS = ROOT / "particles" / "runs" / "hadron" / "rho_operator_basis.json"
DEFAULT_RHO_LEVELS = ROOT / "particles" / "runs" / "hadron" / "rho_levels.json"
DEFAULT_QCD_LAMBDA = ROOT / "particles" / "runs" / "qcd" / "lambda_msbar_descendant.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "full_unquenched_correlator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _quark_inputs(forward_yukawas: dict[str, object]) -> dict[str, float]:
    singular_values_u = [float(value) for value in forward_yukawas.get("singular_values_u") or []]
    singular_values_d = [float(value) for value in forward_yukawas.get("singular_values_d") or []]
    if len(singular_values_u) < 1 or len(singular_values_d) < 2:
        raise ValueError("forward_yukawas must provide singular_values_u and singular_values_d")
    m_u = singular_values_u[0]
    m_d = singular_values_d[0]
    m_s = singular_values_d[1]
    m_l = 0.5 * (m_u + m_d)
    return {
        "m_u_gev": m_u,
        "m_d_gev": m_d,
        "m_s_gev": m_s,
        "m_l_gev": m_l,
    }


def _beta_from_a_lambda(a_lambda_msbar3: float) -> float:
    if a_lambda_msbar3 <= 0.0:
        raise ValueError("aLambda_msbar3 must be positive")
    # Simple 1-loop n_f=3 asymptotic-scaling proxy in Lambda-units.
    return 6.0 + (9.0 / (2.0 * math.pi * math.pi)) * math.log(1.0 / a_lambda_msbar3)


def _seed_ensemble_family(rho_l: float | None, rho_s: float | None) -> tuple[dict[str, object] | None, list[dict[str, object]]]:
    if rho_l is None or rho_s is None or rho_l <= 0.0 or rho_s <= 0.0:
        return None, []

    a_lambda_seed = math.sqrt(rho_l * rho_s)
    lambda_l_target = 8.0
    lambda_t_target = 16.0
    family: list[dict[str, object]] = []
    for n in range(3):
        a_lambda_n = a_lambda_seed * (0.5 ** n)
        family.append(
            {
                "ensemble_id": f"qcd_2p1_seed_n{n}",
                "family_index": n,
                "beta": _beta_from_a_lambda(a_lambda_n),
                "L": math.ceil(lambda_l_target / a_lambda_n),
                "T": math.ceil(lambda_t_target / a_lambda_n),
                "aLambda_msbar3": a_lambda_n,
                "am_l": a_lambda_n * rho_l,
                "am_s": a_lambda_n * rho_s,
                "population_status": "seeded_predictive_unquenched_candidate",
                "stable_channel_sequence_population_status": "awaiting_sequence_population",
                "sequence_targets": {
                    "pi_iso": ["corr_t", "am_eff_t"],
                    "N_iso": ["corr_direct_t", "corr_exchange_t", "corr_t", "am_eff_t"],
                },
            }
        )
    return (
        {
            "aLambda_seed": a_lambda_seed,
            "aLambda_seed_formula": "sqrt(rho_l * rho_s)",
            "beta_formula": "6 + (9 / (2*pi^2)) * log(1 / aLambda_msbar3_n)",
            "lambdaL_target": lambda_l_target,
            "lambdaL_target_formula": "8",
            "lambdaT_target": lambda_t_target,
            "lambdaT_target_formula": "16",
            "sea_flavors": "2+1",
            "qed_off": True,
        },
        family,
    )


def build_artifact(
    forward_yukawas: dict[str, object],
    contraction_plan: dict[str, object],
    rho_basis: dict[str, object],
    rho_levels: dict[str, object],
    lambda_msbar_3_gev: float | None,
    qcd_lambda_artifact: dict[str, object] | None,
) -> dict[str, object]:
    qcd_inputs = _quark_inputs(forward_yukawas)
    rho_l = None
    rho_s = None
    if lambda_msbar_3_gev is not None and lambda_msbar_3_gev > 0.0:
        rho_l = qcd_inputs["m_l_gev"] / lambda_msbar_3_gev
        rho_s = qcd_inputs["m_s_gev"] / lambda_msbar_3_gev
    family_targets, seeded_ensemble_family = _seed_ensemble_family(rho_l, rho_s)
    population_seeded = family_targets is not None
    rho_irreps = rho_basis.get("irreps") or rho_levels.get("irreps") or []
    rho_operator_ids = rho_basis.get("operator_ids") or []
    seed_ensemble_id = seeded_ensemble_family[0]["ensemble_id"] if seeded_ensemble_family else "qcd_2p1_declared_seed"

    return {
        "artifact": "oph_hadron_full_unquenched_correlator",
        "generated_utc": _timestamp(),
        "status": (
            "predictive_ensemble_seeded_candidate"
            if population_seeded
            else "predictive_input_family_declared"
        ),
        "proof_status": "candidate_only",
        "predictive_promotion_allowed": False,
        "predictive_inputs_only": True,
        "producer_law_id": "oph_qcd_2p1_unquenched_measure_pushforward",
        "theorem_scope": "stable_channels_plus_rho_scattering_inputs",
        "next_missing_object": (
            "stable_channel_sequence_population"
            if population_seeded
            else (
                "predictive_unquenched_ensemble_seed"
                if lambda_msbar_3_gev is not None
                else "Lambda_MSbar_3_gev"
            )
        ),
        "next_theorem_after_population": "StableChannelForwardWindowConvergence",
        "readout_consumers": {
            "stable_channel_groundstate_readout": "oph_hadron_stable_channel_groundstate_readout",
            "rho_scattering_levels": rho_levels.get("artifact"),
        },
        "sea_flavor_target": "2+1_unquenched",
        "qed_off": True,
        "light_mass_rule": "m_l=(m_u+m_d)/2",
        "source_artifacts": {
            "forward_yukawas": forward_yukawas.get("artifact"),
            "proton_contraction_plan": contraction_plan.get("artifact"),
            "rho_operator_basis": rho_basis.get("artifact"),
            "rho_levels": rho_levels.get("artifact"),
            "qcd_lambda_msbar3": None if qcd_lambda_artifact is None else qcd_lambda_artifact.get("artifact"),
        },
        "qcd_inputs": {
            "Lambda_MSbar_3_gev": lambda_msbar_3_gev,
            "m_u_gev": qcd_inputs["m_u_gev"],
            "m_d_gev": qcd_inputs["m_d_gev"],
            "m_s_gev": qcd_inputs["m_s_gev"],
            "m_l_gev": qcd_inputs["m_l_gev"],
        },
        "predictive_input_status": {
            "Lambda_MSbar_3_gev": (
                "consumed_from_secondary_qcd_descendant"
                if qcd_lambda_artifact is not None
                else (
                    "provided_current_branch_candidate"
                    if lambda_msbar_3_gev is not None
                    else "awaiting_secondary_qcd_lambda_msbar_descendant"
                )
            ),
            "quark_masses": "consumed_from_reference_free_forward_yukawas_candidate",
        },
        "dimensionless_mass_ratios": {
            "rho_l": rho_l,
            "rho_s": rho_s,
            "rho_l_formula": "(m_u_gev + m_d_gev) / (2 * Lambda_MSbar_3_gev)",
            "rho_s_formula": "m_s_gev / Lambda_MSbar_3_gev",
        },
        "population_contract": {
            "law_id": "oph_qcd_2p1_lambda_ratio_ensemble_population",
            "lambda_msbar3_source_artifact": "oph_qcd_lambda_msbar3",
            "lambda_msbar3_rule": "oph_qcd_lambda_msbar3_from_d10_alpha3",
            "quark_source_artifact": forward_yukawas.get("artifact"),
            "status": (
                "predictive_ensemble_seeded_candidate"
                if population_seeded
                else (
                    "lambda_ratio_family_ready_for_population"
                    if rho_l is not None and rho_s is not None
                    else "awaiting_lambda_msbar3_descendant"
                )
            ),
            "seed_status": (
                "populated_from_lambda_ratio_law"
                if population_seeded
                else "awaiting_lambda_msbar3_descendant"
            ),
            "qcd_descendants": {
                "m_u_gev": "singular_values_u[0]",
                "m_d_gev": "singular_values_d[0]",
                "m_s_gev": "singular_values_d[1]",
                "m_l_gev": "(m_u_gev + m_d_gev) / 2",
                "rho_l": "m_l_gev / Lambda_MSbar_3_gev",
                "rho_s": "m_s_gev / Lambda_MSbar_3_gev",
            },
            "ensemble_family_law": {
                "family_index": "n = 0,1,2,...",
                "aLambda_msbar3_n": "aLambda_seed * 2**(-n)",
                "beta_n": "oph_beta_from_aLambda_msbar3(aLambda_msbar3_n)",
                "L_n": "ceil(lambdaL_target / aLambda_msbar3_n)",
                "T_n": "ceil(lambdaT_target / aLambda_msbar3_n)",
                "am_l_n": "aLambda_msbar3_n * rho_l",
                "am_s_n": "aLambda_msbar3_n * rho_s",
            },
            "family_targets": {
                "aLambda_seed": None if family_targets is None else family_targets["aLambda_seed"],
                "aLambda_seed_formula": None if family_targets is None else family_targets["aLambda_seed_formula"],
                "beta_formula": None if family_targets is None else family_targets["beta_formula"],
                "lambdaL_target": None if family_targets is None else family_targets["lambdaL_target"],
                "lambdaL_target_formula": None if family_targets is None else family_targets["lambdaL_target_formula"],
                "lambdaT_target": None if family_targets is None else family_targets["lambdaT_target"],
                "lambdaT_target_formula": None if family_targets is None else family_targets["lambdaT_target_formula"],
                "sea_flavors": "2+1",
                "qed_off": True,
            },
            "measure_formula": (
                "dmu_n(U)=Z_n^{-1} exp(-S_g(U;beta_n)) det D_l(U;am_l_n)^2 det D_s(U;am_s_n) dU"
            ),
            "sequence_emission": {
                "pi_iso": ["corr_t", "am_eff_t"],
                "N_iso": ["corr_direct_t", "corr_exchange_t", "corr_t", "am_eff_t"],
            },
            "per_ensemble_sequence_contract": {
                "pi_iso": {
                    "corr_formula": "C_pi^(n)(t) = sum_x <P^a(x,t) P^{a\\dagger}(0)>_{mu_n}",
                    "am_eff_formula": "am_eff_pi^(n)(t) = log(C_pi^(n)(t) / C_pi^(n)(t+1))",
                },
                "N_iso": {
                    "corr_direct_formula": "C_N,dir^(n)(t)",
                    "corr_exchange_formula": "C_N,ex^(n)(t)",
                    "corr_formula": "C_N^(n)(t) = C_N,dir^(n)(t) - C_N,ex^(n)(t)",
                    "am_eff_formula": "am_eff_N^(n)(t) = log(|C_N^(n)(t)| / |C_N^(n)(t+1)|)",
                },
            },
        },
        "measure_formula": (
            "dmu_e(U)=Z_e^{-1} exp(-S_g(U;beta_e)) det D_l(U;am_l^(e))^2 det D_s(U;am_s^(e)) dU"
        ),
        "ensemble_family": [
            *seeded_ensemble_family
        ] if seeded_ensemble_family else [
            {
                "ensemble_id": "qcd_2p1_declared_seed",
                "beta": None,
                "L": None,
                "T": None,
                "aLambda_msbar3": None,
                "am_l": None,
                "am_s": None,
                "am_l_formula": "aLambda_msbar3 * rho_l",
                "am_s_formula": "aLambda_msbar3 * rho_s",
                "population_status": "awaiting_predictive_unquenched_ensemble_choice",
            }
        ],
        "channel_payloads": {
            "pi_iso": {
                "kind": "zero_momentum_two_point",
                "operator_id": "P_local_isovector",
                "operator_formula": "sum_x qbar gamma5 tau^a q",
                "correlator_formula": "C_pi^(e)(t) = sum_x <P^a(x,t) P^{a\\dagger}(0)>_{mu_e}",
                "ensemble_id": seed_ensemble_id,
                "corr_t": [],
                "am_eff_t": [],
                "am_ground": None,
                "ratio_to_lambda_msbar3": None,
                "mass_gev": None,
                "effective_mass_formula": "am_eff_pi(t)=log(C_pi(t)/C_pi(t+1))",
                "mass_rule": "mass_gev = ratio_to_lambda_msbar3 * Lambda_MSbar_3_gev",
                "target_promoted_fields": [
                    "corr_t",
                    "am_eff_t",
                    "am_ground",
                    "ratio_to_lambda_msbar3",
                    "mass_gev",
                ],
            },
            "N_iso": {
                "kind": "direct_minus_exchange_baryon_two_point",
                "operator_id": "N_local_Cgamma5",
                "operator_formula": "sum_x eps^{abc}(u^T C gamma5 d)u",
                "parity_projector": "(1+gamma_0)/2",
                "correlator_formula": "C_N^(e)(t) = C_direct^(e)(t) - C_exchange^(e)(t)",
                "contraction_rules": contraction_plan.get("contraction_rules"),
                "ensemble_id": seed_ensemble_id,
                "corr_direct_t": [],
                "corr_exchange_t": [],
                "corr_t": [],
                "am_eff_t": [],
                "am_ground": None,
                "ratio_to_lambda_msbar3": None,
                "mass_gev": None,
                "effective_mass_formula": "am_eff_N(t)=log(|C_N(t)|/|C_N(t+1)|)",
                "mass_rule": "mass_gev = ratio_to_lambda_msbar3 * Lambda_MSbar_3_gev",
                "target_promoted_fields": [
                    "corr_direct_t",
                    "corr_exchange_t",
                    "corr_t",
                    "am_eff_t",
                    "am_ground",
                    "ratio_to_lambda_msbar3",
                    "mass_gev",
                ],
            },
            "rho_scattering": {
                "kind": "correlation_matrix_family",
                "elastic_channel": ((rho_basis.get("channel") or {}).get("elastic_channel")),
                "operator_ids": rho_operator_ids,
                "irreps": rho_irreps,
                "matrix_formula": "C_ij^(e,d,Lambda)(t) = <O_i^(d,Lambda)(t) O_j^(d,Lambda)dagger(0)>_{mu_e}",
                "correlation_matrices": {},
                "principal_correlators": {},
                "aE_lab": [],
                "aE_cm": [],
                "ak": [],
                "q": [],
                "delta1_rad": [],
                "target_promoted_fields": [
                    "correlation_matrices",
                    "principal_correlators",
                    "aE_lab",
                    "aE_cm",
                    "ak",
                    "q",
                    "delta1_rad",
                ],
                "existing_scaffolding_inputs": {
                    "rho_operator_basis_artifact": rho_basis.get("artifact"),
                    "rho_levels_artifact": rho_levels.get("artifact"),
                    "level_point_schema": rho_levels.get("level_point_schema"),
                },
            },
        },
        "missing_physics": [
            "stable_channel_sequence_population" if population_seeded else "predictive_unquenched_ensemble_seed",
            "stable_channel_forward_window_convergence",
            "finite_volume_control",
            "rho_finite_volume_spectrum_extraction",
        ],
        "notes": [
            "This artifact is an upstream 2+1-flavor QCD measure pushforward that feeds the hadron readouts.",
            "It does not promote hadron masses; it defines the unquenched correlator payloads that must become real before the stable and rho readouts can be predictive.",
            "The producer-side seed is derived in Lambda-units from the OPH quark descendants and the local Lambda_MSbar^(3) descendant.",
            "The next stable-channel gate after seeding is sequence population and then convergence, not another schema change.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the full unquenched correlator artifact.")
    parser.add_argument("--forward-yukawas", default=str(DEFAULT_FORWARD_YUKAWAS))
    parser.add_argument("--contraction-plan", default=str(DEFAULT_CONTRACTION_PLAN))
    parser.add_argument("--rho-basis", default=str(DEFAULT_RHO_BASIS))
    parser.add_argument("--rho-levels", default=str(DEFAULT_RHO_LEVELS))
    parser.add_argument("--qcd-lambda", default=str(DEFAULT_QCD_LAMBDA))
    parser.add_argument("--lambda-msbar-gev", type=float, default=None)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    forward_yukawas = _load(Path(args.forward_yukawas))
    contraction_plan = _load(Path(args.contraction_plan))
    rho_basis = _load(Path(args.rho_basis))
    rho_levels = _load(Path(args.rho_levels))
    qcd_lambda_path = Path(args.qcd_lambda)
    qcd_lambda_artifact = _load(qcd_lambda_path) if qcd_lambda_path.exists() else None
    lambda_msbar_3_gev = args.lambda_msbar_gev
    if lambda_msbar_3_gev is None and qcd_lambda_artifact is not None:
        lambda_msbar_3_gev = float(qcd_lambda_artifact["outputs"]["Lambda_MSbar_3_gev"])
    artifact = build_artifact(
        forward_yukawas,
        contraction_plan,
        rho_basis,
        rho_levels,
        lambda_msbar_3_gev,
        qcd_lambda_artifact,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
