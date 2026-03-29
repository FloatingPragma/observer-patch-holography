#!/usr/bin/env python3
"""Smoke-test the full unquenched correlator artifact."""

from __future__ import annotations

import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "particles" / "runs" / "hadron" / "full_unquenched_correlator.json"


def main() -> int:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_hadron_full_unquenched_correlator":
        print("wrong full unquenched correlator artifact id", file=sys.stderr)
        return 1
    if payload.get("producer_law_id") != "oph_qcd_2p1_unquenched_measure_pushforward":
        print("full unquenched correlator should expose the predictive 2+1 measure pushforward", file=sys.stderr)
        return 1
    if payload.get("predictive_inputs_only") is not True:
        print("full unquenched correlator should stay predictive-inputs-only", file=sys.stderr)
        return 1
    if payload.get("predictive_promotion_allowed") is not False:
        print("full unquenched correlator should remain non-promoted scaffolding", file=sys.stderr)
        return 1
    channel_payloads = payload.get("channel_payloads", {})
    if {"pi_iso", "N_iso", "rho_scattering"} - set(channel_payloads):
        print("full unquenched correlator should expose pi_iso, N_iso, and rho_scattering payloads", file=sys.stderr)
        return 1
    if "current_quenched_proxy_inputs" in channel_payloads["pi_iso"]:
        print("full unquenched correlator must not feed quenched proxy inputs into the predictive producer", file=sys.stderr)
        return 1
    target_fields = set(channel_payloads["N_iso"].get("target_promoted_fields", []))
    if {"corr_direct_t", "corr_exchange_t", "corr_t"} - target_fields:
        print("N_iso payload should track full direct-minus-exchange correlators", file=sys.stderr)
        return 1
    qcd_inputs = payload.get("qcd_inputs", {})
    if {"m_u_gev", "m_d_gev", "m_s_gev"} - set(qcd_inputs):
        print("full unquenched correlator should expose predictive quark-mass inputs", file=sys.stderr)
        return 1
    population_contract = payload.get("population_contract", {})
    if population_contract.get("law_id") != "oph_qcd_2p1_lambda_ratio_ensemble_population":
        print("full unquenched correlator should expose the Lambda-ratio ensemble population contract", file=sys.stderr)
        return 1
    if population_contract.get("lambda_msbar3_rule") != "oph_qcd_lambda_msbar3_from_d10_alpha3":
        print("full unquenched correlator should name the Lambda_MSbar^(3) bridge rule", file=sys.stderr)
        return 1
    sequence_emission = population_contract.get("sequence_emission", {})
    if set(sequence_emission.get("pi_iso", [])) != {"corr_t", "am_eff_t"}:
        print("population contract should expose pi_iso sequence emission targets", file=sys.stderr)
        return 1
    if {"corr_direct_t", "corr_exchange_t", "corr_t", "am_eff_t"} - set(sequence_emission.get("N_iso", [])):
        print("population contract should expose N_iso sequence emission targets", file=sys.stderr)
        return 1
    if payload.get("qcd_inputs", {}).get("Lambda_MSbar_3_gev") is None:
        print("full unquenched correlator should consume the local Lambda_MSbar^(3) descendant when available", file=sys.stderr)
        return 1
    if payload.get("next_theorem_after_population") != "StableChannelForwardWindowConvergence":
        print("full unquenched correlator should point at stable-channel convergence as the next theorem", file=sys.stderr)
        return 1
    if payload.get("status") != "predictive_ensemble_seeded_candidate":
        print("full unquenched correlator should now expose a seeded predictive ensemble family", file=sys.stderr)
        return 1
    if payload.get("next_missing_object") != "stable_channel_sequence_population":
        print("full unquenched correlator should now point at stable-channel sequence population", file=sys.stderr)
        return 1
    family_targets = population_contract.get("family_targets", {})
    if family_targets.get("aLambda_seed") is None or family_targets.get("lambdaL_target") is None or family_targets.get("lambdaT_target") is None:
        print("population contract should populate the Lambda-unit ensemble seed targets", file=sys.stderr)
        return 1
    ensemble_family = payload.get("ensemble_family", [])
    if len(ensemble_family) < 1:
        print("full unquenched correlator should expose at least one seeded ensemble", file=sys.stderr)
        return 1
    if ensemble_family[0].get("beta") is None or ensemble_family[0].get("aLambda_msbar3") is None:
        print("seeded ensemble should carry concrete beta and aLambda values", file=sys.stderr)
        return 1
    if {"pi_iso", "N_iso"} - set((ensemble_family[0].get("sequence_targets") or {})):
        print("seeded ensemble should expose stable-channel sequence targets", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
