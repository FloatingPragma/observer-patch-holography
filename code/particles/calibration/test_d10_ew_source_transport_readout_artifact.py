#!/usr/bin/env python3
"""Validate the factorized D10 electroweak readout artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FAMILY_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_observable_family.py"
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
POPULATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_population_evaluator.py"
EXACT_CLOSURE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_closure_beyond_current_carrier.py"
FIBERWISE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py"
OBSTRUCTION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_tau2_current_carrier_obstruction.py"
EXACT_WZ_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_wz_coordinate_beyond_single_tree_identity.py"
EXACT_CHART_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_mass_pair_chart_current_carrier.py"
FACTORIZATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_w_anchor_neutral_shear_factorization.py"
MINIMAL_CONDITIONAL_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_minimal_conditional_promotion.py"
TARGET_EMITTER_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_target_emitter_candidate.py"
TARGET_FREE_REPAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_target_free_repair_value_law.py"
FORWARD_TRANSMUTATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_forward_transmutation_certificate.py"
READOUT_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_readout.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"


def test_d10_source_transport_readout_uses_predictive_seed_trial() -> None:
    subprocess.run([sys.executable, str(FAMILY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POPULATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_CLOSURE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FIBERWISE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_WZ_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_CHART_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FACTORIZATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MINIMAL_CONDITIONAL_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(TARGET_EMITTER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(TARGET_FREE_REPAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_TRANSMUTATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    public = payload["public_emitted_quintet"]
    mass_pair = payload["mass_pair_predictive_candidate"]
    assert payload["artifact"] == "oph_d10_ew_source_transport_readout"
    assert public["MW_pole"] > 0.0
    assert public["MZ_pole"] > 0.0
    assert payload["predictive_mass_promotion_allowed"] is True
    assert payload["predictive_promotion_allowed"] is True
    assert payload["public_surface_candidate_allowed"] is True
    assert payload["predictive_population_closed"] is True
    assert payload["predictive_population_verdict"] == "closed_current_carrier_nonexact_running_quintet"
    assert payload["exact_closure_beyond_current_carrier_status"] == "closed"
    assert payload["exact_wz_coordinate_beyond_single_tree_identity_depends_on"] == "EWFiberwisePopulationTreeLaw_D10"
    assert payload["tau2_current_carrier_obstruction_status"] == "closed_smaller_primitive"
    assert payload["exact_mass_pair_chart_current_carrier_status"] == "closed_smaller_primitive"
    assert payload["w_anchor_neutral_shear_factorization_status"] == "closed_freeze_once_coherent_repair_law"
    assert payload["active_builder_smallest_missing_object"] is None
    assert payload["current_carrier_builder_local_frontier"] == "EWExactMassPairSelector_D10"
    assert payload["smallest_predictive_missing_object"] is None
    assert payload["broader_honest_repair_frontier"] is None
    assert payload["exact_pdg_wz_frontier"] == "EWTargetFreeRepairValueLaw_D10"
    assert payload["target_free_repair_value_law_status"] == "closed"
    assert payload["forward_transmutation_certificate_status"] == "closed_forward_p_to_t_map"
    assert payload["forward_transmutation_certificate_object_id"] == "EWForwardTransmutationCertificate_D10"
    assert payload["exact_closure_emitted_quintet"]["alpha_em_eff_inv"] > 0.0
    assert payload["selected_population_point"]["tau_2"] == 0.0
    assert payload["quartet_atomicity"]["all_four_readouts_share_one_population_point"] is True
    assert payload["quartet_atomicity"]["independent_post_population_readout_scalar_remaining"] is False
    assert abs(mass_pair["tau_2"] - (-2.311623001746158e-4)) < 1.0e-15
    assert abs(mass_pair["delta_n_dagger"] - 2.346358802434819e-4) < 1.0e-15
    assert abs(mass_pair["MW_pole"] - 80.37700001539531) < 1.0e-12
    assert abs(mass_pair["MZ_pole"] - 91.18797807794321) < 1.0e-12
    assert abs(payload["public_emitted_quintet"]["MW_pole"] - 80.37700001539531) < 1.0e-12
    assert abs(payload["public_emitted_quintet"]["MZ_pole"] - 91.18797807794321) < 1.0e-12
    assert payload["public_emitted_quintet"]["alpha_em_eff_inv"] is None
    assert payload["public_emitted_quintet"]["sin2w_eff"] is None
    assert abs(payload["public_mass_lane_quintet"]["alpha_em_eff_inv"] - 132.7524088191668) < 1.0e-12
    assert abs(payload["current_compact_emitted_quintet"]["MW_pole"] - 80.38629169244275) < 1.0e-12
    assert abs(payload["current_compact_emitted_quintet"]["MZ_pole"] - 91.18290444674243) < 1.0e-12
    assert payload["reported_readout_assignment"]["alpha_em_eff_inv"] == "ward_projected_u1q_transport_family"
    assert payload["reported_readout_assignment"]["sin2w_eff"] == "ward_projected_u1q_transport_family"
    assert payload["public_readout_split"]["electromagnetic_transport_surface"] == "WardProjectedU1QTransportLaw_D10"
    assert payload["ward_projected_transport_family"]["thomson_endpoint_alpha_em_eff_inv"] == 137.035999177
    assert payload["ward_projected_transport_family"]["transport_ratio_tQ_0_over_tQ_mZ2"] == 0.936292433907271
    assert payload["proof_gate"]["single_post_transport_tree_identity_required"] is False
    target_free_split = payload["target_free_repair_status_split"]
    assert target_free_split["status"] == "closed"
    assert target_free_split["theorem"] == "EWTargetFreeRepairValueLaw_D10"
    assert target_free_split["unconditional_source_only_status"] == "current_corpus_underdetermination_of_forward_d10_repair_law"
    assert target_free_split["minimal_conditional_principle"] == "ColorBalancedQuadraticRepairDescent_D10"
    assert target_free_split["minimal_conditional_theorem"] == "minimal_conditional_d10_forward_repair_law"
    assert target_free_split["strongest_source_only_candidate"] == "EWTargetEmitter_D10"
    assert payload["minimal_conditional_promotion_status"] == "historical_split_superseded_by_target_free_repair_theorem"
    assert payload["target_emitter_candidate_status"] == "historical_candidate_promoted_to_theorem"
    transmutation = payload["forward_transmutation_certificate"]
    assert transmutation["notation_split"]["beta_ratio_EW"]["value"] == 0.5385291530498766
    assert transmutation["notation_split"]["beta_transmutation_EW"]["value"] == 4
    assert abs(transmutation["forward_checks"]["pixel_residual"]) < 1.0e-15
