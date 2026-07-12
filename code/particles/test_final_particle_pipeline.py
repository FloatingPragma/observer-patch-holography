#!/usr/bin/env python3
"""Guard the simplified final particle-pipeline command."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_final_particle_pipeline import build_steps  # noqa: E402


def test_final_particle_pipeline_orders_contracts_before_predictions() -> None:
    steps = list(build_steps())
    ids = [step.id for step in steps]

    assert ids[0:6] == [
        "p_closure_trunk",
        "thomson_endpoint_contract",
        "thomson_endpoint_package",
        "screening_invariant_no_go",
        "fine_structure_interval_certificate",
        "rg_matching_threshold_contract",
    ]
    assert ids.index("pipeline_closure_status_bootstrap") < ids.index("blind_prediction_provenance")
    assert ids.index("blind_prediction_provenance") < ids.index("pipeline_closure_status_finalize")
    assert ids.index("pixel_screen_resonance_summary") < ids.index("pipeline_closure_status_bootstrap")
    assert ids.index("pipeline_closure_status_finalize") < ids.index("final_end_to_end_predictions")
    assert ids.index("final_end_to_end_predictions") < ids.index("derivation_chain_closure_matrix")
    assert ids.index("charged_end_to_end_nonclosure") < ids.index("charged_trace_lift_theorem")
    assert ids.index("charged_trace_lift_theorem") < ids.index("derivation_gap_ledger")
    assert ids.index("quark_sigma_source_nonidentifiability_obstruction") < ids.index(
        "quark_sigma_source_boundary"
    )
    assert ids.index("quark_axiom_level_yukawa_moduli_nonidentifiability") < ids.index(
        "quark_sigma_source_boundary"
    )
    assert ids.index("quark_running_mass_scheme_convention_obstruction") < ids.index(
        "quark_public_mass_texture_yukawa_boundary"
    )
    assert ids.index("quark_s3_transposition_heat_shape_theorem") < ids.index(
        "quark_s3_d12_template_postdiction"
    )
    assert ids.index("quark_s3_d12_template_postdiction") < ids.index(
        "quark_s3_d12_template_postdiction_audit"
    )
    assert ids.index("quark_s3_d12_template_postdiction_audit") < ids.index(
        "quark_lane_closure_contract"
    )
    assert ids.index("quark_rscc_module_arithmetic") < ids.index(
        "quark_rscc_completion_candidate"
    )
    assert ids.index("quark_rscc_completion_candidate") < ids.index(
        "quark_rscc_completion_candidate_audit"
    )
    assert ids.index("quark_rscc_completion_candidate_audit") < ids.index(
        "quark_flavor_source_closure_contract"
    )
    assert ids.index("quark_flavor_source_closure_contract") < ids.index(
        "quark_lane_closure_contract"
    )
    assert ids.index("quark_current_family_exact_readout_target_audit") < ids.index(
        "quark_current_family_affine_anchor_target_audit"
    )
    assert ids.index("quark_current_family_affine_anchor_target_audit") < ids.index(
        "quark_current_family_exact_pdg_target_audit"
    )
    assert ids.index("quark_current_family_exact_pdg_target_audit") < ids.index(
        "quark_current_family_target_audit_completion"
    )
    assert ids.index("quark_exact_target_audit_wrapper") < ids.index(
        "quark_public_mass_texture_yukawa_boundary"
    )
    assert ids.index("quark_current_family_target_audit_completion") < ids.index(
        "quark_current_family_mass_textures"
    )
    assert ids.index("quark_current_family_mass_textures") < ids.index(
        "quark_current_family_target_audit_chain"
    )
    assert ids.index("quark_public_mass_texture_yukawa_boundary") < ids.index(
        "quark_public_strengthened_sigma_frontier"
    )
    assert ids.index("quark_public_strengthened_sigma_frontier") < ids.index(
        "quark_public_exact_yukawa_promotion_frontier"
    )
    assert ids.index("quark_public_mass_texture_yukawa_boundary") < ids.index(
        "direct_top_bridge_contract"
    )
    assert ids.index("quark_lane_closure_contract") < ids.index("derivation_gap_ledger")
    assert ids[-1] == "mass_derivation_svg"
    scripts = {step.script for step in steps}
    assert "P_derivation/screening_invariant_no_go.py" in scripts
    assert "P_derivation/thomson_endpoint_interval_certificate.py" in scripts
    assert "P_derivation/measured_endpoint_calibration.py" in scripts
    assert "particles/hadron/derive_ward_projected_spectral_measure_contract.py" in scripts
    assert "particles/flavor/derive_quark_class_uniform_public_frame_descent_obstruction.py" in scripts
    assert "particles/flavor/derive_quark_sigma_source_nonidentifiability_obstruction.py" in scripts
    assert "particles/flavor/derive_quark_axiom_level_yukawa_moduli_nonidentifiability.py" in scripts
    assert "particles/flavor/derive_quark_running_mass_scheme_convention_obstruction.py" in scripts
    assert "particles/flavor/verify_s3_transposition_heat_shape.py" in scripts
    assert "particles/flavor/quark_s3_d12_template_postdiction.py" in scripts
    assert "particles/flavor/audit_quark_s3_d12_template_postdiction.py" in scripts
    assert "particles/flavor/verify_quark_flavor_source_closure.py" in scripts
    assert "particles/flavor/verify_quark_rscc_module_arithmetic.py" in scripts
    assert "particles/flavor/quark_rscc_completion_candidate.py" in scripts
    assert "particles/flavor/audit_quark_rscc_completion_candidate.py" in scripts
    assert "particles/flavor/derive_quark_current_family_exact_readout.py" in scripts
    assert "particles/flavor/derive_quark_current_family_affine_anchor_theorem.py" in scripts
    assert "particles/flavor/derive_quark_current_family_exact_pdg_theorem.py" in scripts
    assert "particles/flavor/derive_quark_current_family_transport_frame_exact_pdg_completion.py" in scripts
    assert "particles/flavor/derive_quark_current_family_transport_frame_exact_forward_yukawas.py" in scripts
    assert "particles/flavor/derive_quark_current_family_transport_frame_exact_yukawa_theorem.py" in scripts
    assert "particles/flavor/derive_quark_exact_yukawa_end_to_end_theorem.py" in scripts
    assert "particles/leptons/derive_charged_end_to_end_impossibility_theorem.py" in scripts
    assert "particles/leptons/derive_charged_trace_lift.py" in scripts
    assert "particles/hierarchy/verify_pixel_screen_resonance_summary.py" in scripts
    assert "particles/scripts/build_derivation_chain_closure_matrix.py" in scripts
    pixel_screen = steps[ids.index("pixel_screen_resonance_summary")]
    assert pixel_screen.args == (
        "--check",
        "--output",
        "particles/hierarchy/certificates/R_pixel_screen_resonance_summary.json",
    )
