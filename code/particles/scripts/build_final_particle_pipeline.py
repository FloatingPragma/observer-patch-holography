#!/usr/bin/env python3
"""Run the simplified final particle-prediction pipeline.

This is the source-tree companion to compute_current_output_table.py.  It does
not try to regenerate every historical exploratory artifact.  It refreshes the
current constructive contracts, the public status/provenance ledgers, and the
final prediction bundle in dependency order.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
import sys


CODE_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class PipelineStep:
    id: str
    script: str
    description: str
    args: tuple[str, ...] = ()

    def command(self) -> list[str]:
        return [sys.executable, self.script, *self.args]


FINAL_PIPELINE_STEPS: tuple[PipelineStep, ...] = (
    PipelineStep(
        "p_closure_trunk",
        "P_derivation/emit_p_closure_trunk.py",
        "Refresh the compressed candidate P-trunk.",
    ),
    PipelineStep(
        "thomson_endpoint_contract",
        "P_derivation/thomson_endpoint_contract.py",
        "Refresh the Ward-projected Thomson endpoint contract.",
    ),
    PipelineStep(
        "thomson_endpoint_package",
        "P_derivation/thomson_endpoint_package.py",
        "Refresh the computed conditional Thomson endpoint package.",
    ),
    PipelineStep(
        "screening_invariant_no_go",
        "P_derivation/screening_invariant_no_go.py",
        "Refresh the screening-invariant non-identifiability certificate.",
    ),
    PipelineStep(
        "fine_structure_interval_certificate",
        "P_derivation/thomson_endpoint_interval_certificate.py",
        "Refresh the blocked fine-structure interval certificate and R_Q contract.",
    ),
    PipelineStep(
        "rg_matching_threshold_contract",
        "P_derivation/rg_matching_threshold_contract.py",
        "Refresh the RG/matching/threshold contract.",
    ),
    PipelineStep(
        "hadron_spectral_measure_contract",
        "particles/hadron/derive_ward_projected_spectral_measure_contract.py",
        "Refresh the source-backend hadronic spectral-measure contract.",
    ),
    PipelineStep(
        "source_spectral_theorem",
        "P_derivation/source_spectral_theorem.py",
        "Refresh the Ward-projected source-spectral theorem artifact.",
    ),
    PipelineStep(
        "measured_endpoint_calibration",
        "P_derivation/measured_endpoint_calibration.py",
        "Refresh the OPH plus empirical hadron closure endpoint artifact.",
    ),
    PipelineStep(
        "quark_sigma_source_nonidentifiability_obstruction",
        "particles/flavor/derive_quark_sigma_source_nonidentifiability_obstruction.py",
        "Refresh the target-free two-modulus quark-spread obstruction.",
    ),
    PipelineStep(
        "quark_axiom_level_yukawa_moduli_nonidentifiability",
        "particles/flavor/derive_quark_axiom_level_yukawa_moduli_nonidentifiability.py",
        "Refresh the no-extra-axiom MAR/Yukawa counterfamily theorem.",
    ),
    PipelineStep(
        "quark_sigma_source_boundary",
        "particles/flavor/derive_quark_sigma_source_datum_no_target_leak_required.py",
        "Project the quark-spread obstruction onto the stable promotion-gate artifact.",
    ),
    PipelineStep(
        "quark_running_mass_scheme_convention_obstruction",
        "particles/flavor/derive_quark_running_mass_scheme_convention_obstruction.py",
        "Refresh the light/heavy running-coordinate and physical-Yukawa scheme obstruction.",
    ),
    PipelineStep(
        "quark_s3_transposition_heat_shape_theorem",
        "particles/flavor/verify_s3_transposition_heat_shape.py",
        "Refresh the exact S3 transposition-Cayley spectrum and heat-gap theorem.",
    ),
    PipelineStep(
        "quark_s3_d12_template_postdiction",
        "particles/flavor/quark_s3_d12_template_postdiction.py",
        "Refresh the explicitly quarantined, dimensionless S3/D12 template diagnostic.",
        args=("--allow-template-ancestry",),
    ),
    PipelineStep(
        "quark_s3_d12_template_postdiction_audit",
        "particles/flavor/audit_quark_s3_d12_template_postdiction.py",
        "Refresh the target-separated ancestry, unit, and formula-grammar audit.",
    ),
    PipelineStep(
        "quark_rscc_module_arithmetic",
        "particles/flavor/verify_quark_rscc_module_arithmetic.py",
        "Refresh exact arithmetic and non-theorem boundaries for the RSCC ledger.",
    ),
    PipelineStep(
        "quark_rscc_completion_candidate",
        "particles/flavor/quark_rscc_completion_candidate.py",
        "Refresh the dimensionless, retrospectively selected RSCC candidate.",
        args=("--allow-retrospective-rscc",),
    ),
    PipelineStep(
        "quark_rscc_completion_candidate_audit",
        "particles/flavor/audit_quark_rscc_completion_candidate.py",
        "Refresh the RSCC ancestry, module-selection, ablation, and mixed-chart audit.",
    ),
    PipelineStep(
        "quark_further_theorem_audit",
        "particles/flavor/audit_quark_further_theorems.py",
        "Refresh the finite-MaxEnt no-go and conditional QFRC rigidity boundary.",
    ),
    PipelineStep(
        "quark_flavor_source_closure_contract",
        "particles/flavor/verify_quark_flavor_source_closure.py",
        "Refresh the exact ray lemmas and sufficient physical source-closure obligations.",
    ),
    PipelineStep(
        "quark_current_family_exact_readout_target_audit",
        "particles/flavor/derive_quark_current_family_exact_readout.py",
        "Refresh the target-anchored current-family coordinate readout.",
    ),
    PipelineStep(
        "quark_current_family_affine_anchor_target_audit",
        "particles/flavor/derive_quark_current_family_affine_anchor_theorem.py",
        "Refresh the target-audit affine coordinate identity.",
    ),
    PipelineStep(
        "quark_current_family_exact_pdg_target_audit",
        "particles/flavor/derive_quark_current_family_exact_pdg_theorem.py",
        "Refresh the mixed-convention current-family target audit.",
    ),
    PipelineStep(
        "quark_off_canonical_p_evaluator_obstruction",
        "particles/flavor/derive_quark_off_canonical_p_evaluator_obstruction.py",
        "Refresh the off-canonical quark P-evaluator no-go certificate.",
    ),
    PipelineStep(
        "quark_class_uniform_public_frame_descent_obstruction",
        "particles/flavor/derive_quark_class_uniform_public_frame_descent_obstruction.py",
        "Refresh the global public quark-frame descent no-go certificate.",
    ),
    PipelineStep(
        "quark_current_family_target_audit_completion",
        "particles/flavor/derive_quark_current_family_transport_frame_exact_pdg_completion.py",
        "Refresh the target-attached mixed-convention coordinate completion.",
    ),
    PipelineStep(
        "quark_current_family_mass_textures",
        "particles/flavor/derive_quark_current_family_transport_frame_exact_forward_yukawas.py",
        "Refresh the legacy-named GeV mass-texture audit matrices.",
    ),
    PipelineStep(
        "quark_current_family_target_audit_chain",
        "particles/flavor/derive_quark_current_family_end_to_end_exact_pdg_derivation_chain.py",
        "Refresh the target-anchored mixed-convention quark audit chain.",
    ),
    PipelineStep(
        "quark_current_family_mass_texture_wrapper",
        "particles/flavor/derive_quark_current_family_transport_frame_exact_yukawa_theorem.py",
        "Refresh the current-family mass-texture compatibility wrapper.",
    ),
    PipelineStep(
        "quark_exact_target_audit_wrapper",
        "particles/flavor/derive_quark_exact_pdg_end_to_end_theorem.py",
        "Refresh the target-anchored mixed-convention quark audit wrapper.",
    ),
    PipelineStep(
        "quark_mass_texture_end_to_end_wrapper",
        "particles/flavor/derive_quark_exact_yukawa_end_to_end_theorem.py",
        "Refresh the legacy-named end-to-end mass-texture target audit.",
    ),
    PipelineStep(
        "quark_public_sigma_descent",
        "particles/flavor/derive_quark_public_physical_sigma_datum_descent.py",
        "Refresh the selected-fiber descent without treating descent as source selection.",
    ),
    PipelineStep(
        "quark_public_mass_texture_yukawa_boundary",
        "particles/flavor/derive_quark_public_exact_yukawa_end_to_end_theorem.py",
        "Refresh the fail-closed mass-texture and physical-Yukawa wrapper.",
    ),
    PipelineStep(
        "quark_public_strengthened_sigma_frontier",
        "particles/flavor/derive_quark_public_strengthened_physical_sigma_lift_frontier.py",
        "Refresh the selected-class sigma frontier.",
    ),
    PipelineStep(
        "quark_public_exact_yukawa_promotion_frontier",
        "particles/flavor/derive_quark_public_exact_yukawa_promotion_frontier.py",
        "Refresh the public quark promotion frontier.",
    ),
    PipelineStep(
        "quark_selected_class_public_exact_evaluator",
        "particles/flavor/derive_quark_selected_class_public_exact_evaluator.py",
        "Refresh the selected-class audit evaluator.",
    ),
    PipelineStep(
        "direct_top_bridge_contract",
        "particles/calibration/derive_direct_top_bridge_contract.py",
        "Refresh the compare-only direct-top codomain obstruction after the quark audit wrappers.",
    ),
    PipelineStep(
        "quark_lane_closure_contract",
        "particles/flavor/derive_quark_lane_closure_contract.py",
        "Refresh the quark lane closure and obstruction contract.",
    ),
    PipelineStep(
        "charged_end_to_end_nonclosure",
        "particles/leptons/derive_charged_end_to_end_impossibility_theorem.py",
        "Refresh the charged-lepton end-to-end nonclosure theorem packet.",
    ),
    PipelineStep(
        "charged_trace_lift_theorem",
        "particles/leptons/derive_charged_trace_lift.py",
        "Audit the sector-isolated charged trace lift and its zero-residual gate.",
    ),
    PipelineStep(
        "derivation_gap_ledger",
        "particles/scripts/build_derivation_gap_ledger.py",
        "Refresh the coupled derivation-gap ledger.",
    ),
    PipelineStep(
        "carrier_mode_acceptance",
        "particles/scripts/build_carrier_mode_acceptance.py",
        "Refresh the fail-closed classical-carrier / quantum-particle gate.",
    ),
    PipelineStep(
        "results_status_table",
        "particles/scripts/build_results_status_table.py",
        "Refresh the current public results table.",
    ),
    PipelineStep(
        "exact_nonhadron_mass_bundle",
        "particles/scripts/build_exact_nonhadron_mass_bundle.py",
        "Refresh the exact non-hadron mass bundle.",
    ),
    PipelineStep(
        "exact_fit_surface",
        "particles/scripts/build_exact_fit_surface.py",
        "Refresh the exact-fit-only surface.",
    ),
    PipelineStep(
        "charged_mcpr_completion_conditional",
        "particles/leptons/derive_charged_mcpr_completion_conditional.py",
        "Refresh the conditional MCPR charged-response completion candidate.",
    ),
    PipelineStep(
        "charged_ward_determinant_line",
        "particles/leptons/derive_charged_ward_determinant_line.py",
        "Refresh the Ward determinant-line inversion theorem receipt.",
    ),
    PipelineStep(
        "wzh_residual_elimination_boundary",
        "particles/calibration/derive_wzh_residual_elimination_boundary.py",
        "Refresh the D10 selector rigidity and discrete repair-law boundary.",
    ),
    PipelineStep(
        "qme_source_action_rigidity_mechanism",
        "particles/calibration/derive_qme_source_action_rigidity_mechanism.py",
        "Refresh the source-action rigidity selection mechanism receipt.",
    ),
    PipelineStep(
        "cs133_feshbach_scalarization",
        "particles/clock/derive_cs133_feshbach_scalarization.py",
        "Refresh the cesium Feshbach scalarization theorem receipt.",
    ),
    PipelineStep(
        "clock_source_energy_closure_audit",
        "particles/clock/derive_clock_source_energy_closure_audit.py",
        "Refresh the fail-closed clock-to-source-energy closure audit.",
    ),
    PipelineStep(
        "d11_criticality_boundary_scan",
        "particles/calibration/derive_d11_criticality_boundary_scan.py",
        "Refresh the double-criticality boundary-scale family scan.",
    ),
    PipelineStep(
        "d11_criticality_comparison",
        "particles/calibration/derive_d11_criticality_comparison.py",
        "Refresh the compare-only criticality distance surface.",
    ),
    PipelineStep(
        "d11_boundary_scale_selection_audit",
        "particles/calibration/derive_d11_boundary_scale_selection_audit.py",
        "Refresh the boundary-scale selection audit and frozen candidate registry.",
    ),
    PipelineStep(
        "d11_boundary_scale_midpoint_selection_theorem",
        "particles/calibration/derive_d11_boundary_scale_midpoint_selection_theorem.py",
        "Refresh the conditional midpoint selection theorem receipt.",
    ),
    PipelineStep(
        "d11_anchor_reconciliation_reduction_theorems",
        "particles/calibration/derive_d11_anchor_reconciliation_reduction_theorems.py",
        "Refresh the premise-reduction theorems for the midpoint selection.",
    ),
    PipelineStep(
        "d11_carrier_census",
        "particles/calibration/derive_d11_carrier_census.py",
        "Refresh the model-level D11 carrier census and Higgs derivation trace.",
    ),
    PipelineStep(
        "down_type_register_clebsch_lane",
        "particles/flavor/derive_down_type_register_clebsch_lane.py",
        "Refresh the conditional down-type sector from the register Clebsch pattern.",
    ),
    PipelineStep(
        "up_type_register_exponent_scan",
        "particles/flavor/derive_up_type_register_exponent_scan.py",
        "Refresh the frozen up-type exponent scan and its verdict.",
    ),
    PipelineStep(
        "lambda_qcd_source_transmutation",
        "particles/qcd/derive_lambda_qcd_source_transmutation.py",
        "Refresh the Lambda_QCD dimensional-transmutation lane.",
    ),
    PipelineStep(
        "nucleon_mass_external_qcd_ratio",
        "particles/hadron/derive_nucleon_mass_external_qcd_ratio.py",
        "Refresh the conditional nucleon-mass lane.",
    ),
    PipelineStep(
        "theorem_gap_register",
        "particles/scripts/build_theorem_gap_register.py",
        "Refresh the theorem/proof gap register.",
    ),
    PipelineStep(
        "source_only_mass_prediction_surface",
        "particles/scripts/build_source_only_mass_prediction_surface.py",
        "Refresh the per-family source-only mass prediction surface.",
    ),
    PipelineStep(
        "pixel_screen_resonance_summary",
        "particles/hierarchy/verify_pixel_screen_resonance_summary.py",
        "Refresh the pixel-screen resonance summary receipt.",
        (
            "--check",
            "--output",
            "particles/hierarchy/certificates/R_pixel_screen_resonance_summary.json",
        ),
    ),
    PipelineStep(
        "pipeline_closure_status_bootstrap",
        "particles/scripts/build_particle_pipeline_closure_status.py",
        "Refresh closure status before provenance rebuild.",
    ),
    PipelineStep(
        "blind_prediction_provenance",
        "particles/scripts/build_blind_prediction_provenance.py",
        "Refresh row-level blind-prediction provenance.",
    ),
    PipelineStep(
        "pipeline_closure_status_finalize",
        "particles/scripts/build_particle_pipeline_closure_status.py",
        "Refresh closure status after provenance rebuild.",
    ),
    PipelineStep(
        "final_end_to_end_predictions",
        "particles/scripts/build_final_end_to_end_predictions.py",
        "Refresh the final current prediction bundle.",
    ),
    PipelineStep(
        "derivation_chain_closure_matrix",
        "particles/scripts/build_derivation_chain_closure_matrix.py",
        "Refresh the derivation-chain closure matrix.",
    ),
    PipelineStep(
        "mass_derivation_svg",
        "particles/scripts/generate_mass_derivation_svg.py",
        "Refresh the rendered mass-derivation graph.",
    ),
)


def build_steps() -> tuple[PipelineStep, ...]:
    return FINAL_PIPELINE_STEPS


def _run_step(step: PipelineStep, *, verbose: bool) -> None:
    completed = subprocess.run(
        step.command(),
        cwd=CODE_ROOT,
        check=False,
        text=True,
        capture_output=not verbose,
    )
    if completed.returncode != 0:
        if completed.stdout:
            sys.stdout.write(completed.stdout)
        if completed.stderr:
            sys.stderr.write(completed.stderr)
        raise subprocess.CalledProcessError(
            completed.returncode,
            step.command(),
            output=completed.stdout,
            stderr=completed.stderr,
        )


def _dry_run_payload() -> dict[str, object]:
    return {
        "artifact": "oph_final_particle_pipeline_plan",
        "cwd": str(CODE_ROOT),
        "steps": [
            {
                "id": step.id,
                "script": step.script,
                "args": list(step.args),
                "description": step.description,
            }
            for step in build_steps()
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the simplified final particle-prediction pipeline.")
    parser.add_argument("--dry-run", action="store_true", help="Print the ordered steps without executing them.")
    parser.add_argument("--verbose", action="store_true", help="Stream each builder's output.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.dry_run:
        print(json.dumps(_dry_run_payload(), indent=2, sort_keys=True))
        return 0

    for step in build_steps():
        print(f"[particle-final] {step.id}: {step.description}")
        _run_step(step, verbose=args.verbose)
    print("[particle-final] complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
