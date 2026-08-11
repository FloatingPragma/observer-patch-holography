#!/usr/bin/env python3
"""Run the documented mandatory scientific suite from a clean checkout (#507).

This is the one command REPRODUCE.md documents and CI enforces:

    python tools/run_mandatory_suite.py                # mandatory suite
    python tools/run_mandatory_suite.py --certificates # + exact certificate suites
    python tools/run_mandatory_suite.py --certificate-smoke-only

The mandatory suite validates the claim registry against its live gates,
validates external inputs, public quantitative surfaces, null-model controls,
and the paper release manifest, then proves those gates reject isolated
false-green mutations. It also proves the full scientific collection imports
cleanly (which is what keeps the optional cloud/hardware lanes fail-closed)
and executes the fast fixture suites. The exact certificate suites (#566
port-current, #314 matter-lift) run in their own CI workflow on their own
triggers; `--certificates` runs them here with the same commands.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MANDATORY_STEPS: list[tuple[str, list[str]]] = [
    (
        "Validate the committed open-problem ledger offline",
        [sys.executable, "tools/build_open_problem_ledger.py", "--check"],
    ),
    ("Validate claim registry", [sys.executable, "tools/check_claim_registry.py"]),
    (
        "Validate the selection ledger and its generated surface",
        [sys.executable, "tools/build_selection_ledger.py", "--check"],
    ),
    (
        "Validate the frozen-prediction ladder and its generated surface",
        [sys.executable, "tools/build_fz_registry.py", "--check"],
    ),
    (
        "Validate the comparison-value-free quantum-carrier status packet",
        [
            sys.executable,
            "code/particles/scripts/build_quantum_carrier_status.py",
            "--validate-only",
        ],
    ),
    (
        "Replay the quantum-carrier status independently",
        [
            sys.executable,
            "code/particles/verify_quantum_carrier_status_independent.py",
        ],
    ),
    (
        "Execute the quantum-carrier status mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/particles/test_quantum_carrier_status.py",
        ],
    ),
    (
        "Validate the deterministic postdiction ledger",
        [
            sys.executable,
            "code/particles/scripts/build_postdiction_ledger.py",
            "--check",
        ],
    ),
    (
        "Replay the B12 common-reference obstruction independently",
        [
            sys.executable,
            "code/thermodynamics/common_reference_obstruction/verify_common_reference_obstruction.py",
        ],
    ),
    (
        "Execute the B12 common-reference mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/thermodynamics/common_reference_obstruction/test_common_reference_obstruction.py",
        ],
    ),
    (
        "Replay the B13 algebraic phase-lift boundary hermetically",
        [
            sys.executable,
            "code/born_context_phase_lift/verify_source_phase_lift.py",
        ],
    ),
    (
        "Execute the B13 phase-lift provenance and mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/born_context_phase_lift/test_source_phase_lift.py",
        ],
    ),
    (
        "Replay the post-hoc repair-current diagnostic independently",
        [
            sys.executable,
            "code/thermodynamics/repair_current_orientation/verify_repair_current_orientation.py",
        ],
    ),
    (
        "Execute the repair-current provenance and mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/thermodynamics/repair_current_orientation/test_verify_repair_current_orientation.py",
        ],
    ),
    (
        "Validate the oriented-face three-norm certificate",
        [
            sys.executable,
            "code/b14_jacobi/oriented_face_bracket_selector.py",
            "--check",
        ],
    ),
    (
        "Replay the oriented-face three-norm certificate independently",
        [
            sys.executable,
            "code/b14_jacobi/verify_oriented_face_bracket_selector.py",
        ],
    ),
    (
        "Execute the oriented-face certificate mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/b14_jacobi/test_oriented_face_bracket_selector.py",
        ],
    ),
    (
        "Validate the retrospective alpha/HVP accounting verdict",
        [
            sys.executable,
            "code/particles/alpha_hvp_audit/build_alpha_hvp_verdict.py",
            "--verify",
        ],
    ),
    (
        "Execute the alpha/HVP provenance and mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/particles/alpha_hvp_audit/tests/test_alpha_hvp_verdict.py",
        ],
    ),
    (
        "Execute the source-derived fixed-capacity receipt gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/capacity_readback/test_source_derived_public_checkpoint_packet.py",
        ],
    ),
    (
        "Validate the bounded capacity-indexed counterfamily",
        [
            sys.executable,
            "code/capacity_readback/capacity_indexed_source_family.py",
            "--verify",
        ],
    ),
    (
        "Execute the bounded capacity-family mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/capacity_readback/test_capacity_indexed_source_family.py",
        ],
    ),
    (
        "Validate the complete-packet capacity lift receipts",
        [
            sys.executable,
            "code/capacity_readback/complete_packet_capacity_lift.py",
            "--check",
        ],
    ),
    (
        "Replay the complete-packet lift with the independent verifier",
        [
            sys.executable,
            "code/capacity_readback/verify_complete_packet_lift_independent.py",
        ],
    ),
    (
        "Execute the complete-packet lift mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/capacity_readback/test_complete_packet_capacity_lift.py",
        ],
    ),
    (
        "Validate the bounded direct N source verdict",
        [
            sys.executable,
            "code/capacity_readback/direct_n_closure_verdict.py",
            "--verify",
        ],
    ),
    (
        "Validate the typed retrospective N-closure branch packet",
        [
            sys.executable,
            "code/capacity_readback/n_closure_branch_certificate.py",
            "--validate-only",
        ],
    ),
    (
        "Execute the N-closure branch arithmetic and promotion gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/capacity_readback/test_n_closure_branch_certificate.py",
        ],
    ),
    (
        "Validate the source-only closure preflight inventory",
        [
            sys.executable,
            "code/closure_preflight/source_only_closure_preflight.py",
            "--verify",
        ],
    ),
    (
        "Replay the source-only closure preflight independently",
        [
            sys.executable,
            "code/closure_preflight/verify_source_only_closure_preflight_independent.py",
        ],
    ),
    (
        "Execute the source-only closure preflight mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/closure_preflight/test_source_only_closure_preflight.py",
        ],
    ),
    (
        "Validate the finite-cut to global-capacity attachment verdict",
        [
            sys.executable,
            "code/capacity_readback/global_capacity_attachment.py",
            "--validate-only",
        ],
    ),
    (
        "Replay the global-capacity attachment verdict independently",
        [
            sys.executable,
            "code/capacity_readback/verify_global_capacity_attachment_independent.py",
        ],
    ),
    (
        "Execute the global-capacity attachment mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/capacity_readback/test_global_capacity_attachment.py",
        ],
    ),
    (
        "Validate the named-law N closure verdict",
        [
            sys.executable,
            "code/capacity_readback/named_law_n_closure_verdict.py",
            "--verify",
        ],
    ),
    (
        "Replay the named-law N closure verdict independently",
        [
            sys.executable,
            "code/capacity_readback/verify_named_law_n_closure_verdict_independent.py",
        ],
    ),
    (
        "Execute the named-law N closure verdict mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/capacity_readback/test_named_law_n_closure_verdict.py",
        ],
    ),
    (
        "Validate the named-law horizon bridge short-circuit verdict",
        [
            sys.executable,
            "code/capacity_readback/named_law_horizon_bridge_verdict.py",
            "--verify",
        ],
    ),
    (
        "Execute the named-law horizon bridge mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/capacity_readback/test_named_law_horizon_bridge_verdict.py",
        ],
    ),
    (
        "Replay the named-law horizon bridge verdict independently",
        [
            sys.executable,
            "code/capacity_readback/verify_named_law_horizon_bridge_verdict_independent.py",
        ],
    ),
    (
        "Validate the discriminator stratum verdicts",
        [
            sys.executable,
            "code/discriminator_strata/stratum_verdicts.py",
            "--verify",
        ],
    ),
    (
        "Execute the discriminator stratum verdict tests",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/discriminator_strata/test_stratum_verdicts.py",
        ],
    ),
    (
        "Validate the angular template receipt",
        [
            sys.executable,
            "code/angular_sprint/angular_interpolant_certificate.py",
            "--verify",
        ],
    ),
    (
        "Execute the angular template tests",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/angular_sprint/test_angular_interpolant_certificate.py",
        ],
    ),
    (
        "Validate the kinetic ray receipt",
        [
            sys.executable,
            "code/angular_sprint/kinetic_ray_certificate.py",
            "--verify",
        ],
    ),
    (
        "Execute the kinetic ray tests",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/angular_sprint/test_kinetic_ray_certificate.py",
        ],
    ),
    (
        "Validate the carrier specificity receipt",
        [
            sys.executable,
            "code/discriminator_strata/carrier_specificity.py",
            "--verify",
        ],
    ),
    (
        "Execute the carrier specificity tests",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/discriminator_strata/test_carrier_specificity.py",
        ],
    ),
    (
        "Validate the horizon-record attachment verdict",
        [
            sys.executable,
            "code/capacity_readback/horizon_record_attachment_verdict.py",
            "--verify",
        ],
    ),
    (
        "Execute the direct N verdict mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/capacity_readback/test_direct_n_closure_verdict.py",
        ],
    ),
    (
        "Validate the invariant-mining pre-generation source lock",
        [
            sys.executable,
            "code/invariant_mining/tools/verify_pregeneration_freeze_independent.py",
        ],
    ),
    (
        "Execute the invariant-mining pre-generation mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/invariant_mining/tests/test_pregeneration_lock.py",
        ],
    ),
    (
        "Validate the forecast-contract generated state",
        [
            sys.executable,
            "code/particles/forecast_contract/forecast_contract.py",
            "--verify",
        ],
    ),
    (
        "Execute the forecast-contract adversarial gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/particles/forecast_contract/tests/test_forecast_contract.py",
        ],
    ),
    (
        "Validate the gravity premise ladder and its generated surface",
        [sys.executable, "tools/build_gravity_ladder.py", "--check"],
    ),
    (
        "Execute the gravity premise ladder mutation tests",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/collar_alignment/test_msa_characterizations.py",
            "code/geometry/test_collar_recoverability_receipts.py",
            "code/geometry/test_einstein_closure_receipts.py",
            "code/maxent/test_maxent_closure_acceptance.py",
            "code/geometry/test_modular_clock_instrumentation.py",
            "code/geometry/test_quotient_cap_readout.py",
            "code/particles/clock/test_clock_source_energy_closure_audit.py",
        ],
    ),
    (
        "Validate external-data provenance pins and declared source gaps",
        [sys.executable, "tools/check_external_data_provenance.py"],
    ),
    (
        "Validate clean-clone receipt path portability",
        [sys.executable, "tools/check_receipt_portability.py"],
    ),
    (
        "Validate public quantitative claim surfaces",
        [sys.executable, "tools/check_public_surface_claims.py"],
    ),
    (
        "Validate the null-model scorecard",
        [sys.executable, "tools/check_null_models.py", "--check"],
    ),
    ("Verify receipt promotion", [sys.executable, "tools/verify_receipt_promotion.py"]),
    (
        "Validate source-bound canonical book PDF assets",
        [sys.executable, "tools/book_pdf_assets.py"],
    ),
    (
        "Regression-test canonical book PDF assets",
        [sys.executable, "-m", "pytest", "-q", "tools/test_book_pdf_assets.py"],
    ),
    ("Validate paper release manifest", [sys.executable, "tools/validate_paper_release_manifest.py"]),
    ("Regression-test the manifest validator", [sys.executable, "-m", "pytest", "-q", "tools/test_paper_release_manifest.py"]),
    (
        "Regression-test offline Phase-0 scientific gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_public_surface_claims.py",
            "tools/test_null_models.py",
            "tools/test_receipt_portability.py",
            "tools/test_github_release_channel.py",
            "tools/test_gates_actually_fail.py",
        ],
    ),
    (
        "Regression-test the deterministic publication gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "paper/tools/test_check_build_warnings.py",
            "tools/test_reproducible_build_env.py",
            "tools/test_open_problem_ledger.py",
        ],
    ),
    (
        "Execute the Phase-0 proof and non-identifiability receipts",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/consensus/test_issue_517_proof_obligations.py",
            "code/particles/hierarchy/test_antecedent_only_nonidentifiability.py",
            "code/particles/hierarchy/test_hierarchy_bundle.py",
        ],
    ),
    (
        "Validate the exact discrete-refinement theorem packet",
        [
            sys.executable,
            "code/refinement/discrete_refinement_certificate.py",
            "--check",
        ],
    ),
    (
        "Replay the discrete-refinement packet independently",
        [
            sys.executable,
            "code/refinement/verify_discrete_refinement_independent.py",
        ],
    ),
    (
        "Execute the discrete-refinement mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/refinement/test_discrete_refinement_certificate.py",
        ],
    ),
    (
        "Validate the conditional dimension-six baryon-operator census",
        [
            sys.executable,
            "code/a5_closure/baryon_dimension_six_census.py",
            "--verify",
        ],
    ),
    (
        "Replay the dimension-six baryon census independently",
        [
            sys.executable,
            "code/a5_closure/verify_baryon_dimension_six_census_independent.py",
        ],
    ),
    (
        "Execute the dimension-six baryon census mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_baryon_dimension_six_census.py",
        ],
    ),
    (
        "Replay the bounded source-current capability certificate",
        [
            sys.executable,
            "code/a5_closure/source_current_capability_certificate.py",
            "verify",
            "--projection",
            "code/a5_closure/manifests/source_current_capability_projection.json",
            "--receipt",
            "code/a5_closure/receipts/source_current_capability.receipt.json",
        ],
    ),
    (
        "Verify the bounded source-current capability packet independently",
        [
            sys.executable,
            "code/a5_closure/verify_source_current_capability_independent.py",
            "--projection",
            "code/a5_closure/manifests/source_current_capability_projection.json",
            "--receipt",
            "code/a5_closure/receipts/source_current_capability.receipt.json",
        ],
    ),
    (
        "Execute the bounded source-current capability mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_source_current_capability_certificate.py",
        ],
    ),
    (
        "Execute the target-free current bracket search-space gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/issue_566_bracket_space_stage1/test_stage1.py",
            "code/a5_closure/issue_566_bracket_space_stage2/test_stage2.py",
        ],
    ),
    (
        "Execute the A5 fingerprint and frozen primitive-port prediction suites",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_fingerprint/test_a5_multipole_fixed_point_certificate.py",
            "code/a5_fingerprint/test_a5_multipole_fixed_point_hardening_certificate.py",
            "code/a5_fingerprint/test_a5_multipole_persistence_certificate.py",
            "code/a5_fingerprint/test_a5_multipole_persistence_independent_verifier.py",
            "code/a5_fingerprint/test_a5_multipole_persistence_stage4_certificate.py",
            "code/a5_fingerprint/test_a5_multipole_persistence_stage4_independent_verifier.py",
            "code/a5_fingerprint/test_spin_six_universality_certificate.py",
            "code/a5_fingerprint/test_spin_six_primitive_port_prediction.py",
            "code/a5_fingerprint/test_seam_current_edge_prediction.py",
            "code/a5_fingerprint/test_fz12_observation_map_certificate.py",
            "code/a5_fingerprint/test_fz12_full_symbol_remainder_certificate.py",
            "code/a5_fingerprint/test_fz12_synthetic_recovery_coverage.py",
            "code/a5_fingerprint/test_fz12_free_photon_hamiltonian.py",
            "code/a5_fingerprint/test_fz12_auger_threshold_diagnostic.py",
            "code/a5_fingerprint/test_carrier_scale_bound_diagnostic.py",
            "code/a5_fingerprint/test_carrier_frequency_speed_certificate.py",
        ],
    ),
    ("Check the claims scoreboard is regenerated", [sys.executable, "tools/build_scoreboard.py", "--check"]),
    ("Collect the mandatory scientific suite", [sys.executable, "-m", "pytest", "--collect-only", "-q", "code"]),
    ("Execute the audit fixture suite", [sys.executable, "-m", "pytest", "-q", "code/audit"]),
    ("Audit A5 closure ledgers", [sys.executable, "code/a5_closure/test_audit.py"]),
    (
        "Execute the three-axiom campaign certificate suites",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_equal_state_weights_certificate.py",
            "code/a5_closure/tests/test_response_grammar_completeness_certificate.py",
            "code/a5_closure/tests/test_a3_scheduler_kernel_certificate.py",
            "code/a5_closure/tests/test_routed_seam_grammar_certificate.py",
            "code/a5_closure/tests/test_optimizer_pushforward_certificate.py",
            "code/a5_closure/tests/test_multiplicity_window_certificate.py",
            "code/capacity_readback/test_capacity_semantics_menu_certificate.py",
            "code/a5_closure/tests/test_load_fiber_readback_certificate.py",
            "code/a5_closure/tests/test_noncentral_seam_reduction_certificate.py",
            "code/a5_closure/tests/test_record_counting_mechanism_certificate.py",
            "code/a5_closure/tests/test_scalar_sector_blindness_certificate.py",
            "code/a5_closure/tests/test_seam_grammar_matter_classification_certificate.py",
            "code/a5_closure/tests/test_family_band_attachment_certificate.py",
            "code/particles/flavor/test_entropy_w5_shape_certificate.py",
            "code/particles/leptons/test_koide_balance_comparison_certificate.py",
            "code/particles/scripts/test_build_postdiction_ledger.py",
            "code/a5_closure/tests/test_flux_defect_criterion_certificate.py",
            "code/cosmology/test_edge_center_clock_certificate.py",
            "code/consensus/test_compiled_lattice_settling_certificate.py",
            "code/thermodynamics/test_conditional_repair_certificate.py",
        ],
    ),
    (
        "Validate the conditional-repair thermodynamics certificate",
        [
            sys.executable,
            "code/thermodynamics/conditional_repair_certificate.py",
        ],
    ),
]

CERTIFICATE_STEPS: list[tuple[str, list[str]]] = [
    ("Execute the conditional port-current certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_port_current_inner_certificate.py"]),
    ("Execute the conditional matter-lift certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_super_tannakian_matter_lift_certificate.py"]),
    ("Execute the axis-center-descent certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_axis_center_descent_certificate.py"]),
    ("Execute the matter-menu spectral-ledger certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_matter_menu_spectral_ledger_certificate.py"]),
]

CERTIFICATE_SMOKE_STEPS: list[tuple[str, list[str]]] = [
    (
        "Recompute and verify the canonical port-current certificate",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_port_current_inner_certificate.py"
            "::PortCurrentInnerCertificateTests::test_reference_receipt_is_exactly_recomputable",
        ],
    ),
    (
        "Recompute and verify the canonical matter-lift certificate",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_super_tannakian_matter_lift_certificate.py"
            "::SuperTannakianMatterLiftTests::test_reference_receipt_is_exactly_recomputable",
        ],
    ),
    (
        "Recompute and verify the canonical axis-center-descent certificate",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_axis_center_descent_certificate.py"
            "::AxisCenterDescentTests::test_reference_receipt_is_exactly_recomputable",
        ],
    ),
    (
        "Recompute and verify the canonical matter-menu spectral ledger",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_matter_menu_spectral_ledger_certificate.py"
            "::MatterMenuSpectralLedgerTests::test_manifest_is_deterministic_and_matches_stored",
        ],
    ),
]


def run_steps(steps: list[tuple[str, list[str]]]) -> None:
    for title, command in steps:
        print(f"==> {title}", flush=True)
        result = subprocess.run(command, cwd=ROOT)
        if result.returncode != 0:
            raise SystemExit(f"FAILED: {title} (exit {result.returncode})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificates",
        action="store_true",
        help="also execute the exact certificate suites (adds ~26 minutes)",
    )
    parser.add_argument(
        "--certificates-only",
        action="store_true",
        help="execute only the exact certificate suites",
    )
    parser.add_argument(
        "--certificate-smoke-only",
        action="store_true",
        help="recompute and verify only the two canonical exact certificate receipts",
    )
    args = parser.parse_args()

    if args.certificates_only and args.certificate_smoke_only:
        parser.error("--certificates-only and --certificate-smoke-only are mutually exclusive")

    steps = [] if args.certificates_only or args.certificate_smoke_only else list(MANDATORY_STEPS)
    if args.certificates or args.certificates_only:
        steps += CERTIFICATE_STEPS
    if args.certificate_smoke_only:
        steps += CERTIFICATE_SMOKE_STEPS
    run_steps(steps)
    if args.certificates_only:
        print("certificate suites OK")
    elif args.certificate_smoke_only:
        print("certificate smoke suite OK")
    else:
        print("mandatory suite OK")


if __name__ == "__main__":
    main()
