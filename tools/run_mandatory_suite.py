#!/usr/bin/env python3
"""Run the documented mandatory scientific suite from a clean checkout (#507).

This is the one command REPRODUCE.md documents and CI enforces:

    python tools/run_mandatory_suite.py                # standard suite (per-push CI)
    python tools/run_mandatory_suite.py --full         # + the heavy replay/mutation steps
    python tools/run_mandatory_suite.py --certificates # + exact certificate suites
    python tools/run_mandatory_suite.py --certificate-smoke-only

The standard suite runs every register validation, certificate check,
independent replay, and fast mutation gate. Five long-running replay and
mutation-scan steps (listed in HEAVY_STEP_TITLES, together most of the
suite's runtime) are excluded from the default run and enforced by --full,
which the nightly full job and the pre-release checklist execute; the
default run prints exactly which steps it skipped.

The mandatory suite validates the scientific claim registry, external inputs,
public quantitative surfaces, null-model controls, and the paper release
manifest, then proves those gates reject isolated false-green mutations. It
proves the full scientific collection imports cleanly (which is what keeps the
optional cloud/hardware lanes fail-closed) and executes the fast fixture
suites. The exact certificate suites (#566 port-current, #314 matter-lift) run
in their own CI workflow on their own triggers;
`--certificates` runs them here with the same commands.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MANDATORY_STEPS: list[tuple[str, list[str]]] = [
    (
        "Execute strict JSON duplicate-key gates",
        [sys.executable, "-m", "pytest", "-q", "tools/test_strict_json.py"],
    ),
    (
        "Execute protected-obstruction finite explorer",
        [
            sys.executable,
            "Lean/ObservableNormalForms/tools/verify_protected_obstruction_models.py",
        ],
    ),
    (
        "Regression-test the active-surface inventory guard",
        [sys.executable, "-m", "pytest", "-q", "tools/test_check_axiom_consistency.py"],
    ),
    (
        "Validate the active-surface inventory without rewriting it",
        [
            sys.executable,
            "tools/check_axiom_consistency.py",
            "--check-inventory",
        ],
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
        "Replay the B20 random-scan preflight offline algebra layer",
        [
            sys.executable,
            "code/b20_random_scan/validate_random_scan.py",
        ],
    ),
    (
        "Execute the B20 random-scan mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/b20_random_scan/test_b20_preflight.py",
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
        "Verify the recorded-decision phase-operation receipt independently",
        [
            sys.executable,
            "code/phase_operation_producer/verify_phase_operation_receipt.py",
        ],
    ),
    (
        "Execute the phase-operation producer determinism and mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/phase_operation_producer/test_phase_operation_receipt.py",
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
        "Replay the B14 compact-locus certificate independently",
        [sys.executable, "code/b14_jacobi/verify.py"],
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
        "Validate the V3 observation ledger surface",
        [
            sys.executable,
            "tools/build_observation_ledger.py",
            "--check",
        ],
    ),
    (
        "Execute the observation ledger gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_observation_ledger.py",
        ],
    ),
    (
        "Validate the V3 premise register surface",
        [
            sys.executable,
            "tools/build_premise_register.py",
            "--check",
        ],
    ),
    (
        "Execute the premise register gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_premise_register.py",
        ],
    ),
    (
        "Execute the Issue 750 adaptive-repair cross-surface gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_adaptive_repair_stratification_surfaces.py",
        ],
    ),
    (
        "Execute the Issue 750 cumulative attempt-capacity cross-surface gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_cumulative_attempt_capacity_surfaces.py",
        ],
    ),
    (
        "Execute the Issue 763 source-derived causal-order cross-surface gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_source_derived_causal_order_surfaces.py",
        ],
    ),
    (
        "Execute the Issue 750 fixed-federation cross-surface gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_fixed_federation_progress_surfaces.py",
        ],
    ),
    (
        "Execute the Issue 750 fixed-federation execution-classification gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_fixed_federation_execution_classification_surfaces.py",
        ],
    ),
    (
        "Kernel-check the typed fixed-federation execution audit",
        [sys.executable, "tools/check_fixed_federation_execution_audit.py"],
    ),
    (
        "Execute the fixed-unitary scattering obstruction controls",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/qft/test_finite_unitary_scattering_no_go.py",
        ],
    ),
    (
        "Execute the fixed-unitary scattering cross-surface gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_finite_unitary_scattering_surfaces.py",
        ],
    ),
    (
        "Execute the modal Maxwell factorization and mutation controls",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/electromagnetism/test_modal_maxwell_factorization.py",
        ],
    ),
    (
        "Execute the modal Maxwell cross-surface scope gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_modal_maxwell_factorization_surfaces.py",
        ],
    ),
    (
        "Execute the discrete Coulomb-Green replay, verifier, and mutation guards",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/electromagnetism/test_discrete_coulomb_green.py",
        ],
    ),
    (
        "Validate the V3 constants ancestry surface",
        [
            sys.executable,
            "tools/build_constants_ancestry.py",
            "--check",
        ],
    ),
    (
        "Execute the constants ancestry gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_constants_ancestry.py",
        ],
    ),
    (
        "Validate the emergent-instrument register surface",
        [
            sys.executable,
            "tools/build_instrument_register.py",
            "--check",
        ],
    ),
    (
        "Execute the instrument register gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_instrument_register.py",
        ],
    ),
    (
        "Validate the SM Lagrangian correspondence surface",
        [
            sys.executable,
            "tools/build_sm_correspondence.py",
            "--check",
        ],
    ),
    (
        "Execute the SM correspondence gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_sm_correspondence.py",
        ],
    ),
    (
        "Replay the invariant-metric phase certificate independently",
        [
            sys.executable,
            "code/b14_jacobi/verify_invariant_metric_phase.py",
        ],
    ),
    (
        "Execute the invariant-metric phase mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/b14_jacobi/test_invariant_metric_phase.py",
        ],
    ),
    (
        "Validate the port-dual metric-selector certificate",
        [
            sys.executable,
            "code/b14_jacobi/port_dual_metric_selector.py",
            "--check",
        ],
    ),
    (
        "Replay the port-dual metric-selector certificate independently",
        [
            sys.executable,
            "code/b14_jacobi/verify_port_dual_metric_selector.py",
        ],
    ),
    (
        "Execute the port-dual metric-selector mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/b14_jacobi/test_port_dual_metric_selector.py",
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
            "paper/tools/test_gen_warning_allowlist.py",
            "tools/test_reproducible_build_env.py",
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
    ("Collect the mandatory scientific suite", [sys.executable, "-m", "pytest", "--collect-only", "-q", "code"]),
    (
        # The E4 stale-absence guards need the pinned Mathlib source tree on
        # disk and fail closed without it, by design: an unevaluable citation
        # is a finding, never a green. This job never provisions Lean, so the
        # guards do not run here; they run in the Lean CI build job
        # (.github/workflows/lean-ci.yml), where the cited paths exist. The
        # --ignore flag keeps that exclusion visible in the executed command.
        "Execute the scientific validation fixtures (E4 stale-absence guards excluded here; they run in Lean CI)",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/audit",
            "--ignore=code/audit/test_e4_absence_guards.py",
        ],
    ),
    ("Validate A5 closure ledgers", [sys.executable, "code/a5_closure/test_audit.py"]),
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
    (
        "Execute the core-physics precursor cross-surface gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_core_physics_precursor_surfaces.py",
        ],
    ),
    (
        "Execute the mass-energy selection cross-surface gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_mass_energy_selection_surfaces.py",
        ],
    ),
    (
        "Execute the second-wave transport and collar cross-surface gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_second_wave_surfaces.py",
        ],
    ),
    (
        "Execute the third-wave force-law, speed-limit, clock-rate, and golden-split gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_third_wave_surfaces.py",
        ],
    ),
    (
        "Execute the fourth-wave clocked-chain, irreducibility, in-block, and flow gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_fourth_wave_surfaces.py",
        ],
    ),
    (
        "Execute the fourth-wave clocked-chain, irreducibility, in-block, and flow gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_fifth_wave_surfaces.py",
        ],
    ),
    (
        "Execute the post-r2029 no-go scope and nonlinear-obstruction gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/particles/calibration/test_d10_ew_tau2_current_carrier_obstruction.py",
            "code/particles/leptons/test_charged_entropic_branch_no_go.py",
            "code/particles/flavor/test_entropy_w5_shape_certificate.py",
            "tools/test_post_r2029_audit_surfaces.py",
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


# The heavy replay/mutation steps excluded from the default run and enforced
# by --full. Measured 2026-08-14 on the reference laptop these five steps take
# about 260s, 220s, 45s, 35s, and 22s; every other step finishes in a few
# seconds. Each entry must match a MANDATORY_STEPS title exactly; main()
# fails closed if one does not, so a renamed step cannot silently drop out of
# both modes.
HEAVY_STEP_TITLES: frozenset[str] = frozenset(
    {
        "Replay the B20 random-scan preflight offline algebra layer",
        "Execute the B20 random-scan mutation gates",
        "Execute the fixed-unitary scattering obstruction controls",
        "Execute the invariant-metric phase mutation gates",
        "Execute the complete-packet lift mutation gates",
    }
)


def run_steps(steps: list[tuple[str, list[str]]]) -> None:
    for title, command in steps:
        print(f"==> {title}", flush=True)
        result = subprocess.run(command, cwd=ROOT)
        if result.returncode != 0:
            raise SystemExit(f"FAILED: {title} (exit {result.returncode})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="also execute the heavy replay/mutation steps in HEAVY_STEP_TITLES",
    )
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

    mandatory_titles = {title for title, _ in MANDATORY_STEPS}
    unknown_heavy = HEAVY_STEP_TITLES - mandatory_titles
    if unknown_heavy:
        raise SystemExit(
            "HEAVY_STEP_TITLES entries missing from MANDATORY_STEPS: "
            + ", ".join(sorted(unknown_heavy))
        )

    if args.certificates_only or args.certificate_smoke_only:
        steps = []
    elif args.full:
        steps = list(MANDATORY_STEPS)
    else:
        steps = [
            step for step in MANDATORY_STEPS if step[0] not in HEAVY_STEP_TITLES
        ]
        print(
            "standard mode: skipping "
            f"{len(HEAVY_STEP_TITLES)} heavy steps (run with --full to include):"
        )
        for title in sorted(HEAVY_STEP_TITLES):
            print(f"  - {title}")
    if args.certificates or args.certificates_only:
        steps += CERTIFICATE_STEPS
    if args.certificate_smoke_only:
        steps += CERTIFICATE_SMOKE_STEPS
    run_steps(steps)
    if args.certificates_only:
        print("certificate suites OK")
    elif args.certificate_smoke_only:
        print("certificate smoke suite OK")
    elif args.full:
        print("mandatory suite OK (full)")
    else:
        print("mandatory suite OK (standard; heavy steps deferred to --full)")


if __name__ == "__main__":
    main()
