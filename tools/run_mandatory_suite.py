#!/usr/bin/env python3
"""Run the documented mandatory scientific suite from a clean checkout (#507).

This is the one command REPRODUCE.md documents and CI enforces:

    python tools/run_mandatory_suite.py                # standard suite (per-push CI)
    python tools/run_mandatory_suite.py --full         # + the heavy replay/mutation steps
    python tools/run_mandatory_suite.py --certificates # + exact certificate suites
    python tools/run_mandatory_suite.py --certificate-smoke-only
    python tools/run_mandatory_suite.py --shard-index 0 --shard-count 2

CI can partition the standard suite into nonempty contiguous shards. The
zero-based shard index and shard count must be supplied together; sharding
cannot be combined with full or certificate modes. Concatenating the shard
command lists reproduces the unchanged standard list in order.

The standard suite runs every register validation, certificate check,
independent replay and fast mutation gate; CI runs it in six shards beside
the replay lanes. HEAVY_STEP_TITLES is empty, so --full equals the standard
run. The exhaustive certificate suites run on demand with --certificates.

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
        "Replay the retrospective ACT tilt posterior audit",
        [sys.executable, "code/cosmology/tilt_posterior/audit.py", "--check"],
    ),
    (
        "Replay the retrospective ACT tilt likelihood slice",
        [sys.executable, "code/cosmology/tilt_posterior/likelihood_slice.py", "--check"],
    ),
    (
        "Verify proof-domain, sampling, spectrum and radial-inverse counterexamples",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/audit/test_recurrent_proof_boundaries.py",
            "code/P_derivation/test_thomson_spectral_transport.py",
            "code/P_derivation/test_source_spectral_theorem.py",
            "code/P_derivation/test_thomson_endpoint_interval_certificate.py",
        ],
    ),
    (
        "Verify sampled causal-manifold arithmetic and mutation controls",
        [sys.executable, "evidence/source_net_causal_poset/verify_sampled_manifold.py"],
    ),
    (
        "Replay the local Standard Model jet action and Grassmann Ward controls",
        [sys.executable, "-m", "pytest", "-q", "code/sm_local_action/test_local_action.py"],
    ),
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
    (
        "Validate reader-facing prose style across papers, essays and registers",
        [sys.executable, "tools/check_reader_style.py"],
    ),
    (
        "Replay the Maxwell-lane receipts against their pinned Lean sources",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/electromagnetism/test_discrete_coulomb_green.py",
            "code/electromagnetism/test_local_face_maxwell_action.py",
        ],
    ),
    (
        "Validate Lean comment style and the in-file kernel-axiom audit blocks",
        [sys.executable, "tools/check_lean_docstring_style.py"],
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
        "Validate the seam-Maxwell continuum receipt",
        [sys.executable, "code/electromagnetism/seam_maxwell_continuum.py", "--check"],
    ),
    (
        "Execute independent continuum, forcing and false-green controls",
        [sys.executable, "-m", "pytest", "-q",
         "code/electromagnetism/test_seam_maxwell_continuum.py"],
    ),
    (
        "Execute serial feedback, decoded Maxwell action and provenance controls",
        [sys.executable, "-m", "pytest", "-q",
         "code/electromagnetism/test_serial_maxwell_readout.py"],
    ),
    (
        "Execute refining causal-cone paths and authenticated local-law controls",
        [sys.executable, "-m", "pytest", "-q",
         "code/causal_refinement/test_refining_cone.py",
         "code/causal_refinement/test_source_net_causet.py",
         "code/causal_refinement/test_source_count_clock.py"],
    ),
    (
        "Execute cone Whitney interpolation, full action derivatives and custody controls",
        [sys.executable, "-m", "pytest", "-q",
         "code/electromagnetism/test_cone_whitney_bridge.py"],
    ),
    (
        "Execute stationary Whitney volume, exact stability, readback and charged-data controls",
        [sys.executable, "-m", "pytest", "-q",
         "code/electromagnetism/test_whitney_maxwell_dynamics.py",
         "code/electromagnetism/test_whitney_charged_matter.py"],
    ),
    (
        "Execute coupled Whitney spatial, interacting Hilbert and charged evolution controls",
        [sys.executable, "-m", "pytest", "-q",
         "code/electromagnetism/test_whitney_spatial_consistency.py",
         "code/electromagnetism/test_whitney_interacting_quantum.py",
         "code/electromagnetism/test_whitney_quantum_state.py",
         "code/electromagnetism/test_whitney_charged_dynamics.py"],
    ),
    (
        "Execute same-action real continuum, charged observer, dynamical clock and quantum trial controls",
        [sys.executable, "-m", "pytest", "-q",
         "code/electromagnetism/test_whitney_real_continuum.py",
         "code/electromagnetism/test_whitney_charged_instrument.py",
         "code/electromagnetism/test_whitney_ephemeris_clock.py",
         "code/electromagnetism/test_whitney_quantum_history.py",
         "code/electromagnetism/test_whitney_charged_enclosure.py",
         "code/electromagnetism/test_whitney_charged_checkpoint.py",
         "code/electromagnetism/test_whitney_magnetic_continuum.py",
         "code/electromagnetism/test_whitney_quantum_packet.py"],
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
        "Validate the angular template receipt",
        [
            sys.executable,
            "code/angular_sprint/angular_interpolant_certificate.py",
            "--verify",
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
            "code/particles/hierarchy/test_interval_witness.py",
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
        "Verify the source-current tomography Stage 0 baseline pins",
        [sys.executable, "code/a5_closure/source_current_tomography_stage0.py", "verify"],
    ),
    (
        "Verify the source-current tomography Stage 1 admissibility contract",
        [
            sys.executable,
            "code/a5_closure/source_current_tomography_stage1_contract.py",
            "verify",
        ],
    ),
    (
        "Execute the source-current tomography mutation gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_source_current_tomography_stage0.py",
            "code/a5_closure/tests/test_source_current_tomography_stage1_contract.py",
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
            "code/a5_fingerprint/test_a5_multipole_fixed_point_hardening_certificate.py",
            "code/a5_fingerprint/test_a5_multipole_persistence_certificate.py",
            "code/a5_fingerprint/test_a5_multipole_persistence_independent_verifier.py",
            "code/a5_fingerprint/test_a5_multipole_persistence_stage4_certificate.py",
            "code/a5_fingerprint/test_a5_multipole_persistence_stage4_independent_verifier.py",
            "code/a5_fingerprint/test_carrier_scale_bound_diagnostic.py",
            "code/a5_fingerprint/test_carrier_frequency_speed_certificate.py",
            "code/a5_fingerprint/test_carrier_class_dispersion_certificate.py",
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
            "code/a5_closure/tests/test_routed_seam_grammar_certificate.py",
            "code/a5_closure/tests/test_optimizer_pushforward_certificate.py",
            "code/capacity_readback/test_capacity_semantics_menu_certificate.py",
            "code/capacity_readback/test_toy_readback.py",
            "code/capacity_readback/test_cp2_p4_premise_reduction.py",
            "code/a5_closure/tests/test_load_fiber_readback_certificate.py",
            "code/a5_closure/tests/test_noncentral_seam_reduction_certificate.py",
            "code/a5_closure/tests/test_record_counting_mechanism_certificate.py",
            "code/a5_closure/tests/test_coset_carrier_certificate.py",
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
        "Independently replay the protected-record memory certificate",
        [
            sys.executable,
            "code/thermodynamics/protected_memory/verify_protected_memory.py",
        ],
    ),
    (
        "Execute the protected-record memory mutation controls",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/thermodynamics/protected_memory/test_protected_memory.py",
        ],
    ),
    (
        "Execute the post-r2029 no-go scope and nonlinear-obstruction gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_post_r2029_audit_surfaces.py",
        ],
    ),
    (
        "Verify the triple-carrier noncommutative operator join",
        [sys.executable, "-m", "pytest", "-q", "code/source_operator_join/test_join.py"],
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


# Steps excluded from the standard run and enforced by --full. None: every step runs in the CI shards.
HEAVY_STEP_TITLES: frozenset[str] = frozenset()


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
    parser.add_argument(
        "--shard-index", type=int,
        help="zero-based index of a standard-suite shard (requires --shard-count)",
    )
    parser.add_argument(
        "--shard-count", type=int,
        help="number of nonempty standard-suite shards (requires --shard-index)",
    )
    args = parser.parse_args()

    if args.certificates_only and args.certificate_smoke_only:
        parser.error("--certificates-only and --certificate-smoke-only are mutually exclusive")
    sharded = args.shard_index is not None or args.shard_count is not None
    if sharded:
        if args.shard_index is None or args.shard_count is None:
            parser.error("--shard-index and --shard-count must be supplied together")
        if args.full or args.certificates or args.certificates_only or args.certificate_smoke_only:
            parser.error("sharding is available only for the standard suite")

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
    if sharded:
        total = len(steps)
        if not 1 <= args.shard_count <= total:
            parser.error(f"--shard-count must be between 1 and {total}; shards must be nonempty")
        if not 0 <= args.shard_index < args.shard_count:
            parser.error("--shard-index must satisfy 0 <= index < --shard-count")
        start = total * args.shard_index // args.shard_count
        stop = total * (args.shard_index + 1) // args.shard_count
        steps = steps[start:stop]
        print(
            f"standard shard {args.shard_index + 1}/{args.shard_count}: "
            f"steps {start + 1}..{stop} of {total}", flush=True,
        )
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
    elif sharded:
        print(
            f"mandatory shard OK ({args.shard_index + 1}/{args.shard_count}; "
            "standard mode)"
        )
    else:
        print("mandatory suite OK (standard; heavy steps deferred to --full)")


if __name__ == "__main__":
    main()
