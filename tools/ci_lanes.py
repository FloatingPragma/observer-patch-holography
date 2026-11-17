#!/usr/bin/env python3
"""Lane table for the push CI (ci.yml) and the nightly replays (nightly.yml).

Each lane is one independent replay: shell commands run from the repository
root with the lane's environment. A lane belongs to one tier ("push" runs on
every push and pull request, "nightly" runs on the schedule) and one install
profile (the requirements file its job installs). The runner packs the lanes
of one tier and profile into a fixed number of buckets by their measured
minutes, so a workflow matrix job only names its bucket:

    python tools/ci_lanes.py run --tier push --profile root --bucket 0 --buckets 8
    python tools/ci_lanes.py run --lane m1_fixed_noise
    python tools/ci_lanes.py list --tier push --profile root --buckets 8

Every lane in a bucket runs even after an earlier one fails; the bucket exits
non-zero when any lane failed and prints a summary with durations.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

BLAS1 = {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
CODE = {"PYTHONPATH": "code"}


@dataclass(frozen=True)
class Lane:
    name: str
    tier: str
    profile: str
    minutes: float
    commands: tuple[str, ...]
    env: dict[str, str] = field(default_factory=dict)


def lane(name, tier, profile, minutes, *commands, env=None):
    return Lane(name, tier, profile, minutes, tuple(commands), dict(env or {}))


def pytest_lane(name, minutes, path, *, tier="push", profile="root", env=None, flags="-q"):
    return lane(name, tier, profile, minutes, f"python -m pytest {flags} {path}", env=env)


COHERENT_CONSTRAINED = (
    "two_cell", "periodic_family", "field_limit", "periodic_spectrum", "higgs_composites",
    "spectral_edges", "oriented_links", "neutral_mobility", "family_cell", "neutral_local_spectrum",
    "neutral_composite", "periodic_boson_records", "peter_weyl_links", "family_transport",
    "local_reference", "matter_induced_ring", "matter_induced_plaquette", "finite_clock",
    "atomic_reference_control", "atomic_even_remainder", "charged_ground_exchange",
    "atomic_spin_geometry", "atomic_momentum_response", "ground_composite_response",
    "ground_scalar_response", "collective_spin_response", "interacting_ground_preparation",
    "interacting_charged_response", "interacting_charged_transport",
    "interacting_charge_localization", "interacting_charge_conjugation",
)
COHERENT_CONSTRAINED_MINUTES = {
    "periodic_spectrum": 6, "family_transport": 5, "neutral_composite": 4,
    "periodic_boson_records": 4, "higgs_composites": 3, "peter_weyl_links": 3,
}
COHERENT_COLLECTIVE = (
    "interacting_block_propagation", "interacting_block_quadrature", "interacting_block_scattering",
    "affine_gauss_source", "affine_gauss_mobile", "affine_gauss_spectrum",
    "affine_gauss_full_spectrum", "affine_gauss_spectrum_validation", "affine_gauss_charged_clock",
    "affine_gauss_clock_observables", "native_abelian_field_join", "native_abelian_energy_records",
    "native_block_dirac", "native_charge_resolution", "native_higgs_spectral_join",
    "native_charge_centroid", "native_weak_current_interface", "native_higgs_filtered_state",
    "native_weak_current_response", "native_higgs_current_response", "native_wilson_maxwell",
    "native_wilson_maxwell_refinement", "native_quark_response", "native_wilson_maxwell_ground",
    "native_wilson_maxwell_ground_refinement", "native_local_current_readout",
    "native_quark_higgs_join", "interacting_block_proca",
)
COHERENT_COLLECTIVE_MINUTES = {
    "native_quark_higgs_join": 8, "native_quark_response": 6, "native_local_current_readout": 3,
    "native_wilson_maxwell": 2,
}
COHERENT_DEFECT_WAVE = (
    "energy", "dispersion", "physical_attachment", "mobile_flux", "native_routing_join",
    "many_body_route", "mass_instrument", "propagating_native_join", "field_identity_audit",
    "chiral_source", "chiral_regulator", "chiral_higgs", "chiral_particle_map", "chiral_joint",
)

LANES: tuple[Lane, ...] = (
    # ------------------------------------------------------------------
    # Push tier, root requirements: code/ replays and their controls.
    pytest_lane("capacity_readback", 1, "code/capacity_readback -W error", env={"OPENBLAS_NUM_THREADS": "1"}),
    pytest_lane("collar_gap", 1, "code/yang_mills/test_collar_gap_certificate.py "
                "code/yang_mills/tests/test_finite_collar_gap_certificate.py "
                "code/yang_mills/tests/test_collar_certificate_accounting.py tools/test_strict_json.py",
                flags="-q -W error"),
    lane("consensus", "push", "root", 1.5, "python -W error -m pytest -q code/consensus", env=BLAS1),
    # Red on main since r2066: the ledger requires anchor text in a cosmology paper.
    lane("cosmology_accounting", "nightly", "root", 1,
         "python -m pytest -q code/cosmology/fixed_capacity_wlaw code/cosmology/postdiction_ledger -W error",
         "python code/cosmology/postdiction_ledger/verify_cosmology_postdiction_ledger_independent.py"),
    lane("dark_source_lensing", "push", "root", 2,
         "python -m dark_source_lensing.verify", "python -m pytest -q code/dark_source_lensing",
         env={**CODE, "OPENBLAS_NUM_THREADS": "1"}),
    lane("dispersion_feasibility", "push", "root", 0.5,
         "python -m dispersion_feasibility.check", "python -m pytest -q code/dispersion_feasibility", env=CODE),
    lane("federation_recovery", "push", "root", 2,
         "python -m pytest -q evidence/federation_recovery", "python evidence/federation_recovery/verify.py"),
    lane("finite_gauge_transfer", "push", "root", 2,
         "python -W error -m pytest -q code/yang_mills/tests/test_z2_finite_transfer_receipt.py "
         "code/yang_mills/tests/test_z2_finite_transfer_verifier.py code/yang_mills/tests/test_z2_transfer_spectrum.py "
         "code/yang_mills/tests/test_z2_kogut_precision.py code/yang_mills/tests/test_z2_full_space_controls.py "
         "code/yang_mills/tests/test_kogut_susskind_fiber_rate_instances.py",
         "python -W error -m pytest -q code/edge_sectors tools/test_post_r2029_audit_surfaces.py", env=BLAS1),
    pytest_lane("geometry", 2, "code/geometry", env={"OPENBLAS_NUM_THREADS": "1"}),
    lane("independent_postdictions", "push", "root", 1.5,
         "python -m pytest -q code/independent_postdictions",
         "python -m pytest -q code/particles/calibration/test_d11_criticality_comparison.py "
         "code/particles/flavor/test_audit_conditional_quark_ratio.py "
         "code/particles/flavor/test_compare_conditional_quark_masses.py "
         "code/particles/flavor/test_quark_lane_closure_contract.py "
         "code/particles/flavor/test_quark_p_driven_shared_evaluator_contract.py "
         "code/particles/leptons/test_koide_circulant_identity.py",
         env={**CODE, "OPENBLAS_NUM_THREADS": "1"}),
    lane("k7_gauge_kinetic_cross", "push", "root", 1,
         "python code/angular_sprint/k7_gauge_kinetic_cross_control.py --check",
         "python code/angular_sprint/verify_k7_gauge_kinetic_cross_control_independent.py "
         "--receipt code/angular_sprint/runtime/k7_gauge_kinetic_cross_control.json",
         "python -m pytest -q code/angular_sprint/test_k7_gauge_kinetic_cross_control.py",
         "python code/angular_sprint/kinetic_form_selection_certificate.py --verify",
         "python code/angular_sprint/verify_kinetic_form_selection_independent.py "
         "--receipt code/angular_sprint/runtime/kinetic_form_selection_receipt.json"),
    pytest_lane("m1_expander_archive", 1.5, "code/m1_expander_archive"),
    pytest_lane("m1_fermionic_source", 2, "code/m1_fermionic_source"),
    pytest_lane("m1_fixed_noise", 2.5, "code/m1_fixed_noise"),
    pytest_lane("m1_interfaces", 1, "code/m1_interfaces"),
    pytest_lane("m1_leakage", 1.5, "code/m1_leakage"),
    pytest_lane("m1_necessity", 1.5, "code/m1_necessity"),
    pytest_lane("m1_noisy_records", 1, "code/m1_noisy_records"),
    pytest_lane("m1_operational_clocks", 1, "code/m1_operational_clocks"),
    pytest_lane("m1_quantum_transport", 2.5, "code/m1_quantum_transport"),
    pytest_lane("m1_sea_energy", 1.5, "code/m1_sea_energy", env={**CODE, **BLAS1}),
    pytest_lane("m1_source_realization", 1.5, "code/m1_source_realization"),
    pytest_lane("m1_vacuum_fidelity", 1, "code/m1_vacuum_fidelity", flags="-q -W error"),
    lane("mckay_golden_field", "push", "root", 1.5,
         'python -B code/a5_closure/mckay_golden_field_certificate.py --output "$LANE_TMP/mckay-golden-raw.json"',
         'python -B code/a5_closure/verify_mckay_golden_field_independent.py --input "$LANE_TMP/mckay-golden-raw.json" '
         '--output "$LANE_TMP/mckay-golden-verified.json"',
         "python -B code/a5_closure/test_mckay_golden_field_certificate.py",
         "python -B code/a5_closure/test_mckay_golden_field_cross_control.py",
         "python -B code/a5_closure/verify_mckay_golden_field_cross_control.py"),
    lane("neutral_packet", "push", "root", 2,
         "python -m pytest -q -W error code/electromagnetism/test_neutral_packet_projection.py "
         "code/electromagnetism/test_neutral_packet_observables.py code/electromagnetism/test_whitney_quantum_packet.py",
         "python -m pytest -q -W error code/electromagnetism/test_whitney_interacting_quantum.py "
         "code/electromagnetism/test_whitney_kinetic_precision.py code/electromagnetism/test_whitney_kinetic_inputs.py "
         "code/electromagnetism/test_whitney_kinetic_replay.py code/electromagnetism/test_whitney_kinetic_scale.py "
         "code/electromagnetism/test_whitney_public_scale.py",
         env={"OPENBLAS_NUM_THREADS": "1"}),
    pytest_lane("paired_gravity", 2, "code/paired_gravity", env={**CODE, "OPENBLAS_NUM_THREADS": "1"}),
    lane("particle_mixing", "push", "root", 1,
         "python -m pytest -q code/particles/test_mixing_readout.py -W error",
         "python -m pytest -q code/particles/neutrino/test_neutrino_precision_candidate_lock.py "
         "code/particles/neutrino/test_neutrino_weighted_cycle_repair.py "
         "code/particles/neutrino/test_neutrino_weighted_cycle_shared_basis_representation.py "
         "code/particles/neutrino/test_physical_majorana_phase_theorem.py "
         "code/particles/flavor/test_quark_d12_mass_branch_and_ckm_residual.py "
         "code/particles/flavor/test_quark_local_basis_orbit_diagnostic.py "
         "code/particles/flavor/test_quark_transport_frame_diagnostic_orbit.py "
         "code/particles/flavor/test_quark_sigma_ud_orbit.py"),
    pytest_lane("quantum_information", 2.5,
                "code/quantum_information code/collar_alignment code/maxent code/native_response",
                env={**CODE, "OPENBLAS_NUM_THREADS": "1"}),
    lane("radial_inverse", "push", "root", 1,
         "python -W error -m pytest -q code/cosmology/test_oph_radial_lift_330.py "
         "code/cosmology/test_radial_inverse_regressions.py code/cosmology/test_radial_inverse_controls.py", env=BLAS1),
    pytest_lane("rg_principle", 1, "code/rg_principle"),
    lane("ringdown_arithmetic", "push", "root", 1, "python -W error -m pytest -q code/gravitation/ringdown_comb"),
    lane("source_checkpoint_selection", "push", "root", 3,
         "python -m source_checkpoint_selection.verify", "python -m pytest -q code/source_checkpoint_selection", env=CODE),
    pytest_lane("source_encoded_memory", 1, "code/source_encoded_memory"),
    pytest_lane("source_gravity_admissibility", 1, "code/source_gravity_admissibility",
                env={**CODE, "OPENBLAS_NUM_THREADS": "1"}),
    pytest_lane("source_native_accumulator", 1.5, "code/source_native_accumulator"),
    pytest_lane("source_native_programs", 1, "code/source_native_programs"),
    lane("source_native_program_families", "push", "root", 7,
         "python -m source_native_programs.family_build --workdir temp/native-program-witnesses "
         "--output temp/native-program-families.json",
         "python -m source_native_programs.family_verify --workdir temp/native-program-witnesses "
         "--families temp/native-program-families.json", env=CODE),
    pytest_lane("source_native_updates", 1, "code/source_native_updates"),
    pytest_lane("source_passive_memory", 1, "code/source_passive_memory"),
    lane("source_publication_selection", "push", "root", 1.5,
         "python -m source_publication_selection.verify", "python -m pytest -q code/source_publication_selection", env=CODE),
    pytest_lane("source_read_acceptance", 3.5, "code/source_read_acceptance"),
    lane("source_read_routing", "push", "root", 2,
         "python -m pytest -q code/source_read_routing code/source_feedback_transport/test_transport.py code/source_routing/test_routing.py",
         "python code/source_read_routing/support_comparison.py",
         'g++ -std=c++17 -O3 -Wall -Wextra code/source_read_routing/produce.cpp -o "$LANE_TMP/produce-routing"',
         'g++ -std=c++17 -O3 -Wall -Wextra code/source_read_routing/verify.cpp -o "$LANE_TMP/verify-routing"',
         'python code/source_read_routing/check_reproduction.py --producer "$LANE_TMP/produce-routing" --verifier "$LANE_TMP/verify-routing"',
         'ROUTING_VERIFIER="$LANE_TMP/verify-routing" python -m pytest -q code/source_read_routing/test_regeneration.py'),
    *(
        lane(f"source_read_routing_replay_q{q}_{variant}", "push", "root", 2,
             'g++ -std=c++17 -O3 -Wall -Wextra code/source_read_routing/produce.cpp -o "$LANE_TMP/produce-routing"',
             'g++ -std=c++17 -O3 -Wall -Wextra code/source_read_routing/verify.cpp -o "$LANE_TMP/verify-routing"',
             *(('python code/source_read_routing/native_negative_controls.py --binary "$LANE_TMP/verify-routing"',)
               if (q, variant) == (13, "baseline") else ()),
             f"python code/source_read_routing/regenerate.py evidence/source_net_causal_poset/routed_read_law/q{q}_{variant}.json "
             '--producer "$LANE_TMP/produce-routing" --verifier "$LANE_TMP/verify-routing"')
        for q in (13, 21) for variant in ("baseline", "source")
    ),
    pytest_lane("source_read_selection", 1, "code/source_read_selection"),
    pytest_lane("source_reusable_bus", 1, "code/source_reusable_bus"),
    pytest_lane("source_routing_refinement", 1, "code/source_routing_refinement"),
    lane("source_selection_model", "push", "root", 6.5,
         "python -m source_selection_model.verify", "python -m source_record_gluing.verify",
         "python -m pytest -q code/source_selection_model code/source_record_gluing", env=CODE),
    lane("source_temporal_acceptance", "push", "root", 1.5,
         "python -m source_temporal_acceptance.verify", "python -m pytest -q code/source_temporal_acceptance", env=CODE),
    lane("arithmon_interface_contracts", "push", "root", 5,
         "python -B code/a5_closure/verify_marked_triplet_attachment_independent.py",
         "python -B code/a5_closure/port_support_source_audit.py --check",
         "python -B code/a5_closure/verify_port_support_source_audit.py",
         "python -B code/a5_closure/equivariant_degree_crosscheck.py --check",
         "python -B code/a5_closure/verify_equivariant_degree_crosscheck.py",
         "python -B code/a5_closure/proper_code_attachment_probe.py --check",
         "python -B code/a5_closure/verify_proper_code_attachment_independent.py",
         "python -B code/a5_closure/c12_current_attachment_probe.py --check",
         "python -B code/a5_closure/verify_c12_current_attachment_independent.py",
         "python -B code/a5_closure/coherent_native_access_probe.py --check",
         "python -B code/a5_closure/verify_coherent_native_access_independent.py",
         "python -B code/a5_closure/coherent_decoder_probe.py --check",
         "python -B code/a5_closure/verify_coherent_decoder_independent.py",
         "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m pytest -q -W error model/native_interface/test_tomography.py "
         "model/native_interface/test_pointer.py model/compiled_native_interface/test_construction.py",
         "python -B code/a5_closure/source_current_order_sensitive_inventory.py verify",
         "python -B code/a5_closure/verify_source_current_order_sensitive_inventory.py "
         "--inventory code/a5_closure/manifests/source_current_order_sensitive_inventory.json",
         "python -B code/angular_sprint/k7_gauge_scale_contract.py --check",
         "python -B code/angular_sprint/verify_k7_gauge_scale_contract_independent.py",
         "python -B -m pytest -q code/a5_closure/tests/test_marked_triplet_attachment_certificate.py "
         "code/a5_closure/tests/test_port_support_source_audit.py code/a5_closure/tests/test_equivariant_degree_crosscheck.py "
         "code/a5_closure/tests/test_proper_code_attachment_probe.py code/a5_closure/tests/test_c12_current_attachment_probe.py "
         "code/a5_closure/tests/test_coherent_native_access_probe.py code/a5_closure/tests/test_coherent_decoder_probe.py "
         "code/angular_sprint/test_k7_gauge_scale_contract.py"),
    lane("certificate_smoke", "push", "root", 3,
         "python tools/run_mandatory_suite.py --certificate-smoke-only",
         "python -m pytest -q code/a5_closure/tests/test_source_repair_generator_certificate.py",
         "python -m pytest -q code/a5_closure/tests/test_native_defect_stability_gate.py"),

    # ------------------------------------------------------------------
    # Push tier, root requirements: model/ replays from the completion program.
    lane("model_sectors", "push", "root", 7,
         "python -m pytest -q --durations=10 -W error model/exchange_sector model/particle_screening", env=BLAS1),
    lane("model_repair", "push", "root", 5,
         "python -m pytest -q --durations=10 -W error model/local_repair_family model/active_code_energy "
         "model/screened_recovery model/conserved_sector_action model/screened_transport",
         "python tools/renew_diagnostic_parent_custody.py --check", env=BLAS1),
    lane("model_fields", "push", "root", 5,
         "python -m pytest -q --durations=10 -W error model/moving_defect model/record_action model/particle_obstruction model/interacting_field",
         "python -m model.moving_defect.experiment --check",
         "python -m model.record_action.action --check",
         "python -m model.particle_obstruction.classification --check",
         "python -m model.interacting_field --check model/interacting_field/results.json",
         "python -m pytest -q --durations=10 -W error model/bound_pair model/particle_clock model/field_regions "
         "model/field_stress model/field_tastes model/biased_reset",
         *(f"python -m model.{name} --check model/{name}/results.json"
           for name in ("bound_pair", "particle_clock", "field_regions", "field_stress", "field_tastes", "biased_reset")),
         env=BLAS1),
    lane("model_native", "nightly", "root", 7,
         "python -m pytest -q -W error model/band_clock model/native_localization model/native_program",
         "python -m pytest -q -W error model/pair_cyclic_control model/recovery_elision model/biased_code_universality",
         "python -m pytest -q -W error model/energy_identifiability model/regional_modular model/band_clock_budget model/compiled_native_history",
         "python -m pytest -q -W error model/regional_clock",
         "python -m model.regional_clock.verify --check",
         "python -m pytest -q -W error model/reference_observation model/native_fixed_algebra model/gauged_field model/full_fock_field",
         "python -m pytest -q -W error model/field_thermal model/support_refinement",
         "python -m model.support_refinement --check model/support_refinement/results.json",
         "python -m model.band_clock_budget --check model/band_clock_budget/results.json",
         "python -m model.band_clock --check model/band_clock/results.json",
         "python -m model.native_localization.experiment --check model/native_localization/results.json",
         env=BLAS1),
    lane("model_particle_services", "nightly", "root", 4,
         "python -m pytest -q --durations=10 -W error model/reference_recovery model/particle_response model/contact_scattering "
         "model/particle_comparison model/charged_dynamics model/particle_transport",
         "python -m pytest -q --durations=10 -W error model/modular_flow model/full_fock_energy model/flux_response "
         "model/field_protocols model/dimer_family",
         "python -m pytest -q --durations=10 -W error model/field_symmetry",
         "python -m model.field_symmetry --check model/field_symmetry/results.json",
         "python -m model.flux_response --check model/flux_response/results.json",
         "python -m model.particle_comparison.compare --check model/particle_comparison/results.json",
         "python -m model.charged_dynamics.packet --check model/charged_dynamics/results.json",
         "python -m model.particle_transport.packet --check model/particle_transport/results.json",
         env=BLAS1),
    lane("model_protected_defects", "push", "root", 2,
         "python -m pytest -q --durations=10 -W error model/defect_protection model/defect_identity model/defect_energy",
         "python -m model.defect_protection --check",
         "python -m model.defect_identity.identity --check",
         "python -m model.defect_energy --check", env=BLAS1),
    lane("model_native_field_readout", "push", "root", 2,
         "python -m pytest -q --durations=10 -W error model/defect_field_bridge model/defect_field_benchmark",
         *(f"python -m model.{name} --check" for name in (
             "defect_field_bridge.bridge", "defect_field_bridge.colored", "defect_field_bridge.dynamics_packet",
             "defect_field_benchmark", "defect_field_benchmark.readout", "defect_field_bridge.joint",
             "defect_field_bridge.coherent_extension", "defect_field_bridge.full_domain")),
         env=BLAS1),
    lane("model_native_field_dynamics", "push", "root", 1,
         "python -m model.defect_field_interactions.packet --check",
         "python -m model.defect_field_clock --check",
         "python -m model.defect_field_propagation --check",
         "python -m pytest -q --durations=10 -W error model/defect_field_propagation model/defect_field_clock model/defect_field_interactions",
         env=BLAS1),
    lane("model_native_field_thermodynamics", "nightly", "root", 5,
         "python -m model.defect_field_modular --check",
         "python -m model.defect_field_energy.packet --check",
         "python -m pytest -q --durations=10 -W error model/defect_field_energy model/defect_field_modular",
         "python -m pytest -q --durations=10 -W error model/spin_goldstone_sector model/native_spin_observer "
         "model/lattice_weyl_response model/cubic_particle_channels model/lattice_modular_limit model/spin_particle_scaling "
         "model/native_ordered_ground model/spin_shell_diagnostic tools/test_replay_model_receipts.py "
         "tools/test_read_ordered_spin_energies.py tools/test_particle_mass_table.py",
         *(f"python -W error -m model.{name} --check" for name in (
             "spin_goldstone_sector", "native_spin_observer", "lattice_weyl_response", "cubic_particle_channels",
             "lattice_modular_limit", "spin_particle_scaling", "native_ordered_ground", "spin_shell_diagnostic")),
         "python -W error tools/read_ordered_spin_energies.py --check model/completion/ordered_spin_energy_readout.json",
         "python -W error -m model.spin_canonical_energy --check",
         "python -W error -m model.shell_entropy_response --check",
         "python -W error tools/build_particle_mass_table.py --check",
         env=BLAS1),
    lane("model_finite_observables", "push", "root", 1,
         "python -m pytest -q -W error model/charged_composite_search model/quadratic_conservation "
         "model/prepared_regional_entropy",
         "python -W error -m model.charged_composite_search --check",
         "python -W error -m model.quadratic_conservation --check",
         "python -W error -m model.prepared_regional_entropy --check", env=BLAS1),
    lane("model_native_observations", "push", "root", 7,
         "python -m pytest -q -W error model/composite_spatial_response model/native_composite_probe "
         "model/interacting_observation_process model/native_regional_modular",
         "python -W error -m model.composite_spatial_response --check",
         "python -W error -m model.native_composite_probe --check",
         "python -W error -m model.interacting_observation_process --check",
         "python -W error -m model.native_regional_modular --check", env=BLAS1),

    # ------------------------------------------------------------------
    # Push tier, model requirements (model/requirements.txt).
    lane("coherent_seam_repair", "push", "model", 8,
         "python -m model.coherent_seam_repair --check",
         "python -m model.coherent_seam_repair.agreement --check",
         "python -m model.coherent_seam_repair.protection --check",
         "python -m pytest -q -W error model/coherent_seam_repair",
         *(cmd for pkg in ("particle_catalogue_correspondence", "finite_family_clock", "defect_particle_classes")
           for cmd in (f"python -m model.{pkg} --check", f"python -m pytest -q -W error model/{pkg}")),
         "python -m model.lepton_register_round --check",
         *(f"python -m model.coherent_defect_wave.{name} --check" for name in COHERENT_DEFECT_WAVE),
         "python -m pytest -q -W error model/coherent_defect_wave "
         + " ".join(f"--ignore=model/coherent_defect_wave/chiral_joint/{pkg}" for pkg in COHERENT_CONSTRAINED),
         env={**BLAS1, "PYTHONDONTWRITEBYTECODE": "1"}),
    lane("field_record_addenda", "nightly", "model", 7.5,
         "python -B -m model.field_record_addenda --check",
         "python -B -m pytest -q model/field_record_addenda -W error",
         "python -B tools/check_field_record_addenda_claims.py",
         "python -B -m pytest -q tools/test_field_record_addenda_claims.py -W error",
         "python -B tools/build_finite_field_records_atlas.py --check",
         env={**BLAS1, "PYTHONDONTWRITEBYTECODE": "1"}),
    lane("particle_extensions_hashes", "push", "model", 2,
         "python -m model.verify_particle_extensions --hashes-only",
         "python -m pytest -q -W error model/test_verify_particle_extensions.py"),

    # ------------------------------------------------------------------
    # Push tier, model requirements: per-package coherent-seam replays.
    *(
        lane(f"chiral_joint_{pkg}", "push", "model", COHERENT_CONSTRAINED_MINUTES.get(pkg, 1),
             f"python -m model.coherent_defect_wave.chiral_joint.{pkg} --check",
             f"python -m pytest -q -W error model/coherent_defect_wave/chiral_joint/{pkg}",
             env={**BLAS1, "PYTHONDONTWRITEBYTECODE": "1"})
        for pkg in COHERENT_CONSTRAINED
    ),
    *(
        lane(f"collective_{pkg}", "push", "model", COHERENT_COLLECTIVE_MINUTES.get(pkg, 1),
             f"python -m model.{pkg} --check", f"python -m pytest -q -W error model/{pkg}",
             env={**BLAS1, "PYTHONDONTWRITEBYTECODE": "1"})
        for pkg in COHERENT_COLLECTIVE
    ),

    # ------------------------------------------------------------------
    # Nightly tier, root requirements.
    lane("native_program_receipt", "nightly", "root", 3,
         "python -m model.native_program.receipt --check model/native_program/results.json", env=BLAS1),
    lane("support_wiring", "nightly", "root", 45,
         "python -m pytest -q code/support_wiring/test_support_wiring.py",
         "python code/support_wiring/verify.py",
         "python code/support_wiring/experiment.py trace --output temp/support-wiring-fresh",
         "python code/support_wiring/experiment.py controls --output temp/support-wiring-fresh",
         "python code/support_wiring/readouts.py provenance --output temp/support-wiring-fresh",
         "python code/support_wiring/readouts.py q13 --output temp/support-wiring-fresh",
         "python code/support_wiring/verify.py --part trace --archive temp/support-wiring-fresh",
         "python code/support_wiring/check_reproduction.py temp/support-wiring-fresh",
         "python evidence/exact_federation_L6_canonical_20260909/verify_archive.py",
         "python evidence/source_net_causal_poset/verify_causal_poset_archive.py", env=BLAS1),
)


def select(tier: str | None, profile: str | None) -> list[Lane]:
    return [x for x in LANES if (tier is None or x.tier == tier) and (profile is None or x.profile == profile)]


def pack(lanes: list[Lane], buckets: int) -> list[list[Lane]]:
    """Longest-processing-time packing; deterministic for a fixed lane table."""
    if buckets < 1:
        raise SystemExit("--buckets must be positive")
    bins: list[list[Lane]] = [[] for _ in range(buckets)]
    loads = [0.0] * buckets
    for item in sorted(lanes, key=lambda x: (-x.minutes, x.name)):
        index = min(range(buckets), key=lambda i: (loads[i], i))
        bins[index].append(item)
        loads[index] += item.minutes
    return bins


def run_lane(item: Lane, tmp_root: Path) -> tuple[bool, float]:
    tmp = tmp_root / item.name
    tmp.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update(item.env)
    env["LANE_TMP"] = str(tmp)
    start = time.monotonic()
    ok = True
    print(f"::group::lane {item.name}", flush=True)
    for command in item.commands:
        print(f"$ {command}", flush=True)
        result = subprocess.run(["bash", "-c", "set -euo pipefail\n" + command], cwd=ROOT, env=env)
        if result.returncode != 0:
            print(f"::error::lane {item.name} failed (exit {result.returncode}): {command}", flush=True)
            ok = False
            break
    print("::endgroup::", flush=True)
    return ok, time.monotonic() - start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=("run", "list"))
    parser.add_argument("--tier", choices=("push", "nightly"))
    parser.add_argument("--profile", choices=("root", "model"))
    parser.add_argument("--bucket", type=int)
    parser.add_argument("--buckets", type=int, default=1)
    parser.add_argument("--lane", action="append", default=[])
    args = parser.parse_args()

    names = [x.name for x in LANES]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise SystemExit("duplicate lane names: " + ", ".join(duplicates))

    if args.lane:
        by_name = {x.name: x for x in LANES}
        unknown = [n for n in args.lane if n not in by_name]
        if unknown:
            raise SystemExit("unknown lanes: " + ", ".join(unknown))
        groups = [[by_name[n] for n in args.lane]]
    else:
        lanes = select(args.tier, args.profile)
        if not lanes:
            raise SystemExit("no lanes match the selection")
        groups = pack(lanes, args.buckets)
        if args.bucket is not None:
            if not 0 <= args.bucket < args.buckets:
                raise SystemExit("--bucket must satisfy 0 <= bucket < --buckets")
            groups = [groups[args.bucket]]

    if args.action == "list":
        for index, group in enumerate(groups):
            load = sum(x.minutes for x in group)
            print(f"bucket {index}: {load:.1f} min: " + ", ".join(x.name for x in group))
        return

    tmp_root = Path(os.environ.get("RUNNER_TEMP") or ROOT / "temp") / "ci-lanes"
    results = []
    for group in groups:
        for item in group:
            ok, seconds = run_lane(item, tmp_root)
            results.append((item.name, ok, seconds))
    print("\nlane summary:")
    for name, ok, seconds in results:
        print(f"  {'PASS' if ok else 'FAIL'} {seconds:7.1f}s {name}")
    failed = [name for name, ok, _ in results if not ok]
    if failed:
        raise SystemExit("failed lanes: " + ", ".join(failed))


if __name__ == "__main__":
    main()
