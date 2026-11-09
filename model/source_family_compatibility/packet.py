"""Exact same-source clock, grading and gauge-attachment controls."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from . import algebra,bounds,rotations,hole_records,program_selection,defect_budget,local_matter,charge_bands,hole_dispersion,composite,charge_scaling,boost_current,cube_response,cube_spectrum
from . import momentum_response, charge_escape, uv_tail, local_spectral
from . import kinetic_nodes, nodal_memory, ground_marginals, cube_memory
from . import nodal_channels, spectral_residue, spectral_window_access
from . import nodal_locality, nodal_moment_bounds, mass_scale_selection

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def serialize(value):
    if type(value)is F:return str(value)
    if isinstance(value,dict):return {k:serialize(v)for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serialize(v)for v in value]
    return value


def bindings():
    paths=set(HERE.glob('*.py'));queue=list(paths)
    roots=(ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim',
           ROOT/'code/capacity_readback',ROOT)
    def include(path):
        if path.is_file()and path not in paths:
            paths.add(path)
            if path.suffix=='.py':queue.append(path)
    # These source modules are loaded by the capacity adapter using its explicit path.
    for name in program_selection.record_capacity.PUBLIC_SOURCES:include(ROOT/name)
    for name in kinetic_nodes.free_parent.bindings():include(ROOT/name)
    for name in ground_marginals.NATIVE_INPUTS:include(ROOT/name)
    while queue:
        path=queue.pop();base_root=next(root for root in roots if path.is_relative_to(root))
        package='.'.join(path.relative_to(base_root).with_suffix('').parts[:-1])
        for node in ast.walk(ast.parse(path.read_text())):
            names=[]
            if isinstance(node,ast.Import):names=[a.name for a in node.names]
            elif isinstance(node,ast.ImportFrom):
                name='.'*node.level+(node.module or '')
                base=importlib.util.resolve_name(name,package)if node.level else name
                names=[base]+[base+'.'+a.name for a in node.names]
            for name in names:
                parts=name.split('.')
                for source_root in roots:
                    include(source_root.joinpath(*parts).with_suffix('.py'))
                    for j in range(1,len(parts)+1):include(source_root.joinpath(*parts[:j])/'__init__.py')
    paths.update(ROOT/p for p in (
        'model/strong_coupling_refinement/README.md','model/protected_spin_limit/README.md',
        'model/bipartite_spin_limit/README.md','model/bipartite_spin_limit/LOCAL_LIMIT.md',
        'model/particle_mass_readout_gate/README.md','model/spin_canonical_field/README.md',
        'model/spin_canonical_energy/README.md','model/ballistic_free_boundary/README.md',
        'model/standard_model_source_map/PROOF.md','model/gauge_defect_transport/PROOF.md',
        'model/color_effective_action/PROOF.md','model/particle_gauge_band/QUANTUM_LINK_PROOF.md',
        'model/common_source/chiral_candidate/README.md','model/interacting_extension/chiral_electric/PROOF.md',
        'Lean/Screen/GaugeSectorBracketCompletion.lean','Lean/Screen/AssembledActionComposition.lean',
        'Lean/Screen/QuantumMatterIntegration.lean','Lean/Screen/FermionSectorAssembly.lean',
        'Lean/Screen/MatterGrammarIndexBridge.lean','Lean/Screen/ExteriorComponentBridge.lean',
        'Lean/Screen/SL2F5PortCover.lean','Lean/QFT/SectorInheritance.lean',
        'paper/tex_fragments/DERIVATION_TECHNICAL_SUPPLEMENT_PORT.tex',
        'model/native_program/README.md','model/observer/CAPACITY.md',
        'model/full_fock_field/README.md','model/centered_source_ground/README.md',
        'model/thermodynamic_ground/README.md',
        'model/retarded_spatial_locality/README.md','model/field_symmetry/README.md',
        'Lean/ObserverPatchHolography/RepairWordSchedule.lean',
        'Lean/ObserverPatchHolography/CapacityClosurePrinciple.lean',
        'model/spin_joint_field_limit/PROOF.md'))
    # The literal representation additionally reads its source geometry fixtures.
    from model.common_source import native_matter_attachment as matter
    paths.update(ROOT/name for name in matter.source_bindings())
    paths.update((HERE/'README.md',HERE/'PROOF.md',HERE/'REPAIR_TRAJECTORY.md',
                  HERE/'DEFECT_BUDGET.md',HERE/'LOCAL_MATTER_READOUT.md',
                  HERE/'CHARGED_SPECTRAL_BANDS.md',HERE/'ALIGNED_HOLE_DISPERSION.md',
                  HERE/'COMPOSITE_CAR.md',HERE/'CUBE_RESPONSE.md',
                  HERE/'STRONG_CHARGE_SCALING.md',HERE/'ENERGY_DIPOLE.md',HERE/'CUBE_SPECTRUM.md',
                  HERE/'MOMENTUM_RESPONSE.md',HERE/'LOCAL_CHARGE_ESCAPE.md',
                  HERE/'UV_TAIL.md',HERE/'LOCAL_SPECTRAL_READOUT.md',HERE/'KINETIC_NODE_MAP.md',
                  ROOT/'model/weyl_geometric_modular/README.md',
                  ROOT/'model/native_mass_channels/README.md',ROOT/'model/source_mass_shell/README.md',
                  HERE/'NODAL_MEMORY.md',HERE/'GROUND_MARGINALS.md',HERE/'CUBE_MEMORY.md',
                  HERE/'NODAL_CHANNELS.md',HERE/'SPECTRAL_RESIDUE.md',HERE/'SPECTRAL_WINDOW_ACCESS.md',
                  HERE/'NODAL_LOCALITY.md',HERE/'NODAL_MOMENT_BOUNDS.md',HERE/'MASS_SCALE_SELECTION.md'))
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted(paths)}


def claims():
    return {'same_source_full_ground_gap_reused':True,
        'all_U_at_least_four_cE_clock_constraint':True,
        'complete_charged_weight_and_original_odd_probability_bound':True,
        'chemical_shift_conjugate_sum_invariant':True,
        'finite_source_odd_and_spin_excitations_both_retained':True,
        'actual_internal_currents_and_Gauss_countercontrol_replayed':True,
        'selected_spin_source_proper_cubic_scalar_pullback':True,
        'local_even_source_creators_disjointly_commute':True,
        'actual_full_domain_one_hole_spin_instrument':True,
        'complete_C_E_coherence_and_declared_Poisson_generator':True,
        'conditional_counting_schedule_exact_native_attempt_decay':True,
        'global_attempt_complete_record_disturbance_bound':True,
        'actual_aligned_Hubbard_hole_sector_and_generic_leakage_checked':True,
        'actual_program_equal_capacity_different_coherent_response':True,
        'complete_repair_error_from_full_interacting_defect_budget':True,
        'fixed_coupling_complete_local_interacting_CAR_histories':True,
        'original_Weyl_action_with_bare_repair_locally_instantiated':True,
        'bare_repair_equated_across_hidden_spin_frames':False,
        'centered_ground_both_original_charge_weights_one_half':True,
        'volume_independent_positive_original_charged_band_weight':True,
        'source_generated_dimer_composite_clock_checked':True,
        'full_ground_source_composite_CAR_and_spectral_residual_constructed':True,
        'default_cube_composite_residual_strictly_positive_on_all_modes':True,
        'strong_contact_original_charge_weight_bound_uniform_in_volume':True,
        'full_centered_source_energy_dipole_identity':True,
        'complete_cube_wavepacket_energy_record_enclosures':True,
        'cubic_ground_composite_momentum_symbol_from_three_local_correlations':True,
        'complete_local_clock_recovers_original_energy_measures':True,
        'original_CAR_UV_tails_uniform_in_volume_and_packet':True,
        'default_ballistic_entire_unfiltered_unit_mass_regular_limit_excluded':True,
        'default_filtered_finite_physical_energy_weight_proved':False,
        'kinetic_node_winding_attached_to_proved_free_Weyl_family':True,
        'default_interacting_memory_positive_near_all_kinetic_nodes':True,
        'default_symmetric_ground_one_site_faithfulness_uniform':True,
        'complete_native_prepared_cube_local_occupation_records_executed':True,
        'full_default_cube_frequency_dependent_memory_certified':True,
        'cube_profile_labels_identified_with_momentum_eigenvalues':False,
        'energy_dipole_establishes_interacting_geometric_boost':False,
        'canonical_composite_family_appends_independent_source_flavours':False,
        'retarded_response_error_is_arbitrary_input_physical_channel_bound':False,
        'aligned_hole_source_position_translation_and_nonrelativistic_limit':True,
        'local_readout_equality_includes_remote_raw_archive':False,
        'full_infinite_sharp_energy_PVM_implemented_by_local_history_bound':False,
        'fixed_coupling_centered_state_identified_with_ordered_spin_limit':False,
        'classical_record_dephasing_is_global_dilation_information_loss':False,
        'capacity_and_service_counts_select_program':False,
        'repair_only_program_selects_spin_ground_and_clock':False,
        'one_hole_sector_is_interacting_Hubbard_ground_band':False,
        'excluded_all_OPH_fermionic_sectors':False,
        'local_dressed_theorem_assumed_full_ground_gap':False,
        'neutral_spin_fields_identified_with_original_odd_CAR':False,
        'negative_reference_excitation_relabelled_as_ground_gap':False,
        'internal_even_occupation_identified_as_physical_statistics':False,
        'same_source_nonabelian_chiral_spin_continuum_derived':False,
        'capacity_law_or_couplings_selected_by_gauge_covariance':False,
        'physical_particle_mass_or_Lorentz_spin_assigned':False,
        'cubic_scalar_profile_alone_proves_continuum_spin_zero':False,
        'source_DHR_category_or_physical_chirality_instantiated':False,
        'conditioning_used_to_restore_vanishing_original_weight':False}


def build():
    before=bindings();clocks=[]
    for edges in (1,2,12,96):
        ratio=F(4*edges)
        for p in (F(1,4),F(1,16),F(1,256)):
            # This is a theorem-bound specialization, not a selected trajectory.
            # c=R,U=R² gives J_eff=1; q=1/p preserves its ballistic clock.
            clocks.append(bounds.common_clock(ratio,ratio*ratio,edges,p,1/p))
    value={'schema':'oph.source_family_compatibility.v1','all_pass':True,
        'source_bindings':before,'same_clock_controls':serialize(clocks),
        'chemical_reference_controls':serialize([bounds.chemical_pair(F(3),F(4),F(n),F(2))for n in (-5,-1,0,2,5)]),
        'actual_full_CAR_dimer':algebra.dimer(),
        'actual_internal_current_and_complete_gauge_control':algebra.internal_and_transport(),
        'actual_source_spatial_rotation_and_exchange_control':rotations.certificate(),
        'actual_complete_one_hole_repair':hole_records.certificate(),
        'actual_capacity_and_native_program_control':program_selection.certificate(),
        'full_interacting_defect_budget':defect_budget.certificate(),
        'local_interacting_matter_history':local_matter.certificate(),
        'centered_original_charged_bands':charge_bands.certificate(),
        'source_generated_composite_CAR_response':serialize(composite.certificate()),
        'default_cube_composite_response':serialize(cube_response.certificate()),
        'volume_uniform_charged_scaling':serialize(charge_scaling.certificate()),
        'full_source_energy_dipole':serialize(boost_current.certificate()),
        'complete_cube_wavepacket_spectrum':serialize(cube_spectrum.certificate()),
        'actual_cubic_momentum_response':serialize(momentum_response.certificate()),
        'default_original_charge_escape':serialize(charge_escape.certificate()),
        'uniform_original_energy_tails':serialize(uv_tail.certificate()),
        'complete_local_source_energy_readout':serialize(local_spectral.certificate()),
        'kinetic_node_to_free_field_map':serialize(kinetic_nodes.certificate()),
        'actual_nodal_response_memory':serialize(nodal_memory.certificate()),
        'ground_marginal_preparation_and_native_records':ground_marginals.certificate(),
        'full_cube_frequency_dependent_memory':serialize(cube_memory.certificate()),
        'actual_nodal_scalar_channels':serialize(nodal_channels.certificate()),
        'actual_spectral_residue_and_subcontact_band':serialize(spectral_residue.certificate()),
        'complete_Lorentzian_window_access':serialize(spectral_window_access.certificate()),
        'actual_nodal_locality_and_window_survival':serialize(nodal_locality.certificate()),
        'uniform_full_CAR_star_force_bound':serialize(nodal_moment_bounds.certificate()),
        'fixed_source_scaling_and_mass_shell_scope':serialize(mass_scale_selection.certificate()),
        'aligned_hole_position_dispersion_and_continuum':serialize(hole_dispersion.certificate()),
        'claims':claims(),
        'replay_scope':{'actual_finite_graph_bound_producer_called':True,
            'actual_full_sixteen_state_source_reconstructed':True,
            'literal_M15_currents_and_color_transport_recomputed':True,
            'new_Hamiltonian_or_empirical_fit_introduced':False,
            'all_ancestral_native_words_replayed':False}}
    if before!=bindings():raise ValueError('source changed during source-family checks')
    return value


def validate(value):
    expected=build()
    if json.dumps(value,sort_keys=True,allow_nan=False)!=json.dumps(expected,sort_keys=True,allow_nan=False):
        raise ValueError('source, exact compatibility bound or physical scope changed')
    return expected


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',nargs='?',type=Path,const=HERE/'results.json')
    g.add_argument('--check',nargs='?',type=Path,const=HERE/'results.json')
    a=p.parse_args()
    if a.write:a.write.write_text(json.dumps(build(),indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:validate(strict_json.load(a.check))
    print('PASS: common-clock charge bounds and scoped gauge/matter source attachments')


if __name__=='__main__':main()
