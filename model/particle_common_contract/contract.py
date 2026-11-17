"""Construct and authenticate the primary finite observer-field program.

The program is a declared architecture. Compatibility is checked before
admitting a consequence; a source hash alone does not prove that consequence.
"""
from __future__ import annotations
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from . import provenance
from model.particle_gauss_ground_state import gauss
from model.particle_gauge_characteristics import characteristics
from model.particle_spatial_ground_state import ground

ROOT=Path(__file__).resolve().parents[2]
PROGRAM='oph.open_wilson_gauss_observer.v1'
SOURCE_INPUTS=(
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_gauss_ground_state/gauss.py',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_spatial_ground_state/ground.py',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_gauge_characteristics/characteristics.py',
    'model/particle_spatial_modular_response/regional.py',
    'model/particle_spatial_modular_response/PROOF.md',
    'model/particle_physical_boundary_state/boundary.py',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_primary_locality/locality.py',
    'model/particle_primary_locality/PROOF.md',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_local_wilson_density/density.py',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/interacting_extension/_source.py',
    'model/interacting_extension/source_manifest.json',
    'model/common_source/overlap_regulator.py',
    'model/observer_particles/results/representation_identification.json',
    'model/observer_particles/results/representation_identification_operators.json',
)
SOURCE_INPUTS=provenance.definition_inputs(set(SOURCE_INPUTS)|set(gauss.SOURCE_INPUTS)
                          |set(ground.SOURCE_INPUTS)|set(characteristics.SOURCE_INPUTS))


def serial(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v) for v in value]
    return value


def digest(value):
    return hashlib.sha256(json.dumps(serial(value),sort_keys=True,
        separators=(',',':'),allow_nan=False).encode()).hexdigest()


def member(level=0):
    """A fixed finite sequence. Its joint physical limit is a separate claim."""
    if type(level) is not int or level<0:
        raise ValueError('a nonnegative integer level is required')
    side=2**(2*level+1);spacing=F(1,2**level);degree=level+1
    target=F(1,2**(level+4))
    plan=gauss.spatial_plan(side=side,spacing=spacing,link_degree=degree,
        electric='105',magnetic='1/7',target_excess=target)
    source={
        'program':PROGRAM,'level':level,'side':side,'spacing':spacing,
        'implementation_sha256':provenance.definition_hashes(SOURCE_INPUTS),
        'per_link_word_degree':degree,'CAR_modes_per_cell':60,
        'gauge_representation':'native M15, four Dirac copies, all occupations',
        'electric_coefficient':F(105),'magnetic_coefficient':F(1,7),
        'regulator':'finite_range_Wilson','boundary':'open_cube_zero_external_charge',
        'Gauss_constraint':'all owned-link endpoints, including ghost endpoints',
        'word_projection':'tensor_product_per_link_Haar_orthogonal',
        'word_basis':'length-then-lexicographic native entries/conjugates, identity first; exact positive-norm Haar Gram-Schmidt',
        'Fock_basis':'site lexicographic, Dirac row0..3, native internal row0..14',
        'link_basis':'owner site lexicographic, direction0..2',
        'energy_zero':'unchanged negative plaquette trace; no added scalar',
        'onsite_one_body':'3 beta_D/a',
        'interior_forward_hop':'(-beta_D+beta_D gamma_i) R(U_i)/(2a), plus Hermitian adjoint',
        'electric_action':'105 sum_link C_trace/(2a)',
        'magnetic_action':'-sum_interior_plaquettes ReTr_M15 R(U_p)/(105a)',
        'compression_order':'form complete parent products, then compress',
        'owned_links':3*side**3,'all_gauge_vertices':side**3+3*side**2,
        'interior_hopping_bonds':3*side**2*(side-1),
        'interior_plaquettes':3*side*(side-1)**2,
        'coframe':'identity; source coordinates are a declared dictionary',
    }
    cooling=plan['cooling']
    state={
        'source_id':digest(source),'target_ground_excess':target,
        'beta':cooling['beta'],'Taylor_order':cooling['Taylor_order'],
        'dyadic_scale':cooling['dyadic_scale'],
        'spectral_lower':plan['source']['spectral_lower'],
        'spectral_upper':plan['source']['spectral_upper'],
        'amplitude':'Q=P_G [T_m(beta(U I-H)/(2s))]^s',
        'density':'Q^2/Tr(Q^2)','purification':'vec(Q)/sqrt(Tr(Q^2))',
        'seed':'empty Fock and constant-link blank, with conjugate blank',
        'preparation':'canonical Householder; identity if target equals blank',
        'preparation_primitive':'exact algebraic Householder map declared as a program service',
        'source_only_Gauss_support':True,'postselection':False,
        'exact_Gibbs_state':False,'source_stationary':True,
        'exact_matrix_construction_executed_at_full_spatial_size':False,
        'compiled_approximation_stationarity_claimed':False,
    }
    reference={'source_id':state['source_id'],'state_id':digest(state),
        'carrier':'full conjugate source carrier','idle_H':0,
        'separately_Gauss_invariant':True,'retained_in_complete_outputs':True,
        'joint_purification_stationary_under_source_only_clock':False,
        'counter_evolution':'only an explicitly declared apparatus gate'}
    clock={'source_id':state['source_id'],'energy_unit':'source action units; hbar=1',
        'Schrodinger':'exp(-i t H)','Heisenberg':'exp(+i t H) A exp(-i t H)',
        'modular_generator_on_support':'log(Tr Q^2) I - 2 log(Q)',
        'modular_equals_minus_beta_times_source_time':False,
        'GeV_or_proper_time_conversion':None}
    ids={'source_id':state['source_id'],'state_id':digest(state),
         'reference_id':digest(reference),'clock_id':digest(clock)}
    result={'source':source,'state':state,'reference':reference,'clock':clock,
            'ids':ids,'preparation_plan':plan,
            'approximate_realization':{'half_trace_allowance':F(1,2**(level+20)),
                'is_the_exact_canonical_preparation':False,
                'deterministic_compiler_and_tape_must_be_bound_separately':True},
            'matched_characteristics':characteristics.clock_family('1/7','105')}
    result['member_id']=digest(result)
    return serial(result)


def verify_member(value):
    """Recompute the specified source and state; do not trust supplied ids."""
    try: expected=member(value['source']['level'])
    except (KeyError,TypeError) as error:raise ValueError('incomplete common member') from error
    if value!=expected:raise ValueError('source, state, reference or clock differs from the declared program')
    return expected


def expected_result(value, claim):
    """Only registered, recomputed consumers can supply admitted results."""
    verified=verify_member(value)
    if claim=='complete local boundary records on the selected physical source':
        from model.particle_primary_locality import locality
        return locality.common_contract_payload(verified)
    if claim=='physical boundary-charge preparation and complete central readout':
        from model.particle_physical_boundary_state import boundary
        return boundary.common_contract_payload(verified)
    if claim=='faithful regional support of the selected source-Gauss polynomial state':
        from model.particle_spatial_modular_response import regional
        return regional.common_contract_payload(verified)
    producers={
        'full-carrier source dimension and energy bounds':
            verified['preparation_plan']['source'],
        'specified source-only Gauss polynomial preparation':
            verified['preparation_plan'],
        'gauge and matter principal characteristic match':
            verified['matched_characteristics'],
    }
    if claim not in producers:
        raise ValueError('unregistered consequence; a verified producer is required')
    return producers[claim]


def bind_result(value, claim, result):
    if serial(result)!=expected_result(value,claim):
        raise ValueError('result differs from the declared consequence producer')
    return {'program':PROGRAM,'member_id':value['member_id'],**value['ids'],
            'claim':claim,'result':serial(result),'result_sha256':digest(result)}


def verify_result(value, receipt):
    verify_member(value)
    for key,expected in {'program':PROGRAM,'member_id':value['member_id'],**value['ids']}.items():
        if receipt.get(key)!=expected:raise ValueError('incompatible consequence: '+key)
    if receipt.get('result_sha256')!=digest(receipt.get('result')):
        raise ValueError('result content changed after binding')
    if receipt.get('result')!=expected_result(value,receipt.get('claim')):
        raise ValueError('result differs from the declared consequence producer')


# Complete certificate inventory. These controls keep their own program/state.
# None of their detector laws is implicitly admitted as a primary-field law.
CONTROL_DOMAINS={
    'neutral_scalar':'scalar_Qstar0', 'charged_scalar':'scalar_Qstar0',
    'complete_apparatus':'scalar_Qstar0_and_declared_compiler',
    'thermal_clock':'scalar_Qstar0_exact_Gibbs',
    'gauge_field':'gauge_parent_theorems_and_directional_control',
    'holonomy_transport':'finite_sign_link_controller',
    'consensus_records':'finite_record_comparison',
    'gauge_input':'multiple_declared_gauge_controls',
    'bare_effective':'scalar_Qstar0_coefficients',
    'controller_readout':'three_cycle_controller',
    'defect_field_attachment':'declared_memory_or_holonomy_scalar_couplings',
    'field_causality':'parameterized_parent_operator_bound',
    'energy_and_gravity_response':'declared_scalar_and_rotor_action_responses',
    'full_scalar_ground':'scalar_Qstar6_full_ground',
    'regional_modular_response':'scalar_Qstar0_internal_family_region',
    'gauge_field_preparation':'periodic_total_degree_word_seed',
    'full_ground_charged_records':'scalar_Qstar6_compensated_charged_ground',
    'prepared_spatial_field_records':'periodic_total_degree_short_time_control',
    'local_gauge_density_family':'mixed_ambient_spatial_density',
    'local_wilson_spectral_density':'mixed_ambient_spectral_density',
    'fixed_time_local_polynomial_records':'mixed_polynomial_field_family',
    'exact_sign_response_boundary':'native_sign_sensitivity_counterfamily',
    'smooth_interacting_local_records':'positive_smoothing_action',
    'joint_position_translation_records':'periodic_translation_and_reference_family',
    'refinement_uniform_causality':'parameterized_parent_operator_bound',
    'interacting_field_thermal_records':'native208_warm_beta1over64_full',
    'declared_spatial_geometry_response':'parameterized_coframe_and_native208_warm_control',
    'spatial_low_energy_preparation':'ambient_product_cutoff_and_native208_cold_control',
    'ground_energy_spectral_records':'native208_cold_beta4_full',
    'gauge_matter_characteristics':'quadratic_parent_family_multiple_coefficients',
    'Gauss_constrained_preparation':'parameterized_open_Gauss_family_and_native78_control',
    'accessible_charged_threshold':'native208_ground_trial_and_residue',
    'spatial_regional_modular_response':'parameterized_region_theorem_and_warm_internal_controls',
    'fixed_volume_link_cutoff_limit':'fixed_volume_parent_and_rotor_cutoff_comparisons',
    'physical_boundary_state':'selected_primary_member_physical_boundary_readout',
    'primary_locality':'selected_primary_complete_boundary_records_and_local_influence',
    'continuum_gap_transfer':'conditional_same_clock_decay_and_finite_transfer_controls',
}


def validate_inventory(report):
    actual={k for k,v in report.items() if isinstance(v,dict)}-{'input_sha256','common_source_contract'}
    if actual!=set(CONTROL_DOMAINS):
        raise ValueError('unclassified or omitted certificate domains: '+str(sorted(actual^set(CONTROL_DOMAINS))))
    if report.get('physical_particle_assignments') or report.get('source_predicted_GeV_masses'):
        raise ValueError('physical identification requires a separately proved dictionary')


def certificate():
    from model.particle_spatial_modular_response import regional
    from model.particle_physical_boundary_state import boundary
    from model.particle_primary_locality import locality
    selected=member(0); next_member=member(1)
    bounds=selected['preparation_plan']['source']
    checks={
        'source_bounds':bind_result(selected,'full-carrier source dimension and energy bounds',bounds),
        'physical_preparation':bind_result(selected,'specified source-only Gauss polynomial preparation',selected['preparation_plan']),
        'common_quadratic_clock':bind_result(selected,'gauge and matter principal characteristic match',selected['matched_characteristics']),
        'physical_regional_support':regional.common_contract_bound(selected),
        'physical_boundary_preparation':boundary.common_contract_bound(selected),
        'complete_local_records':locality.common_contract_bound(selected),
    }
    for receipt in checks.values():verify_result(selected,receipt)
    return {'schema':'oph.common_source_contract.v1','program':PROGRAM,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_INPUTS},
        'fixed_primary_member':selected,'next_declared_member':next_member,
        'bound_consequences':checks,
        'comparison_certificate_domains':{k:v for k,v in CONTROL_DOMAINS.items() if k not in ('physical_boundary_state','primary_locality')},
        'primary_certificate_domains':{k:CONTROL_DOMAINS[k] for k in ('physical_boundary_state','primary_locality')},
        'item1_fixed_construction_specification':True,
        'item2_finite_physical_preparation_with_boundary_charges':True,
        'item3_finite_complete_boundary_records_with_uniform_locality':True,
        'architecture_uniqueness_required':False,
        'narrow_native_repair_interpreter_equivalence_proved':False,
        'full_spatial_matrix_or_gate_tape_executed':False,
        'joint_physical_refinement_limit_proved':False,
        'observed_particle_or_Einstein_identification':False,
        'empirical_calibrations_used_in_primary_program':[]}
