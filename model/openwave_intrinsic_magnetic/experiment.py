"""Convergent finite stationary solves and independent observable controls."""
from functools import lru_cache
import numpy as np
from .solver import native_reduction,solve
from .checks import derivative_controls,gauge_control,dipole_profile,external_probe


def summary(result):
    assert result['eigen_residual']<1e-10
    assert result['fixed_point_residual']<2e-11
    assert result['maxwell_residual']<1e-9
    assert result['current_divergence_norm']<1e-10
    assert result['total_current_norm']<1e-10
    assert result['field_gauge_divergence_norm']<1e-10
    return {**{key:result[key] for key in ('n','spacing','field_n','coupling','mass','trap','wilson','spin_sign','external_B')},
        'unreduced_internal_spatial_dimension':16*result['n']**3,
        'positive_D_eigenvalue':round(result['eigenvalue'],10),
        'spin_expectation':round(result['spin'],10),
        'stationary_current_moment':np.round(result['moment'],9).tolist(),
        'field_energy':round(result['field_energy'],12),
        'empty_vacuum_relative_excitation_plus_field_energy':round(result['vacuum_relative_rest_energy'],10),
        'g_using_supplied_mass_coefficient':round(float(result['g_supplied_mass']),8),
        'g_using_actual_bound_excitation_energy':round(float(result['g_actual_rest_energy']),8),
        'current_norm':round(result['current_norm'],10),
        'numerical_residual_controls':{'eigenstate_below_1e_minus_10':True,
            'fixed_point_below_2e_minus_11':True,'Maxwell_below_1e_minus_9':True,
            'current_continuity_and_net_current_below_1e_minus_10':True,
            'transverse_gauge_below_1e_minus_10':True,
            'interval_certificate_of_the_computed_solution':False},
        'relaxed_fixed_point_steps':result['steps'],
        'observed_tail_contraction_ratio':round(float(max(result['history'][i+1]/result['history'][i]
            for i in range(3,len(result['history'])-1))),5)}


def clean_floats(value):
    if isinstance(value,dict):return {k:clean_floats(v) for k,v in value.items()}
    if isinstance(value,list):return [clean_floats(v) for v in value]
    if isinstance(value,(float,np.floating)):return float(f'{float(value):.8g}')
    return value


@lru_cache(maxsize=1)
def run():
    reduction=native_reduction()
    family=[solve(n,a,N) for n,a,N in ((5,1.,19),(7,.75,25),(9,.6,31))]
    base=family[0]
    assert all(family[i+1]['g_actual_rest_energy']>family[i]['g_actual_rest_energy'] for i in (0,1))
    low=solve(coupling=.05);high=solve(coupling=.2)
    assert abs(low['moment'][2]-high['moment'][2])>1e-5
    reversed_state=solve(spin_sign=-1)
    assert np.linalg.norm(reversed_state['current']+base['current'])<2e-9
    assert np.linalg.norm(reversed_state['field']+base['field'])<2e-9
    perturbed=solve(initial_field=1.4*base['field'],relaxation=.7)
    assert abs(perturbed['vacuum_relative_rest_energy']-base['vacuum_relative_rest_energy'])<1e-10
    assert np.linalg.norm(perturbed['field']-base['field'])<2e-9
    no_wilson=solve(wilson=0)
    assert no_wilson['current_norm']<1e-10
    large_box=solve(field_n=31)
    assert abs(large_box['moment'][2]-base['moment'][2])<1e-5
    controls=derivative_controls();gauge=gauge_control(base);probe=external_probe(base)
    return {'schema':'oph.intrinsic_stationary_magnetic.v1',
        'native_complete_algebra':reduction,
        'declared_construction':{'matter_kernel':'K(A)=|D_W(A)| on all native internal and spatial rows; number-conserving full Fock extension',
            'Dirac_operator':'sum gamma_j i(T_j-T_j^dag)/(2a)+beta[m+eta |x|^2+r sum_j(2-T_j-T_j^dag)/(2a)]',
            'gauged_link':'T_(x,x+e_j)=exp(i a A_e)',
            'charge_convention':'q=-1 with conventional p=-i partial and alpha_physical=-gamma; J=-delta<K>/delta A',
            'empty_Fock_vacuum_energy':'0, at zero transverse field',
            'mass_coefficient':8,'mass_well_coefficient':.01,'Wilson_coefficient':.5,
            'self_field_energy':'a^3 ||C A||^2/(2 kappa), kappa>0',
            'field_boundary':'periodic containing cube, mean-zero transverse A; matter has open Dirichlet links',
            'source_law':'semiclassical stationary expectation current coupled to the finite transverse Maxwell field',
            'binding_and_electrostatics':'External scalar mass well and open box; longitudinal electrostatic energy and its feedback are excluded by the declared transverse-only architecture. No neutralization or screening mechanism is derived.',
            'physical_units_or_measured_g_used_in_coefficients':False,
            'original_PW119_or_number_nonconserving_Hamiltonian_identified_with_this_kernel':False},
        'exact_finite_results':{'all_sixteen_rows_and_eight_folded_valleys_retained':True,
            'both_signed_D_branches_retained_by_absolute_value':True,
            'positive_kernel_gap_lower':'m for every real link potential A',
            'empty_Fock_vacuum_unique_at_fixed_A':True,
            'finite_one_particle_stationary_minimizer_exists':'By compact matter sphere, positive-kernel continuity and coercive mean-zero transverse field energy; this does not identify or certify the displayed numerical branch.',
            'source_current_includes_Wilson_link_derivative':True,
            'constant_current_harmonic_and_divergence_vanish_at_exact_stationarity':True,
            'current_derived_without_inserted_Pauli_or_g_parameter':True},
        'fixed_matter_box_refinement':[summary(row) for row in family],
        'coupling_family':[summary(row) for row in (low,base,high)],
        'independent_matrix_derivative_controls':clean_floats(controls),
        'gauge_controls':clean_floats(gauge),
        'intrinsic_dipole_profile':clean_floats(dipole_profile(base)),
        'larger_field_box_profile':clean_floats(dipole_profile(large_box)),
        'weak_external_energy_probe':clean_floats(probe),
        'stationary_controls':{'opposite_spin_reverses_actual_current_and_self_field':True,
            'perturbed_start_and_changed_relaxation_recover_same_branch':True,
            'self_field_feedback_changes_current_when_coupling_changes':True,
            'naive_central_difference_mass_well_current_norm_below_1e_minus_10':True,
            'naive_zero_current_control_is_same_as_Wilson_source':False},
        'scope':{'stationary_intrinsic_current_and_its_self_consistent_finite_magnetic_field':True,
            'coupled_current_and_Maxwell_residuals_numerically_small':True,
            'specific_nonzero_dipole_branch_interval_certified':False,
            'unconfined_self_bound_electron':False,
            'full_longitudinal_transverse_quantum_EM_dynamics':False,
            'source_selected_mass_or_electromagnetic_coupling':False,
            'precision_g2_or_anomalous_moment_prediction':False,
            'coupled_time_dependent_Larmor_process':False,
            'primitive_repair_law_implementation_transfer':False},
        'verdict':'PASS_DECLARED_STATIONARY_INTRINSIC_MAGNETIC_SOURCE'}
