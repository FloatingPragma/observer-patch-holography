"""Exact rational majorants for the actual-source reference theorem."""
from fractions import Fraction as F
from math import isqrt

from ..bipartite_spin_limit.bounds import encode
from ..bipartite_spin_limit.local import local_bounds
from ..strong_coupling_refinement.bounds import rational


def sqrt_upper(value):
    """A rational upper square root, exact on rational squares."""
    value = rational(value)
    product = value.numerator*value.denominator
    root = isqrt(product)
    if root*root != product:
        root += 1
    return F(root, value.denominator)


def reference_bounds(c, U, beta=1, support=1, volume=None):
    c, U = rational(c, positive=True), rational(U, positive=True)
    beta = rational(beta, positive=True)
    if type(support) is not int or support < 1:
        raise ValueError('positive integer support required')
    if volume is not None and (type(volume) is not int or volume < support):
        raise ValueError('integer volume at least the support required')
    local = local_bounds(c, U, support=support)
    B, nu = local['local_interaction_majorant'], local['onsite_penalty']
    eta, r = local['spin_generator_local_error'], local['off_diagonal_local_error']
    ground_m = 36*c/nu
    thermal_m = ground_m+sqrt_upper(2/(beta*nu))
    ground_d = 3*ground_m+8*B/nu
    thermal_d = 3*thermal_m+8*B/nu
    ground_p, thermal_p = ground_d**2, thermal_d**2
    if thermal_p > F(1, 4):
        raise ValueError('thermal charge majorant must be at most 1/4')
    entropy = 3*thermal_d
    lower_defect = eta+r+4*B*thermal_p+entropy/beta
    epsilon = eta+lower_defect
    return {
        'exchange': local['exchange'],
        'ground_bare_defect_density_upper': ground_m**2,
        'ground_dressed_onsite_defect_upper': ground_p,
        'ground_quadratic_form_error': support*(8*B*ground_p+2*eta+2*r),
        'thermal_positive_energy_density_upper': 2/beta,
        'thermal_bare_defect_density_upper': thermal_m**2,
        'thermal_dressed_onsite_defect_upper': thermal_p,
        'charge_entropy_per_site_upper': entropy,
        'spin_free_energy_defect_per_site': epsilon,
        'source_minus_spin_free_energy_density_lower': -lower_defect,
        'source_minus_spin_free_energy_density_upper': eta,
        'relative_entropy_per_site_upper': beta*epsilon,
        'finite_volume_trace_norm_error': None if volume is None else min(F(2), sqrt_upper(2*beta*volume*epsilon)),
        'local_ground_and_pressure_bounds_independent_of_volume': True,
        'thermal_global_trace_bound_independent_of_volume': False,
    }


def certificate():
    return encode({
        'scale_example': reference_bounds(2**128, 2**256, support=8, volume=64),
        'actual_ground_states_use_full_Fock_source': True,
        'ground_input_averaged_under_even_translations': True,
        'ground_limit_satisfies_spin_ground_condition': True,
        'full_source_Gibbs_reference_faithful': True,
        'modular_clock_source_time_factor': '-beta',
        'fixed_volume_thermal_spin_convergence': True,
        'uniform_volume_pressure_convergence': True,
        'unique_thermodynamic_phase_identified': False,
        'infinite_volume_KMS_identification_proved': False,
        'finite_c_Gibbs_invariant_under_native_repair': False,
        'source_state_selection_derived_from_repair': False,
        'interacting_relativistic_continuum': False,
        'integrated_physical_gate_accepted': False,
    })
