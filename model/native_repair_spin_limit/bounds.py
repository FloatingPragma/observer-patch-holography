"""Rational, volume-independent error bounds proved in PROOF.md."""
from fractions import Fraction as F
from math import isqrt

from ..bipartite_spin_limit.bounds import encode
from ..bipartite_spin_limit.local import local_bounds
from ..strong_coupling_refinement.bounds import rational


def sqrt_upper(value):
    value = rational(value)
    product = value.numerator*value.denominator
    root = isqrt(product)
    return F(root+(root*root < product), value.denominator)


def repair_bounds(c, U, rate=1, time=1, support=1, observed_edges=1,
                  positive_energy_density=None):
    """Unit-norm local readouts, complete records and a bounded-energy input.

    The source marginal must have the period-two symmetry in PROOF.md.
    A None energy upper bound means J=c*c/U, covering the explicit faithful
    input family. Zero covers bare protected inputs and averaged grounds.
    """
    c, U = rational(c, positive=True), rational(U, positive=True)
    rate, time = rational(rate, positive=True), rational(time)
    if type(observed_edges) is not int or observed_edges < 0:
        raise ValueError('nonnegative integer observed edge count required')
    base = local_bounds(c, U, time, support)
    B, nu, J = base['local_interaction_majorant'], base['onsite_penalty'], base['exchange']
    energy = J if positive_energy_density is None else rational(positive_energy_density)
    m = (36+1000*rate*time)*c/nu+sqrt_upper(energy/nu)
    p = (3*m+8*B/nu)**2
    eta, r = base['spin_generator_local_error'], base['off_diagonal_local_error']
    second_order = 48*B*B/nu
    repaired_map_error = 2**20*rate*B/nu
    delta = (2*B/3+6*second_order+96*rate)*p+eta+r+repaired_map_error
    cone = (1+3072*J*time)**3
    quantum = 2**18*support*time*cone*delta
    record = rate*observed_edges*time*(8*p+64*B/nu)
    output_dressing = 8*support*B/nu
    initial_dressing = 2**19*support*cone*B/nu
    faithful_weight = (c/nu)**2
    return {
        'exchange': J,
        'repair_rate_per_edge': rate,
        'positive_energy_density_upper': energy,
        'bare_defect_density_upper': m*m,
        'dressed_onsite_defect_upper': p,
        'diagonal_second_order_norm_upper': second_order,
        'repair_generator_conjugation_error': repaired_map_error,
        'comparison_interaction_error': delta,
        'spin_propagation_cone_factor': cone,
        'dressed_quantum_error_uncapped': quantum,
        'complete_record_error_uncapped': record,
        'dressed_joint_error': min(F(2), quantum+record),
        'bare_readout_joint_error': min(F(2), quantum+record+output_dressing),
        'bare_input_and_readout_joint_error': min(F(2), quantum+record+output_dressing+initial_dressing),
        'faithful_input_mixing_weight': faithful_weight,
        'faithful_input_and_bare_readout_joint_error': min(F(2), quantum+record+output_dressing+initial_dressing+2*faithful_weight),
        'penalty_over_exchange': nu/J,
        'depends_on_total_volume': False,
        'rate_tends_to_zero': False,
    }


def growing_time_exponents(alpha):
    """At fixed J, Gamma and energy, T=O(c**alpha), 0<=alpha<1/6."""
    alpha = rational(alpha)
    if alpha >= F(1, 6):
        raise ValueError('the proved growing-time window requires alpha < 1/6')
    return {
        'time_power': alpha,
        'dressed_defect_power': 2*alpha-2,
        'quantum_error_power': 6*alpha-1,
        'observed_record_error_power': alpha-1,
        'initial_dressing_error_power': 3*alpha-1,
        'expected_observed_events_power': alpha,
        'all_errors_vanish': True,
    }


def certificate():
    fixed = repair_bounds(2**256, 2**512, rate=1, time=1, support=8, observed_edges=12)
    # T=c**(1/12), so expected retained events diverge along this family.
    growing = repair_bounds(2**768, 2**1536, rate=1, time=2**64,
                            support=8, observed_edges=12)
    if max(fixed['faithful_input_and_bare_readout_joint_error'],
           growing['faithful_input_and_bare_readout_joint_error']) >= F(1, 10**6):
        raise ArithmeticError('explicit scale certificate failed')
    return encode({
        'fixed_time_scale_certificate': fixed,
        'growing_time_scale_certificate': growing,
        'growing_time_exponents': growing_time_exponents(F(1, 12)),
        'native_repair_rate_fixed_positive_on_every_edge': True,
        'period_two_source_marginal_required': True,
        'untouched_reference_allowed': True,
        'complete_local_coin_raw_time_records_retained': True,
        'coherent_coin_compared_to_matching_plus_state': True,
        'classical_coin_readout_is_explicit_dephasing': True,
        'arbitrary_record_feedback_covered': False,
        'infinite_archive_global_trace_norm_claim': False,
        'volume_independent_local_interacting_limit': True,
        'faithful_full_source_inputs_constructed': True,
        'finite_c_faithful_inputs_stationary': False,
        'finite_c_source_ground_gap_assumed': False,
        'interacting_relativistic_continuum': False,
        'simultaneous_finite_charged_and_spin_scales': False,
        'Einstein_dynamics': False,
        'integrated_physical_gate_accepted': False,
    })
