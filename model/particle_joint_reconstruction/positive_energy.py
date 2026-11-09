"""Selected-source instance of the existing downward-weight transfer theorem.

The source density is stationary; its idle-reference purification need
not be. A continuous limit of ordered correlations is a stated premise.
"""
from __future__ import annotations

from fractions import Fraction as F

import sympy as s

from model.particle_common_contract import contract
from model.particle_ground_spectral_records import records
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_ground_spectral_records/records.py',
    'model/particle_ground_spectral_records/PROOF.md',
    'Lean/QFT/GroundRecordTail.lean',
)


def selected_tail_bound(member, threshold='1/4', operator_norm_squared='1'):
    """An exact bound for any source-Gauss-preserving probe with this norm cap."""
    contract.verify_member(member)
    delta = ground._q(threshold, 'positive spectral threshold', True)
    cap = ground._q(operator_norm_squared, 'squared operator norm bound')
    epsilon = F(member['state']['target_ground_excess'])
    expected = F(1, 2**(member['source']['level']+4))
    if epsilon != expected or not member['state']['source_stationary']:
        raise ValueError('the selected stationary cooling target differs')
    bound = cap*records.negative_transfer_bound(epsilon, delta)
    return contract.serial({
        'member_id': member['member_id'], **member['ids'],
        'source_level': member['source']['level'],
        'mean_physical_ground_excess_upper': epsilon,
        'negative_transfer_threshold': delta,
        'operator_norm_squared_upper': cap,
        'original_negative_spectral_weight_upper': bound,
        'formula': 'min(C_squared, C_squared*epsilon/delta)',
        'correlation': 'omega(A_dagger alpha_t(A)), alpha_t(A)=exp(i t H) A exp(-i t H)',
        'frequency': 'final source energy minus initial source energy',
        'probe_domain': 'source algebra; A maps ran(P_G) into ran(P_G)',
        'source_density_stationary': True,
        'idle_reference_purification_assumed_stationary': False,
        'original_weights_normalized_by_probe_success': False,
        'bare_charged_probe_uses_this_floor': False,
        'ordered_correlation_is_sequential_pointer_probability': False,
        'fixed_time_correlation_convergence_proved': False,
        'limiting_nonzero_finite_energy_weight_proved': False,
    })


def bind_member(member, threshold='1/4', operator_norm_squared='1'):
    payload = selected_tail_bound(member, threshold, operator_norm_squared)
    return {'schema': 'oph.selected_positive_energy.v1',
            'result': payload, 'result_sha256': contract.digest(payload)}


def verify_bound(member, receipt, threshold='1/4', operator_norm_squared='1'):
    if receipt != bind_member(member, threshold, operator_norm_squared):
        raise ValueError('positive-energy receipt differs from the selected source or result')


def triangular_filter_bounds(member, half_width='1', operator_norm_squared='1', threshold='1/4'):
    """Finite time-filter moment bounds; no filtered apparatus execution claim."""
    T = ground._q(half_width, 'positive temporal half width', True)
    cap = ground._q(operator_norm_squared, 'squared operator norm bound')
    tail = selected_tail_bound(member, threshold, cap)
    return contract.serial({
        'selected_tail': tail, 'temporal_half_width': T,
        'filter': 'f_T(s)=(1-abs(s)/T)_+/T',
        'filter_L1_norm': F(1), 'filter_derivative_L1_norm': 2/T,
        'spectral_multiplier': 'sinc(omega*T/2)^4',
        'spectral_total_mass_upper': cap,
        'spectral_first_absolute_moment_upper': 2*cap/T,
        'correlation_uniform_Lipschitz_constant': 2*cap/T,
        'continuous_response_subsequence_for_fixed_T_and_uniform_norm_cap': True,
        'subsequence_effectively_selected': False,
        'nonzero_or_nonzero_energy_response_proved': False,
        'same_point_time_observation_domain': False,
        'joint_field_algebra_or_apparatus_execution_proved': False,
    })


def rotated_three_level_control(excess='1/16'):
    """Independent exact matrix check; this finite control is not the primary."""
    epsilon = ground._q(excess, 'control excess', True)
    if epsilon > 1:
        raise ValueError('the finite control uses 0<epsilon<=1')
    e = s.Rational(epsilon.numerator, epsilon.denominator)
    rotation = s.Matrix([[s.Rational(3, 5), -s.Rational(4, 5), 0],
                         [s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, 1]])
    energies = (0, 1, 3)
    h = rotation*s.diag(*energies)*rotation.T
    rho = rotation*s.diag(1-e/2, e/4, e/4)*rotation.T
    cycle = s.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    probe = rotation*cycle*rotation.T
    if rho*h != h*rho or probe.H*probe != s.eye(3) or s.trace(rho) != 1:
        raise ArithmeticError('the exact stationary control or norm premise failed')
    projectors = []
    for energy in energies:
        p = s.eye(3)
        for other in energies:
            if other != energy:
                p = p*(h-other*s.eye(3))/s.Integer(energy-other)
        projectors.append(p)
    weights = []
    spectrum = {}
    for i, pi in enumerate(projectors):
        for j, pj in enumerate(projectors):
            weight = s.trace(rho*pi*probe.H*pj*probe*pi)
            if weight < 0:
                raise ArithmeticError('ordered Gram spectral weight is negative')
            weight = F(weight)
            frequency = energies[j]-energies[i]
            spectrum[frequency] = spectrum.get(frequency, F(0))+weight
            weights.append({'initial': i, 'final': j, 'frequency': frequency, 'weight': weight})
    if sum(row['weight'] for row in weights) != 1 or F(s.trace(rho*h)) != epsilon:
        raise ArithmeticError('the complete spectral weights or energy excess differ')
    return {'H': s.ImmutableMatrix(h), 'rho': s.ImmutableMatrix(rho),
            'A': s.ImmutableMatrix(probe), 'projectors': tuple(projectors),
            'all_nine_transition_rows': weights, 'spectral_measure': spectrum,
            'mean_ground_excess': epsilon,
            'idle_purification_energy_variance': F(s.trace(rho*h*h)-s.trace(rho*h)**2),
            'is_selected_spatial_source': False}


def escaping_energy_control(n):
    """Finite ground-probe measures illustrating why continuity at0 matters.

    The normalized probe vector has n² equal amplitudes1/n at energies
    k/n. Its characteristic functions tend to zero at each fixed t!=0,
    while their value at0 is one. This is not a stationary-primary limit.
    """
    n = ground._i(n, 'escape control size', 1)
    return {'excited_energies': [F(k, n) for k in range(1, n*n+1)],
            'each_original_weight': F(1, n*n),
            'total_original_weight': F(1),
            'mean_excitation_energy': F(n*n+1, 2*n),
            'source_preparation_ground_excess': F(0),
            'continuous_limiting_characteristic_function': False,
            'is_selected_spatial_source': False}


def certificate():
    control = rotated_three_level_control()
    rows = control['all_nine_transition_rows']
    checks = []
    for delta in (F(1, 2), F(1), F(2), F(4)):
        actual = sum(row['weight'] for row in rows if row['frequency'] <= -delta)
        upper = records.negative_transfer_bound(control['mean_ground_excess'], delta)
        if actual > upper:
            raise ArithmeticError('the independent finite matrix control violates the reused tail bound')
        checks.append({'threshold': delta, 'actual_original_weight': actual, 'upper': upper})
    return contract.serial({
        'schema': 'oph.selected_positive_energy_certificate.v1',
        'selected_member_bounds': [bind_member(contract.member(n)) for n in (0, 2, 6)],
        'independent_three_level_control': {
            'all_nine_transition_rows': rows,
            'mean_ground_excess': control['mean_ground_excess'],
            'source_is_stationary': True,
            'idle_purification_energy_variance': control['idle_purification_energy_variance'],
            'is_selected_spatial_source': False, 'tail_checks': checks},
        'existing_finite_Lean_theorem': 'OPH.QFT.GroundRecordTail.downwardWeight_le_excess_div',
        'fixed_triangular_filter': triangular_filter_bounds(contract.member(0)),
        'continuous_limit_positive_support_is_conditional': True,
        'joint_refinement_or_noncollapse_proved': False,
        'interacting_relativistic_field_or_mass_gap_proved': False,
    })
