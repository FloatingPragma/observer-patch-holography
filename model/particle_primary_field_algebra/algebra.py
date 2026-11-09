"""Exact premises for a common filtered source algebra and positive-energy GNS.

The limiting subsequence is an existence construction. The finite descriptor
and norm controls do not execute its large spatial density matrices.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import factorial

import sympy as s

from model.particle_common_contract import contract
from model.particle_joint_reconstruction import positive_energy
from model.particle_physical_boundary_state import boundary
from model.particle_primary_locality import locality
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_physical_boundary_state/boundary.py',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_primary_locality/locality.py',
    'model/particle_primary_locality/PROOF.md',
    'model/particle_joint_reconstruction/positive_energy.py',
    'model/particle_joint_reconstruction/POSITIVE_ENERGY.md',
    'model/particle_ground_spectral_records/PROOF.md',
    'model/continuum_field_net/README.md',
    'Lean/QFT/ColimitSelectedStateGNS.lean',
    'Lean/QFT/ColimitGNSRepresentation.lean',
    'Lean/QFT/GroundRecordTail.lean',
)


def _rational(value, name):
    if isinstance(value, (float, bool)):
        raise ValueError(name+' must be an exact rational')
    try:
        return F(value)
    except (TypeError, ValueError, ZeroDivisionError) as error:
        raise ValueError(name+' must be an exact rational') from error


def probe(member, position=('0', '0', '0'), half_width='1/100000',
          kind='charge', axis=0, charge=1, axes=(0, 1)):
    """One actual bounded physical observable section at a centered location.

    The family contains all rational positions, positive rational filter
    widths, integer charges, coordinate axes and coordinate plaquettes.
    Early levels outside the interior and absent charge sectors are zero.
    """
    contract.verify_member(member)
    if not isinstance(position, (tuple, list)) or len(position) != 3:
        raise ValueError('three exact centered coordinates are required')
    position = tuple(_rational(x, 'centered coordinate') for x in position)
    T = ground._q(half_width, 'positive filter half width', True)
    if kind not in ('charge', 'plaquette'):
        raise ValueError('the concrete menu contains charge and plaquette probes')
    if type(axis) is not int or axis not in range(3) or type(charge) is not int:
        raise ValueError('integer charge and coordinate axis are required')
    if (not isinstance(axes, (tuple, list)) or len(axes) != 2
            or any(type(i) is not int or i not in range(3) for i in axes)
            or axes[0] >= axes[1]):
        raise ValueError('two increasing plaquette axes are required')
    a = F(member['source']['spacing']); side = member['source']['side']
    cell = tuple((F(side, 2)+x/a).__floor__() for x in position)
    inside = all(1 <= x <= side-2 for x in cell)
    menu = boundary.charge_spectrum(member['source']['per_link_word_degree'])
    active = inside and (kind != 'charge' or charge in menu)
    support = (cell,) if active else ()
    if active and kind == 'plaquette':
        support = locality.plaquette(side, cell, axes)['owner_support']
    return contract.serial({
        'member_id': member['member_id'], **member['ids'],
        'kind': kind, 'centered_target': position,
        'centered_embedding': 'a*(owner_cell-side/2)',
        'nearest_lower_owner_cell': cell,
        'actual_centered_anchor': tuple(a*(x-F(side, 2)) for x in cell),
        'half_width': T, 'axis': axis, 'integer_Q6_charge': charge,
        'plaquette_axes': tuple(axes), 'active': active,
        'zero_reason': None if active else ('outside_interior' if not inside else 'absent_charge'),
        'unsmeared_owner_support': support,
        'operator': ('spectral projection of outgoing-link right endpoint Q6'
                     if kind == 'charge' else 'compressed ReTr_M15(U_plaquette)/15'),
        'normalization': 'Q6=6Y' if kind == 'charge' else 'native trace divided by15',
        'unsmeared_norm_upper': F(int(active)),
        'filtered_norm_upper': F(int(active)),
        'sequence_uniform_norm_upper': F(1),
        'filter': 'f_T(s)=(1-abs(s)/T)_+/T',
        'sequence_time_Lipschitz_upper': 2/T,
        'commutes_all_local_Gauss_constraints': True,
        'source_even': True, 'stationary_algebra_includes_idle_reference': False,
        'full_spatial_matrix_executed': False,
        'individual_filtered_observable_is_a_PVM_projection': False,
    })


def verify_probe(member, value):
    try:
        expected = probe(member, value['centered_target'], value['half_width'],
                         value['kind'], value['axis'], value['integer_Q6_charge'],
                         value['plaquette_axes'])
    except (KeyError, TypeError) as error:
        raise ValueError('incomplete physical probe descriptor') from error
    if value != expected:
        raise ValueError('the physical probe differs from its source construction')
    return expected


def word_bound(factors):
    """Norm and common-time Lipschitz cap for any finite star word.

    Each factor is (uniform norm cap, positive filter width). Fixed real
    translations and adjoints leave these bounds unchanged. The empty word
    is the identity. Zero caps require no division by another factor cap.
    """
    rows = [(ground._q(c, 'uniform norm cap'),
             ground._q(t, 'positive filter half width', True)) for c, t in factors]
    norm = F(1)
    for c, _ in rows:
        norm *= c
    lipschitz = F(0)
    for i, (c, t) in enumerate(rows):
        summand = 2*c/t
        for j, (other, _) in enumerate(rows):
            if j != i:
                summand *= other
        lipschitz += summand
    return {'norm_upper': norm, 'time_Lipschitz_upper': lipschitz,
            'degree': len(rows)}


def word_negative_tail(member, factors, threshold='1/4'):
    cap = word_bound(factors)
    result = positive_energy.selected_tail_bound(
        member, threshold, cap['norm_upper']**2)
    return {'word': contract.serial(cap), 'spectral_bound': result,
            'covers_all_finite_star_words_on_one_subsequence': True}


def filtered_locality_bound(member, left, right, relative_time='0'):
    """A rational upper bound for the actual two filtered observables.

    The ambient all-state parent cone is restricted to physical operators.
    Fixed temporal filtering enlarges the observation separation condition.
    """
    verify_probe(member, left); verify_probe(member, right)
    source = locality.selected_source_bound(member)
    t = abs(_rational(relative_time, 'relative source time'))
    T = F(left['half_width'])+F(right['half_width'])+t
    a = F(member['source']['spacing'])
    W = source['family']['weighted_rate_times_spacing_upper']
    if not left['active'] or not right['active']:
        return {'commutator_norm_upper': F(0), 'exact_zero': True,
                'zero_reason': 'one source observable is zero'}
    X, Y = left['unsmeared_owner_support'], right['unsmeared_owner_support']
    distance = min(sum(abs(x-y) for x, y in zip(p, q)) for p in X for q in Y)
    margin = F(distance, 2)*a-2*W*T
    # exp(-z)<=2^-floor(z) for z>=0. No floating underflow is used.
    power = max(0, (margin/a).__floor__())
    bound = min(F(2), F(2*len(X), 2**power)) if margin > 0 else F(2)
    return {'commutator_norm_upper': bound, 'exact_zero': False,
            'graph_owner_distance': distance, 'source_spacing': a,
            'total_temporal_reach': T, 'rational_cone_margin': margin,
            'weighted_rate_times_spacing_upper': W,
            'prefactor_owner_count': len(X), 'exponential_power': power,
            'same_primary_source_clock': True,
            'global_Gauss_projected_H_declared_local': False}


def polynomial_support_remainder(spacing, separation, temporal_reach,
                                 count_coefficient='1', count_power=0,
                                 norm_product='1'):
    """Explicit vanishing cone remainder, including centered-cell rounding.

    Assumes |X_a|<=K*a^-p and actual separation>=d-8a. The latter covers
    two centered one-cell plaquettes. More general supports must prove it.
    """
    a = ground._q(spacing, 'positive spacing', True)
    d = ground._q(separation, 'positive assigned separation', True)
    time = ground._q(temporal_reach, 'nonnegative temporal reach')
    K = ground._q(count_coefficient, 'support coefficient')
    p = ground._i(count_power, 'support polynomial power', 0)
    C = ground._q(norm_product, 'product norm cap')
    margin = d/2-2*F(5232, 7)*time
    if margin <= 0 or a > margin/8:
        raise ValueError('a positive cone margin and spacing<=margin/8 are required')
    # (d-8a)/2-2Wt >= margin/2; exp(x)>=x^(p+1)/(p+1)!.
    error = 2*C*K*factorial(p+1)*a/(margin/2)**(p+1)
    return {'commutator_norm_upper': min(2*C, error),
            'explicit_vanishing_upper': error, 'assigned_margin': margin,
            'actual_margin_lower': margin/2,
            'support_power': p, 'rounding_distance_allowance': 8*a}


def exact_two_level_control():
    """Exact all-word Gram and GNS checks in a separate two-state source.

    H=diag(0,pi), rho=|0><0|, T=1. The triangular filter of the plus
    projection has off-diagonal coefficient 2/pi^2. Rational times0,1/2
    make its two translated matrices exact. This is a theorem control.
    """
    I = s.eye(2); h = s.diag(0, s.pi); rho = s.diag(1, 0)
    P = s.ones(2)/2; r = 4/s.pi**2
    B = s.Matrix([[1, r], [r, 1]])/2
    U = s.diag(1, s.I)
    Bt = U*B*U.H
    words = (I, B, Bt, B*Bt, Bt*B, B.H*Bt*B)
    vectors = s.Matrix.hstack(*(w*s.Matrix([1, 0]) for w in words))
    gram = s.Matrix(len(words), len(words),
                    lambda i, j: s.trace(rho*words[i].H*words[j]))
    if (gram-vectors.H*vectors).applyfunc(s.simplify) != s.zeros(len(words)):
        raise ArithmeticError('the complete polynomial Gram matrix is not positive')
    shifted = tuple(U*w*U.H for w in words)
    shifted_gram = s.Matrix(len(words), len(words),
                            lambda i, j: s.trace(rho*shifted[i].H*shifted[j]))
    if (gram-shifted_gram).applyfunc(s.simplify) != s.zeros(len(words)):
        raise ArithmeticError('common source time shifts changed word moments')
    e01 = s.Matrix([[0, 1], [0, 0]]); e10 = e01.H
    null_left = s.trace(rho*e01.H*e01)
    nonnull_right_product = s.trace(rho*(e01*e10).H*(e01*e10))
    positive_weight = s.simplify((B*s.Matrix([1, 0]))[1]**2)
    return {'H': h, 'rho': rho, 'unfiltered_projection': P, 'filtered': B,
            'translated': Bt, 'word_matrices': words, 'word_gram': gram,
            'gram_factor': vectors, 'GNS_generator': h,
            'canonical_ground_vector': s.Matrix([1, 0]),
            'negative_spectral_weight': s.Integer(0),
            'positive_frequency': s.pi, 'positive_frequency_weight': positive_weight,
            'constant_frequency_weight': s.Rational(1, 4),
            'filter_amplitude': r, 'null_left_ideal_example_norm_squared': null_left,
            'same_element_right_product_norm_squared': nonnull_right_product,
            'control_is_primary_spatial_source': False}


def escaping_charge_control(level, fixed_charges=(-1, 0, 1)):
    """A complete finite PVM can lose every fixed-charge outcome in the limit."""
    level = ground._i(level, 'positive level', 1)
    if any(type(q) is not int for q in fixed_charges):
        raise ValueError('integer fixed charges are required')
    return {'finite_PVM_labels': list(range(-level, level+1)),
            'state_charge': level, 'complete_finite_probability_sum': 1,
            'fixed_charge_probabilities': {str(q): int(q == level) for q in fixed_charges},
            'finite_gap_free_ground': True,
            'limiting_fixed_charge_PVM_completeness_follows': False,
            'control_is_primary_spatial_source': False}


def complete_pointer_control():
    """Two complete noncommuting effect readouts with all four histories.

    Exact square-root Kraus maps execute on a separate two-level ground
    control. The final finite pointer retains its off-diagonal entries.
    """
    I = s.eye(2); X = s.Matrix([[0, 1], [1, 0]])
    E = I/2+X/4
    a = (s.sqrt(3)+1)/4; b = (s.sqrt(3)-1)/4
    K = (a*I+b*X, a*I-b*X)
    U = s.diag(1, s.I)
    moved = tuple(U*k*U.H for k in K)
    histories = tuple(second*first for first in K for second in moved)
    completeness = sum((k.H*k for k in histories), s.zeros(2)).applyfunc(s.simplify)
    V = s.Matrix.hstack(*(k*s.Matrix([1, 0]) for k in histories))
    pointer = (V.T*V.conjugate()).applyfunc(s.simplify)
    gate = s.BlockMatrix([[K[0], -K[1]], [K[1], K[0]]]).as_explicit()
    if completeness != I or (gate.H*gate).applyfunc(s.simplify) != s.eye(4):
        raise ArithmeticError('complete finite pointer dilation failed')
    return {'effect': E, 'Kraus': K, 'pointer_gate': gate,
            'history_Kraus': histories, 'completeness': completeness,
            'pointer_density': pointer, 'pointer_Gram_factor': V.conjugate(),
            'all_four_record_probabilities': tuple(pointer[i, i] for i in range(4)),
            'failure_outcome_retained': True, 'postselection': False,
            'control_is_primary_spatial_source': False}


def common_payload(member):
    contract.verify_member(member)
    probes = [probe(member), probe(member, kind='plaquette'),
              probe(member, position=('1', '0', '0'))]
    return contract.serial({
        'member_id': member['member_id'], **member['ids'],
        'construction': 'separable filtered sequence C*-algebra modulo norm-null sequences',
        'actual_probe_examples': probes,
        'word_example': word_negative_tail(member, [(1, '1/100000'), (1, '1/100000')]),
        'filtered_locality_example': filtered_locality_bound(member, probes[0], probes[2]),
        'one_subsequence_for_all_finite_star_moments': True,
        'stationary_state_and_norm_continuous_source_time_group': True,
        'canonical_GNS_U_t_Omega_equals_Omega': True,
        'GNS_source_energy_generator_nonnegative': True,
        'all_fixed_finite_complete_effect_history_probabilities_converge': True,
        'finite_pointer_density_matrices_converge_in_trace_norm': True,
        'growing_canonical_reference_state_limit_claimed': False,
        'compiled_filtered_source_apparatus_execution_claimed': False,
        'limit_exists_for_full_selected_sequence_without_subsequence': False,
        'subsequence_effectively_computed': False,
        'perfect_finite_interlevel_state_embeddings_assumed': False,
        'arbitrary_idle_reference_observables_in_stationary_algebra': False,
        'faithful_nontrivial_field_or_nonzero_energy_limit_proved': False,
        'infinite_charge_PVM_complete_without_tightness': False,
        'Lorentz_symmetry_interactions_or_particles_proved': False,
        'countable_menu_claimed_complete_for_all_physical_observables': False,
        'separate_CAR_continuum_source_transferred': False,
    })


def bind_member(member):
    result = common_payload(member)
    return {'schema': 'oph.primary_filtered_algebra.v1', 'result': result,
            'result_sha256': contract.digest(result)}


def verify_bound(member, receipt):
    if receipt != bind_member(member):
        raise ValueError('the common filtered algebra receipt differs from its source premises')


def certificate():
    control = exact_two_level_control()
    pointer = complete_pointer_control()
    return contract.serial({
        'schema': 'oph.primary_filtered_algebra_certificate.v1',
        'selected_members': [bind_member(contract.member(n)) for n in (1, 6)],
        'two_level_control': {
            'word_count': len(control['word_matrices']),
            'Gram_positive_by_explicit_factorization': True,
            'Gram_rank': int(control['gram_factor'].rank()),
            'positive_frequency': str(control['positive_frequency']),
            'positive_frequency_weight': str(control['positive_frequency_weight']),
            'left_null_ideal_is_not_two_sided': True,
            'is_selected_spatial_source': False,
        },
        'escaping_charge_control': escaping_charge_control(4),
        'complete_pointer_control': {
            'history_count': len(pointer['history_Kraus']),
            'all_record_probabilities': [str(p) for p in pointer['all_four_record_probabilities']],
            'pointer_trace': str(s.trace(pointer['pointer_density'])),
            'complete_isometry_checked': True,
            'off_diagonal_pointer_entries_retained': True,
            'is_selected_spatial_source': False,
        },
        'SOURCE_INPUTS': SOURCE_INPUTS,
        'proof': 'model/particle_primary_field_algebra/PROOF.md',
        'full_interacting_relativistic_field_reconstruction': False,
    })
