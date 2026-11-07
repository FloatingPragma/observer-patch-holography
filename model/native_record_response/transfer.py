"""Positive transfer reconstructed from correctly reversed native records.

The Hilbert space here represents scalar stationary record histories. It is
not the full quantum operator-history space or an implemented real-time field.
"""
from __future__ import annotations

from fractions import Fraction as F

import sympy as s

from . import response as owning
from ..protected_reference import reversal


LOWER = F(296, 495)
DEFICIT = 1-LOWER
SOURCE_FILES = (
    'model/native_record_response/transfer.py',
    'model/native_record_response/TRANSFER.md',
    'model/native_record_response/test_transfer.py',
    'model/protected_reference/__init__.py',
    'model/protected_reference/reference.py',
    'model/protected_reference/reversal.py',
    'model/protected_reference/REVERSAL.md',
    'model/exchange_sector/statistics.py',
    'model/primitive_observer/DERIVATION_TRANSFER.md',
    'model/native_quantum_reconstruction/PROOF.md',
    'Lean/SourceLaw/ReflectionPositivity.lean',
    'Lean/SourceLaw/TransferStructure.lean',
)


def label(row):
    return row.kind, row.law, row.seam, row.coin, row.raw


def reverse_label(value):
    """Reuse the owning 720-slot involution and retain the zero global hold."""
    if type(value) is not tuple or len(value) != 5:
        raise ValueError('one complete native raw label is required')
    kind, law, seam, coin, raw = value
    if value == ('hold', None, None, None, None):
        return value
    if kind != 'repair':
        raise ValueError('stationary completed native repair or hold label required')
    return ('repair', *reversal.reverse_label((law, seam, coin, raw)))


def reverse_word(word):
    if type(word) is not tuple:
        raise ValueError('a finite raw-record tuple is required')
    return tuple(reverse_label(row) for row in reversed(word))


def _index(data):
    result = {label(row): row for row in data['branches']}
    if len(result) != 721:
        raise ArithmeticError('all 721 native raw labels must be retained')
    return result


def forward_word(data, word, vector=None):
    """Unnormalized population vector after every raw event in a word."""
    indexed = _index(data); reverse_word(word)
    result = data['pi'] if vector is None else s.ImmutableMatrix(vector)
    if result.shape != (66, 1):
        raise ValueError('one population-space column is required')
    for value in word:
        row = indexed[value]; following = s.zeros(66, 1)
        for i, j in row.pairs:
            following[j] += row.weight*result[i]
        result = s.ImmutableMatrix(following)
    return result


def conditional_word(data, word):
    """Probability of an entire future word conditional on its initial state."""
    indexed = _index(data); reverse_word(word)
    result = s.ones(66, 1)
    for value in reversed(word):
        row = indexed[value]; following = s.zeros(66, 1)
        for i, j in row.pairs:
            following[i] = row.weight*result[j]
        result = following
    return s.ImmutableMatrix(result)


def reflected_joint(data, left, right, gap=0):
    """Actual stationary word probability across n unread bridge attempts."""
    owning._count(gap)
    state = forward_word(data, reverse_word(left))
    state = data['P']**gap*state
    return sum(forward_word(data, right, state))


def reflected_gram(data, words, gap=0):
    """Indicator-cylinder Gram, keeping reserved and zero words as zero rows."""
    owning._count(gap)
    if type(words) is not tuple:
        raise ValueError('a tuple of finite record words is required')
    functions = [conditional_word(data, word) for word in words]
    bridge = s.diag(*data['pi'])*(data['P'].T**gap)
    return s.ImmutableMatrix([[sum(f.conjugate().multiply_elementwise(bridge*g))
                               for g in functions] for f in functions])


def branch_reversal_certificate(data):
    indexed = _index(data); pi = data['pi']
    reverse_insertions = [s.zeros(66) for _ in range(30)]
    checked = moving = zero = abort = 0
    for row in data['branches']:
        key = label(row); reverse = reverse_label(key)
        if reverse not in indexed or reverse_label(reverse) != key:
            raise ArithmeticError('native raw reversal is not a complete involution')
        other = indexed[reverse]
        if other.pairs != tuple(sorted((j, i) for i, j in row.pairs)):
            raise ArithmeticError('reversed native record has a different partial map')
        zero += not row.pairs
        abort += row.kind == 'repair' and row.label is None
        moving += any(i != j for i, j in row.pairs)
        for i, j in row.pairs:
            if row.weight*pi[i] != other.weight*pi[j]:
                raise ArithmeticError('weighted raw-history adjoint identity failed')
            if row.law == 'E':
                e = owning.local.SEAMS.index(row.seam)
                if other.label != row.label-2*row.current:
                    raise ArithmeticError('contrast reflection lost its time orientation')
                reverse_insertions[e][j, i] += row.weight*other.feature(row.seam)
            checked += 1
    for e, reverse in enumerate(reverse_insertions):
        for i in range(66):
            for j in range(66):
                if reverse[j, i]*pi[i] != data['B'][e][i, j]*pi[j]:
                    raise ArithmeticError('reflected insertion is not the weighted adjoint')
        if (s.ones(1, 66)*reverse != data['G'][e, :]
                or reverse*pi != s.diag(*pi)*data['D'][e, :].T):
            raise ArithmeticError('reflected incoming/outgoing response identity failed')
    return {'complete_raw_slots': len(indexed), 'zero_slots_including_hold': zero,
            'reserved_abort_slots': abort, 'moving_slots': moving,
            'weighted_source_columns_checked': checked,
            'same_owning_raw_label_involution': True,
            'forward_adjoint': 'T_reverse = Pi T_h^T Pi^-1',
            'backward_adjoint_on_L2_pi': 'L_reverse = L_h^dagger_pi',
            'reflected_contrast': 'incoming contrast minus twice the recorded current',
            'all_thirty_reflected_insertions_verified': True}


def reflected_response(data, gap=0):
    """Centered true-reflection contrast Gram with n unread bridge attempts."""
    owning._count(gap)
    return s.ImmutableMatrix(data['D']*data['P']**gap*s.diag(*data['pi'])*data['D'].T
                             -data['mean']*data['mean'].T)


def logarithm_tail(degree):
    """Operator-norm bound in the reconstructed Hilbert space, no time fit."""
    owning._count(degree, 'logarithm polynomial degree')
    return DEFICIT**(degree+1)/((degree+1)*(1-DEFICIT))


def logarithm_action(data, function, degree):
    """Rational action of sum_(k=1)^m (I-Q)^k/k, Q=P^T."""
    owning._count(degree, 'logarithm polynomial degree')
    term = s.ImmutableMatrix(function)
    if term.shape != (66, 1):
        raise ValueError('a function on all 66 source configurations is required')
    result = s.zeros(66, 1); deficit = s.eye(66)-data['P'].T
    for k in range(1, degree+1):
        term = deficit*term
        result += term/k
    return s.ImmutableMatrix(result)


def clock_distinction(data):
    """An actual transition separates transfer H from working-state -log(pi)."""
    x,y = (data['states'].index(tuple(int(i in pair) for i in range(12)))
           for pair in ((0,1),(1,5)))
    pi = data['pi'];p = data['P'][y,x]
    entry = s.simplify(p*s.sqrt(pi[x]/pi[y]))
    commutator = s.simplify(entry*(pi[x]-pi[y]))
    if (p != F(4,495) or entry != s.sqrt(2)/99
            or commutator != 85*s.sqrt(2)/148599):
        raise ArithmeticError('the exact source clock-distinction witness changed')
    return {'input_configuration':[0,1], 'output_configuration':[1,5],
            'source_indices':[x,y], 'source_weights':[str(pi[x]),str(pi[y])],
            'P_output_input':str(p), 'S_output_input':'sqrt(2)/99',
            'S_Pi_commutator_output_input':'85*sqrt(2)/148599',
            'P_Pi_commutator_output_input':str(p*(pi[x]-pi[y])),
            'S_native_modular_K_commutator_output_input':'sqrt(2)*log(8/25)/99',
            'native_modular_K':'-log(diag(pi))',
            'transfer_H_commutes_with_native_modular_K':False,
            'transfer_H_is_scalar_rescaling_and_shift_of_native_modular_K':False}


def certificate(data=None, ranks=None):
    data = owning.data() if data is None else data
    owning.validate(data)
    rows = branch_reversal_certificate(data)
    minimum = min(data['P'][i, i] for i in range(66))
    if minimum != F(791, 990) or 2*minimum-1 != LOWER:
        raise ArithmeticError('the exact positive-transfer lower bound changed')
    ranks = owning.response_ranks(data) if ranks is None else ranks
    if ranks['observable_rank_by_depth'][-1] != 66:
        raise ArithmeticError('actual record cylinders fail to span the transfer Hilbert space')
    naive = owning.covariance(data, 1)
    vector = s.Matrix([-1, 1]+[0]*28)
    negative = (vector.T*naive*vector)[0]
    if negative != -F(83481968, 1858725825):
        raise ArithmeticError('the naive reflection counterexample changed')
    reflected = reflected_response(data)
    positive = (vector.T*reflected*vector)[0]
    if reflected != reflected.T or positive <= 0:
        raise ArithmeticError('true reflected record Gram failed its positive control')
    return {'schema': 'oph.native_record_transfer.v1',
            'raw_reversal': rows,
            'Hilbert_space': 'L2(pi) of scalar stationary raw-record cylinders modulo reflection nulls',
            'dimension': 66, 'ground_function': '1',
            'backward_transfer': 'Q=P^T',
            'standard_transfer': 'S=Pi^-1/2 P Pi^1/2',
            'spectral_interval': [str(LOWER), '1'],
            'positive_generator': 'H=-log(S); uniquely defined positive self-adjoint transfer logarithm',
            'ground_vector_in_standard_coordinates': 'sqrt(pi)',
            'full_history_reflection_form': '<f_F,Q^n f_F>_pi for every finite scalar cylinder and n>=0',
            'actual_record_cylinders_span_dimension': 66,
            'reflected_contrast_response': 'D P^n diag(pi) D^T - m m^T',
            'reflected_residues': 'positive semidefinite rank-one residues on all 65 nonstationary simple modes',
            'naive_same_orientation_control': {
                'score': 'centered Y_(0,5) minus Y_(0,1)',
                'lag_in_attempts': 1, 'unreflected_quadratic_form': str(negative),
                'properly_reflected_quadratic_form': str(positive)},
            'positive_logarithm_polynomial_tail': {'degree': 16, 'operator_norm_upper': str(logarithm_tail(16)),
                'norm': 'L2(pi) or equivalent standard Hilbert norm'},
            'working_state_clock_distinction': clock_distinction(data),
            'scope': {'same_stationary_scalar_complete_record_law': True,
                'raw_attempts_operational_only_at_integer_spacings': True,
                'full_quantum_operator_history_positivity': False,
                'native_quantum_working_state_equals_reconstructed_ground_ket': False,
                'every_fractional_transfer_power_is_a_native_stochastic_channel': False,
                'unitary_analytic_continuation_is_an_executed_native_service': False,
                'universal_operational_effect_menu_or_Born_bridge_completed': False,
                'physical_clock_spatial_field_or_particle_masses_identified': False}}
