"""Actual port-exchange responses on the persistent geometric matter object.

The source input is the A5-covariant matching code and its two-carrier3+3'
sector from native_gauge_object. No desired gauge-current matrix or empirical
mass enters. The complete native port exchanges are those of
native_virtual_exchange, on all10626 binary load-four configurations.

We derive every second port-weight response before imposing the independently
specified requirement that the actual geometric six-dimensional object remain
invariant. Equal-color carriers leave only a scalar preserving response.
Opposite colors leave six Hermitian responses; their Lie closure is tested
without identifying encoded occupation with electric charge. A compression
which leaks out of the object is not counted as a persistent interaction.

Source controls, corrected matching-code protection and a coherent response
parameter remain declared. These finite response and symmetry calculations
do not derive a physical time, interacting3+1D law or particle mass spectrum.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
import sympy as s

from tools import strict_json
from . import native_gauge_object as objects
from . import native_virtual_exchange as virtual

ROOT = Path(__file__).resolve().parents[2]
PORT_PAIRS = tuple(itertools.combinations_with_replacement(range(12), 2))


def _entries(matrix):
    return [[str(x) for x in matrix.row(i)] for i in range(matrix.rows)]


def _rank(matrix):
    return matrix.to_DM(extension=True).rank()


def _independent(matrices):
    if not matrices:
        return [], ()
    columns = s.Matrix.hstack(*(matrix.reshape(matrix.rows*matrix.cols, 1) for matrix in matrices))
    pivots = columns.to_DM(extension=True).rref()[1]
    return [matrices[i] for i in pivots], pivots


def _normalize(matrix):
    nonzero = next((x for x in matrix if x != 0), None)
    return matrix/nonzero if nonzero is not None else matrix


@lru_cache(maxsize=1)
def geometric_coordinates():
    """Rational coordinates for the actual character-projector range."""
    projector = objects.geometric_projectors()['six']
    if any(not value.is_Rational for value in projector):
        raise ArithmeticError('the combined geometric3+3prime projector lost rationality')
    pivots = projector.to_DM().rref()[1]
    frame = s.ImmutableMatrix(projector[:, list(pivots)])
    metric = s.ImmutableMatrix(frame.T*frame)
    if frame.cols != 6 or _rank(metric) != 6:
        raise ArithmeticError('the actual geometric object does not have dimension six')
    scale = s.ilcm(*(value.q for value in frame))
    pscale = s.ilcm(*(value.q for value in projector))
    return {'projector': projector, 'frame': frame, 'metric': metric,
            'left_inverse': metric.inv()*frame.T,
            'integer_frame': np.asarray(scale*frame, dtype=np.int64),
            'frame_scale': int(scale),
            'integer_projector': np.asarray(pscale*projector, dtype=np.int64),
            'projector_scale': int(pscale)}


def _colors(left, right):
    if type(left) is not int or type(right) is not int or left not in (0, 1) or right not in (0, 1):
        raise ValueError('two declared native matching-color choices are required')
    return left, right


@lru_cache(maxsize=4)
def native_port_columns(left=0, right=0):
    """Exact full-native Q V_p U columns for the selected persistent codes.

    U=u/6 for an integer two-carrier matching encoder u. Each port exchange
    has P V_p P=(5/18)I. Thus Q V_p U=(18V_p u-5u)/108. No intermediate native
    load or syndrome is discarded in computing these columns.
    """
    left, right = _colors(left, right)
    masks, _, permutations, _, _ = virtual.geometry()
    index = {mask: i for i, mask in enumerate(masks)}
    left_supports, right_supports = objects.equivariant_matchings(left), objects.equivariant_matchings(right)
    rows = np.array([[index[sum(1 << p for p in a) | (sum(1 << p for p in b) << 12)]
                      for a in lrow for b in rrow]
                     for lrow in left_supports for rrow in right_supports], dtype=np.int64)
    u = sparse.csr_matrix((np.ones(rows.size, dtype=np.int64),
                          (rows.ravel(), np.repeat(np.arange(36), 36))), shape=(10626, 36))
    columns = []
    for permutation in permutations:
        native = u-u[permutation]
        if not np.array_equal((u.T@native).toarray(), 10*np.eye(36, dtype=np.int64)):
            raise ArithmeticError('a first-order native port compression is not scalar')
        columns.append(18*native-5*u)
    return u, tuple(columns)


@lru_cache(maxsize=4)
def local_code_projectors(left=0, right=0):
    left, right = _colors(left, right)
    masks = virtual.geometry()[0]
    index = {mask: i for i, mask in enumerate(masks)}
    native_two = tuple(sum(1 << p for p in pair) for pair in itertools.combinations(range(12), 2))
    supports = [[sum(1 << p for p in pair) for pair in matching]
                for matching in objects.equivariant_matchings(left)]
    rights = [[sum(1 << p for p in pair) for pair in matching]
              for matching in objects.equivariant_matchings(right)]
    ga = np.array([[index[a | b << 12] for a in row] for row in supports for b in native_two])
    gb = np.array([[index[a | b << 12] for b in row] for a in native_two for row in rights])
    return virtual._group_matrix(ga, 10626), virtual._group_matrix(gb, 10626)


@lru_cache(maxsize=4)
def native_interaction(left=0, right=0):
    _, columns = native_port_columns(left, right)
    total = sum(columns)
    return s.ImmutableMatrix(s.Matrix((total.T@total).toarray())/11664)


@lru_cache(maxsize=8)
def response_matrices(left=0, right=0, protection='flat_joint_phase'):
    """All78 exact mixed responses, full-color leakage and actual M6 action.

    For additive local protection the virtual inverse has eigenvalues1/U
    and1/(2U), on one- and two-code-complement sectors respectively. On Q,
    at U=1 it is(6I+Pa_integer+Pb_integer)/12. The flat joint-phase law uses
    the unweighted Gram and has its own finite commutator-word construction.
    """
    _, columns = native_port_columns(left, right)
    if protection not in ('flat_joint_phase', 'additive_local_penalty'):
        raise ValueError('a declared virtual response law is required')
    if protection == 'additive_local_penalty':
        pa, pb = local_code_projectors(left, right)
        weighted = tuple(6*column+pa@column+pb@column for column in columns)
        denominator = 139968
    else:
        weighted = columns
        denominator = 11664
    coordinates = geometric_coordinates()
    li, scale = coordinates['integer_frame'], coordinates['frame_scale']
    gi = coordinates['metric'].inv()
    complement = (coordinates['projector_scale']*np.eye(36, dtype=np.int64)
                  -coordinates['integer_projector'])
    full, reduced, leakage = [], [], []
    for p, q in PORT_PAIRS:
        numerator = (columns[p].T@weighted[q]).toarray()
        if not np.array_equal(numerator, numerator.T):
            raise ArithmeticError('commuting native exchanges lost their symmetric mixed response')
        full.append(s.ImmutableMatrix(s.Matrix(numerator)/denominator))
        reduced.append(s.ImmutableMatrix(gi*s.Matrix(li.T@numerator@li)/(denominator*scale*scale)))
        leakage.append((complement@numerator@li).reshape(216))
    constraints = s.ImmutableMatrix(np.stack(leakage, axis=1))
    return tuple(full), tuple(reduced), constraints


def _kernel(matrix):
    rows = matrix.to_DM().to_field().nullspace(divide_last=True).to_Matrix()
    return s.ImmutableMatrix.vstack(*(_normalize(rows.row(i)) for i in range(rows.rows))) if rows.rows else s.zeros(0, matrix.cols)


@lru_cache(maxsize=8)
def preserving_responses(left=0, right=0, protection='flat_joint_phase'):
    """Solve invariance in all36 color outputs before taking the response span."""
    full, reduced, constraints = response_matrices(left, right, protection)
    coefficients = _kernel(constraints)
    response_map = s.Matrix.hstack(*(matrix.reshape(36, 1) for matrix in reduced))
    responses = response_map*coefficients.T
    pivots = responses.to_DM().rref()[1]
    selected = []
    source_coefficients = []
    for pivot in pivots:
        matrix = s.Matrix(6, 6, responses[:, pivot])
        nonzero = next(x for x in matrix if x != 0)
        selected.append(s.ImmutableMatrix(matrix/nonzero))
        source_coefficients.append(coefficients.row(pivot)/nonzero)
    frame = geometric_coordinates()['frame']
    projector = geometric_coordinates()['projector']
    for coefficients_row, response in zip(source_coefficients, selected, strict=True):
        source = sum((a*h for a, h in zip(coefficients_row, full, strict=True)), s.zeros(36))
        if source*frame != frame*response or (s.eye(36)-projector)*source*frame != s.zeros(36, 6):
            raise ArithmeticError('a selected response does not preserve the complete object')
    return {'constraint_rank': _rank(constraints),
            'coefficient_kernel_dimension': coefficients.rows,
            'response_dimension': len(selected), 'basis': tuple(selected),
            'source_coefficients': tuple(s.ImmutableMatrix(row) for row in source_coefficients)}


def lie_closure(hermitian, metric):
    """Exact graded real Lie closure in the actual positive Gram metric.

    Every element is K+iS, with K real G-skew and S real G-selfadjoint.
    Keeping the two grades separate prevents a complex-rank substitution
    from changing the meaning of real Lie dimension.
    """
    symmetric, _ = _independent([_normalize(s.Matrix(h)) for h in hermitian])
    skew = []
    for matrix in symmetric:
        if matrix.T*metric != metric*matrix:
            raise ArithmeticError('source phase response is not metric self-adjoint')
    for _ in range(36):
        previous = len(symmetric)+len(skew)
        odd = skew+[a*b-b*a for i, a in enumerate(symmetric) for b in symmetric[:i]]
        odd += [a*b-b*a for i, a in enumerate(skew) for b in skew[:i]]
        skew, _ = _independent([_normalize(matrix) for matrix in odd if any(matrix)])
        even = symmetric+[a*b-b*a for a in skew for b in symmetric]
        symmetric, _ = _independent([_normalize(matrix) for matrix in even if any(matrix)])
        if len(symmetric)+len(skew) == previous:
            break
    else:
        raise ArithmeticError('finite six-dimensional Lie closure did not stabilize')
    if any(matrix.T*metric != -metric*matrix for matrix in skew):
        raise ArithmeticError('source ordered response is not metric skew')
    if any(matrix.T*metric != metric*matrix for matrix in symmetric):
        raise ArithmeticError('source symmetric response escaped its metric domain')
    return tuple(skew), tuple(symmetric)


def response_report(left, right, protection='flat_joint_phase'):
    coordinates = geometric_coordinates()
    frame, metric, inverse = coordinates['frame'], coordinates['metric'], coordinates['left_inverse']
    b = native_interaction(left, right)
    if protection == 'additive_local_penalty':
        b = b/2
    full, reduced, constraints = response_matrices(left, right, protection)
    block = s.ImmutableMatrix(inverse*b*frame)
    if b*frame != frame*block:
        raise ArithmeticError('uniform source interaction fails to preserve the actual object')
    preserving = preserving_responses(left, right, protection)
    skew, symmetric = lie_closure(preserving['basis'], metric)
    projectors = objects.geometric_projectors()
    three = inverse*projectors['three']*frame
    three_prime = s.eye(6)-three
    phase_ranks, offblock = [], []
    for projector in (three, three_prime):
        traceless = [(projector*h*projector-s.trace(projector*h)*projector/3).applyfunc(s.expand)
                     for h in preserving['basis']]
        phase_ranks.append(_rank(s.Matrix.hstack(*(h.reshape(36, 1) for h in traceless))))
    for h in preserving['basis']:
        offblock.append(any((h*three-three*h).applyfunc(s.expand)))
    # The scalar phase is invisible to conjugation readout on End(M6).
    currents = [s.I*h for h in symmetric]+list(skew)
    tracefree = [h-s.trace(h)*s.eye(6)/6 for h in currents]
    observable_rank = _rank(s.Matrix.hstack(*(h.reshape(36, 1) for h in tracefree)))
    derivatives = []
    for port in range(12):
        derivative = s.zeros(6)
        for (p, q), h in zip(PORT_PAIRS, reduced, strict=True):
            if p == q == port:
                derivative += 2*h
            elif port in (p, q):
                derivative += 2*h
        derivatives.append(s.ImmutableMatrix(derivative))
    report = {'left_matching_choice': left, 'right_matching_choice': right,
              'virtual_response_law': protection,
              'response_definition': ('P V_p Q V_q P, from joint-code-phase commutator polarization'
                 if protection == 'flat_joint_phase' else
                 'P V_p Q (Q H0 Q)^-1 Q V_q P, H0=(I-P_A)+(I-P_B), U=1'),
              'uniform_interaction_on_geometric_M6_exact': _entries(block),
              'uniform_interaction_preserves_the_object': True,
              'first_port_derivative_span': _rank(s.Matrix.hstack(*(h.reshape(36, 1) for h in derivatives))),
              'raw_mixed_response_count': 78,
              'raw_compressed_response_span': _rank(s.Matrix.hstack(*(h.reshape(36, 1) for h in reduced))),
              'all_36_color_output_leakage_constraint_rank': preserving['constraint_rank'],
              'preserving_coefficient_kernel_dimension': preserving['coefficient_kernel_dimension'],
              'preserving_Hermitian_response_dimension': preserving['response_dimension'],
              'preserving_response_basis_exact': [_entries(h) for h in preserving['basis']],
              'native_mixed_response_coefficients_exact': [_entries(row)[0] for row in preserving['source_coefficients']],
              'Lie_real_skew_dimension': len(skew), 'Lie_imaginary_selfadjoint_dimension': len(symmetric),
              'Lie_real_dimension': len(skew)+len(symmetric),
              'commutator_observable_dimension': observable_rank,
              'commutator_observable_algebra': 'End(geometric_M6), with conjugation readout',
              'common_phase_is_silent_for_every_larger_reference_algebra': False,
              'geometric_triplet_phase_traceless_ranks': phase_ranks,
              'preserving_phase_response_couples_distinct_triplets': any(offblock),
              'physical_gauge_group_or_local_redundancy_derived': False,
              'mixed_response_linear_combinations_are_declared_coherent_controls': True,
              'finite_native_pulse_errors_are_zero': False}
    prefix = f'c{left}{right}_{protection}'
    arrays = {prefix+'_full_interaction': np.asarray(b, float),
              prefix+'_mixed_responses': np.array(reduced, dtype=np.float64),
              prefix+'_preserving_responses': np.array(preserving['basis'], dtype=np.float64),
              prefix+'_ordered_current_basis': np.array(skew, dtype=np.float64).reshape((-1, 6, 6)),
              prefix+'_phase_current_basis': np.array(symmetric, dtype=np.float64).reshape((-1, 6, 6))}
    return report, arrays


def virtual_domain_control(left=0, right=1):
    _, columns = native_port_columns(left, right)
    pa, pb = local_code_projectors(left, right)
    counts = [(pa@column).nnz+(pb@column).nnz for column in columns]
    total = sum(columns)
    if not all(counts) or (pa@total).nnz or (pb@total).nnz:
        raise ArithmeticError('single-complement terms or their uniform cancellation changed')
    return {'individual_port_single_complement_nonzero_entries': counts,
            'uniform_sum_single_complement_entries': 0,
            'uniform_inverse_denominator': '2U',
            'general_additive_inverse_denominators': ['U', '2U'],
            'flat_Gram_menu_is_an_additive_penalty_derivative_without_correction': False}


def rotation_endogeneity():
    """Actual finite port rotations lie in the ideal connected response group.

    No finite source-pulse implementation follows from this exact Cayley
    identity alone. A near-fivefold class sum is a separate native operator,
    whose promotion to a phase control adds a relative central direction.
    """
    coordinates = geometric_coordinates()
    frame, inverse, metric = coordinates['frame'], coordinates['left_inverse'], coordinates['metric']
    preserved = preserving_responses(0, 1)
    skew, phase = lie_closure(preserved['basis'], metric)
    odd_span = s.Matrix.hstack(*(matrix.reshape(36, 1) for matrix in skew))
    phase_span = s.Matrix.hstack(*(matrix.reshape(36, 1) for matrix in phase))
    rows = objects.geometric_projectors()['rotation_rows']
    actions, triples, class_sum = [], [], s.zeros(6)
    for row in rows:
        rotation = objects.permutation_matrix(row['axis_action'])
        full = s.kronecker_product(rotation, rotation)
        action = s.ImmutableMatrix(inverse*full*frame)
        if full*frame != frame*action or action.T*metric*action != metric:
            raise ArithmeticError('actual port rotation lost its geometric intertwiner')
        actions.append(action)
        if row['order'] == 3:
            cayley = (action-s.eye(6))*(action+s.eye(6)).inv()
            if _rank(odd_span.row_join(cayley.reshape(36, 1))) != len(skew):
                raise ArithmeticError('an actual order-three rotation is outside the native current algebra')
            if (s.eye(6)+cayley)*(s.eye(6)-cayley).inv() != action:
                raise ArithmeticError('native Cayley rotation reconstruction failed')
            triples.append(action)
        if row['order'] == 5 and s.sympify(row['character3']) > 0:
            class_sum += action
    generated = {tuple(s.eye(6))} | {tuple(action) for action in triples}
    generated |= {tuple(first*second) for first in triples for second in triples}
    if any(tuple(action) not in generated for action in actions):
        raise ArithmeticError('two order-three source rotations do not cover the finite group')
    relative = class_sum-s.trace(class_sum)*s.eye(6)/6
    if relative.T*metric != metric*relative or relative == s.zeros(6):
        raise ArithmeticError('source fivefold class sum has no relative phase')
    added_rank = _rank(phase_span.row_join(relative.reshape(36, 1)))
    return {'actual_native_rotation_count': len(actions), 'order_three_rotation_count': len(triples),
            'every_rotation_is_at_most_two_order_three_rotations': True,
            'order_three_Cayley_generators_belong_to_actual_real_skew_response_span': True,
            'finite_source_pulse_and_physical_error_certificate_supplied_by_this_identity': False,
            'near_fivefold_class_sum_exact': _entries(class_sum),
            'near_class_sum_triplet_values_exact': ['2+2*sqrt(5)', '2-2*sqrt(5)'],
            'relative_class_sum_phase_in_original_current_span': added_rank == len(phase),
            'current_dimension_after_promoting_class_sum_to_phase_control': len(skew)+added_rank,
            'class_sum_phase_is_a_separate_control_choice': True}


def collective_stabilizer(interaction, metric):
    """Global collective dynamical symmetry, not a local gauge redundancy."""
    n = metric.rows
    if interaction.shape != (n*n, n*n) or metric.shape != (n, n):
        raise ValueError('an interaction on two copies of the stated metric space is required')
    product_metric = s.kronecker_product(metric, metric)
    if interaction.T*product_metric != product_metric*interaction:
        raise ValueError('the interaction is not self-adjoint in its declared metric')
    result = {}
    for symmetric in (False, True):
        basis = []
        for i in range(n):
            for j in range(i if symmetric else i+1, n):
                value = s.zeros(n); value[i, j] = 1
                if i != j:
                    value[j, i] = 1 if symmetric else -1
                basis.append(metric.inv()*value)
        columns = []
        for value in basis:
            collective = s.kronecker_product(value, s.eye(n))+s.kronecker_product(s.eye(n), value)
            columns.append((interaction*collective-collective*interaction).reshape(n**4, 1))
        kernel = _kernel(s.Matrix.hstack(*columns))
        operators = []
        for row in range(kernel.rows):
            value = sum((a*b for a, b in zip(kernel.row(row), basis, strict=True)), s.zeros(n))
            collective = s.kronecker_product(value, s.eye(n))+s.kronecker_product(s.eye(n), value)
            if interaction*collective != collective*interaction:
                raise ArithmeticError('collective stabilizer kernel failed direct substitution')
            operators.append(_entries(_normalize(value)))
        result['imaginary_symmetric' if symmetric else 'real_skew'] = {
            'dimension': len(operators), 'basis_exact': operators}
    result['local_gauge_redundancy_derived'] = False
    return result


@lru_cache(maxsize=1)
def composite_interaction():
    """Four actual matching carriers, with full color-space leakage retained."""
    coordinates = geometric_coordinates()
    li, scale = coordinates['integer_frame'], coordinates['frame_scale']
    metric = s.kronecker_product(coordinates['metric'], coordinates['metric'])
    colors = (0, 1, 0, 1)
    cross_pairs = ((0, 2), (0, 3), (1, 2), (1, 3))
    interactions = [native_interaction(colors[a], colors[c]) for a, c in cross_pairs]
    bscale = s.ilcm(*(value.q for b in interactions for value in b))
    encoded = np.kron(li, li)
    responses, raw = [], []
    for (a, c), b in zip(cross_pairs, interactions, strict=True):
        bi = np.asarray(bscale*b, dtype=np.int64)
        order = [a, c]+[i for i in range(4) if i not in (a, c)]+[4]
        front = np.transpose(encoded.reshape(6, 6, 6, 6, 36), order).reshape(36, -1)
        moved = (bi@front).reshape(6, 6, 6, 6, 36)
        back = np.transpose(moved, np.argsort(order)).reshape(1296, 36)
        compressed = s.Matrix(encoded.T@back)/(bscale*scale**4)
        responses.append(s.ImmutableMatrix(metric.inv()*compressed))
        raw.append(back)
    if responses[0] != responses[3] or responses[1] != responses[2]:
        raise ArithmeticError('source cross-pair types lost their exact exchange symmetry')
    total = sum(responses, s.zeros(36))
    moved = sum(raw)
    leakage = s.Matrix(moved.T@moved)/(bscale*bscale*scale**4)-total.T*metric*total
    witness = next((i for i in range(36) if leakage[i, i] != 0), None)
    if witness is None:
        raise ArithmeticError('the four-carrier interaction unexpectedly preserves both objects')
    norm = s.factor(leakage[witness, witness]/metric[witness, witness])
    return {'matching_color_pattern': list(colors),
            'full_four_carrier_color_dimension': 1296, 'compressed_object_pair_dimension': 36,
            'equal_cross_pair_types_have_equal_compressions': True,
            'full_action_preserves_geometric_object_pair': False,
            'normalized_input_column': witness, 'leakage_norm_squared_exact': str(norm),
            'source_interaction_is_a_local_gauge_law': False}, s.ImmutableMatrix(total)


def full_interface_scope():
    """The restricted common phase has no universal unobservability claim."""
    full, reduced, constraints = response_matrices(0, 1)
    coefficients = _kernel(constraints)
    full_map = s.Matrix.hstack(*(h.reshape(1296, 1) for h in full))
    preserved_full = full_map*coefficients.T
    full_rank = _rank(preserved_full)
    restricted_rank = preserving_responses(0, 1)['response_dimension']
    full_identity = s.eye(36).reshape(1296, 1)
    contains_identity = _rank(preserved_full.row_join(full_identity)) == full_rank
    uniform = native_interaction(0, 1)
    lifts = [sum((a*h for a, h in zip(row, full, strict=True)), s.zeros(36))
             for row in preserving_responses(0, 1)['source_coefficients']]
    noncommuting = sum(uniform*h != h*uniform for h in lifts)
    return {'full_corrected_color_algebra': 'End(C^36)',
            'full_preserving_Hermitian_response_rank': full_rank,
            'restricted_geometric_M6_response_rank': restricted_rank,
            'full_current_Lie_dimension_lower_bound': 17+full_rank-restricted_rank,
            'full_End36_commutator_dimension_lower_bound': 17+full_rank-restricted_rank-int(contains_identity),
            'global_identity36_belongs_to_preserving_response_span': contains_identity,
            'uniform_B_is_scalar_on_complete_color_space': uniform == uniform[0, 0]*s.eye(36),
            'uniform_B_noncommuting_selected_preserving_lifts': noncommuting,
            'complement_coherence_can_observe_uniform_B_relative_phase': True,
            'uniform_B_is_an_enlarged_interface_central_generator': noncommuting == 0,
            'relative_geometric_triplet_phase_supplied_by_common_M6_phase': False}


def bell_reference():
    """Both native singlet lines are tested against the same actual interaction."""
    _, interaction = composite_interaction()
    coordinates = geometric_coordinates()
    frame, inverse, metric = coordinates['frame'], coordinates['left_inverse'], coordinates['metric']
    p = inverse*objects.geometric_projectors()['three']*frame
    q = s.eye(6)-p
    pair_metric = s.kronecker_product(metric, metric)
    singlets = s.Matrix.hstack((p*metric.inv()).reshape(36, 1),
                               (q*metric.inv()).reshape(36, 1))
    clean = lambda matrix: matrix.applyfunc(lambda x: s.expand(s.radsimp(x)))
    if clean(singlets.T*pair_metric*singlets) != 3*s.eye(2):
        raise ArithmeticError('actual geometric singlet normalization failed')
    block = clean(singlets.T*pair_metric*interaction*singlets/3)
    if clean(interaction*singlets-singlets*block) != s.zeros(36, 2):
        raise ArithmeticError('the projected interaction does not preserve its singlet sector')
    x = s.Symbol('x')
    characteristic = s.Poly(interaction.charpoly(x).as_expr(), x)
    factors = s.factor_list(characteristic.as_expr())[1]
    roots = {}
    for polynomial, multiplicity in factors:
        for value, root_multiplicity in s.roots(polynomial, x).items():
            roots[value] = roots.get(value, 0)+multiplicity*root_multiplicity
    if sum(roots.values()) != 36 or any(value.is_real is not True for value in roots):
        raise ArithmeticError('complete exact projected reference spectrum was not resolved')
    ordered = sorted(roots)
    upper_Bell = max(block.eigenvals())
    gap = ordered[-1]-ordered[-2] if roots[ordered[-1]] == 1 else s.Integer(0)
    skew, phase = lie_closure(preserving_responses(0, 1)['basis'], metric)

    # vec(P G^-1) is the invariant Bell tensor on each real triplet.
    # Every block-preserving G-skew current fixes both tensors. An imaginary
    # G-selfadjoint current fixes a mixed Bell ray exactly when both blocks
    # are the same scalar; for a pure ray only that triplet must be scalar.
    # Verify these hypotheses on the actual source closure before using the
    # smaller exact constraint matrices.
    if any(clean(current*p-p*current) != s.zeros(6) for current in skew+phase):
        raise ArithmeticError('reference-ray test requires actual triplet-preserving currents')
    if any(s.expand(s.trace(p*current)-s.trace(q*current)) != 0 for current in phase):
        raise ArithmeticError('reference-ray test lost equal triplet phase traces')

    def ray_stabilizer(pure=False):
        constraints = [clean(p*current*p-s.trace(p*current)*p/3) if pure
                       else current-s.trace(current)*s.eye(6)/6 for current in phase]
        dimensions = [len(skew), len(phase)-_rank(s.Matrix.hstack(*(h.reshape(36, 1) for h in constraints)))]
        return {'real_skew_dimension': dimensions[0],
                'imaginary_selfadjoint_dimension': dimensions[1],
                'raw_current_dimension': sum(dimensions),
                'End_M6_conjugation_observable_dimension': sum(dimensions)-1}

    upper_Bell_vector = next(vector for value, _, vectors in block.eigenvects()
                             if value == upper_Bell for vector in vectors)
    if any(value == 0 for value in upper_Bell_vector):
        raise ArithmeticError('the actual selected Bell reference is not mixed')
    return {'actual_matching_color_pattern': [0, 1, 0, 1],
            'geometric_Bell_block_exact': _entries(block),
            'geometric_Bell_eigenvalues_exact': [str(value) for value in sorted(block.eigenvals())],
            'pure_triplet_Bell_lines_are_eigenrays': block[0, 1] == 0,
            'full_projected_characteristic_factors': [
                {'polynomial': str(polynomial), 'multiplicity': multiplicity}
                for polynomial, multiplicity in factors],
            'number_of_full_projected_eigenvalues_above_upper_Bell_value': sum(
                multiplicity for value, multiplicity in roots.items() if value > upper_Bell),
            'upper_Bell_is_unique_projected_maximum': ordered[-1] == upper_Bell and roots[upper_Bell] == 1,
            'next_projected_Gram_eigenvalue_exact': str(ordered[-2]),
            'projected_Gram_gap_exact': str(s.simplify(gap)),
            'upper_Bell_ray_coordinates_exact': [str(value) for value in upper_Bell_vector],
            'actual_mixed_reference_ray_stabilizer': ray_stabilizer(),
            'supplied_pure_triplet_reference_control_stabilizer': ray_stabilizer(pure=True),
            'sign_of_second_order_energy_shift': '-t^2/(2U) times the source Gram interaction',
            'full_native_four_carrier_ground_state_selected_by_this_compression': False,
            'a_pure_triplet_reference_supplied_to_select_SM_currents': False}


def build():
    same, same_arrays = response_report(0, 0)
    opposite, opposite_arrays = response_report(0, 1)
    additive_same, additive_same_arrays = response_report(0, 0, 'additive_local_penalty')
    additive_opposite, additive_opposite_arrays = response_report(0, 1, 'additive_local_penalty')
    composite, h = composite_interaction()
    coordinates = geometric_coordinates()
    same_stabilizer = collective_stabilizer(native_interaction(0, 0), s.eye(6))
    composite_stabilizer = collective_stabilizer(h, coordinates['metric'])
    report = {'schema': 'oph.native_persistent_response.v1',
              'object_producer': 'native_gauge_object: actual corrected matching codes and geometric3+3prime',
              'same_color_response': same, 'opposite_color_response': opposite,
              'additive_local_same_color_response': additive_same,
              'additive_local_opposite_color_response': additive_opposite,
              'virtual_denominator_domain_control': virtual_domain_control(),
              'native_rotation_endogeneity': rotation_endogeneity(),
              'full_interface_scope': full_interface_scope(),
              'actual_pair_reference_test': bell_reference(),
              'same_color_global_collective_stabilizer': same_stabilizer,
              'four_carrier_composite_interaction': composite,
              'projected_composite_global_collective_stabilizer': composite_stabilizer,
              'empirical_inputs_to_prediction': [],
              'operational_control_scope': {
                  'flat_joint_Gram_polarization': 'Bpq=(B(ep+eq)-B(ep)-B(eq))/2 for p!=q; Bpp=B(ep)',
                  'negative_Gram_term_requires_inverse_commutator_word': True,
                  'signed_synthesis_requires_declared_joint_code_phases_and_coherent_port_rotations': True,
                  'linear_combination_limit_uses_Trotter_and_finite_word_error_bounds': True,
                  'all_selected_signed_schedules_executed_as_complete_native_histories_here': False,
                  'additive_local_response_uses_its_own_protection_gap_and_virtual_inverse': True,
                  'exact_autonomous_finite_current_law_selected_by_CE_repair': False},
              'scope': {'desired_gauge_current_matrices_imported': False,
                        'source_response_matrices_derived_on_actual_persistent_object': True,
                        'complete_continuing_native_gauge_history_executed': False,
                        'one_triplet_only_phase_restriction_selected_by_repair': False,
                        'Standard_Model_gauge_current_law_selected': False,
                        'measured_particle_spectrum_derived': False}}
    arrays = {**same_arrays, **opposite_arrays, **additive_same_arrays, **additive_opposite_arrays,
              'geometric_rational_frame': np.asarray(coordinates['frame'], float),
              'geometric_Gram_metric': np.asarray(coordinates['metric'], float),
              'four_carrier_projected_interaction': np.asarray(h, float)}
    return report, arrays


def source_bindings():
    paths = set(virtual.source_bindings()) | set(objects.source_bindings()) | {
        'model/common_source/native_persistent_response.py',
        'model/common_source/test_native_persistent_response.py',
        'model/common_source/native_gauge_object.py',
        'model/common_source/test_native_gauge_object.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def _compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual):
        raise ValueError(path+' original type differs')
    if type(expected) is dict:
        if set(expected) != set(actual):
            raise ValueError(path+' keys differ')
        for key in expected:
            _compare(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError(path+' length differs')
        for i, (a, b) in enumerate(zip(expected, actual, strict=True)):
            _compare(a, b, path+f'[{i}]')
    elif expected != actual:
        raise ValueError(path+' differs')


def check(path):
    path = Path(path)
    pins = source_bindings()
    report = strict_json.load(path)
    if type(report) is not dict or report.get('source_bindings') != pins:
        raise ValueError('persistent-response source bindings differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('persistent-response binary custody differs')
    fresh, arrays = build()
    if source_bindings() != pins:
        raise ValueError('persistent-response sources changed during replay')
    _compare(fresh, {key: value for key, value in report.items() if key not in ('source_bindings', 'artifact')})
    with np.load(artifact, allow_pickle=False) as stored:
        if len(stored.files) != len(set(stored.files)) or set(stored.files) != set(arrays):
            raise ValueError('persistent-response complete array keys differ')
        for key, value in arrays.items():
            old = stored[key]
            if old.shape != value.shape or old.dtype != value.dtype or not np.isfinite(old).all():
                raise ValueError('persistent-response array schema differs: '+key)
            if not np.array_equal(old, value):
                raise ValueError('persistent-response rational array comparison failed: '+key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check', type=Path)
    mode.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('Persistent native object, actual port response and symmetry controls: PASS')
        return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('choose fresh receipt and artifact paths')
    pins = source_bindings()
    report, arrays = build()
    if source_bindings() != pins:
        raise ValueError('persistent-response sources changed during generation')
    report['source_bindings'] = pins
    with artifact.open('xb') as output:
        np.savez_compressed(output, **arrays)
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x') as output:
        json.dump(report, output, indent=2, sort_keys=True, allow_nan=False)
        output.write('\n')


if __name__ == '__main__':
    main()
