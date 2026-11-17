"""Compare the actual exchange response and the declared current on one M6.

The native six Hermitian responses are derived before compression from all
10626 load-four configurations and the full36 color-output invariance test.
Their real Lie closure is transported into the independently constructed
orthonormal geometric3+3prime frame used by native_sm_current. Exact change
of-span identities, rather than names of Lie algebras, establish the
comparison. Fixed real preparations and Pauli-Y readouts distinguish its
five extra phase directions and missing relative central direction.

Readouts here are exact operational comparators. No signed port-Gram word,
source pointer or laboratory current is claimed to have been executed by
this module. The separately compiled coherent apparatus is not evidence that
recorded C/E repair selects either current grammar.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
from pathlib import Path

import sympy as s

from . import native_persistent_response as response
from . import native_sm_current as declared

ROOT = Path(__file__).resolve().parents[2]


def _clean(matrix):
    return s.ImmutableMatrix(matrix.applyfunc(lambda x: s.expand(s.radsimp(x))))


def _zero(matrix):
    return all(x == 0 for x in _clean(matrix))


def _flat(matrix):
    return s.Matrix([s.re(x).expand() for x in matrix]+[s.im(x).expand() for x in matrix])


def _columns(basis):
    return s.Matrix.hstack(*(_flat(matrix) for matrix in basis))


def _rank(matrix):
    return matrix.to_DM(extension=True).rank()


def _coordinates(basis, targets):
    """Exact real coordinates, including a check on every original entry."""
    columns = _columns(basis)
    rows = columns.T.to_DM(extension=True).rref()[1]
    if len(rows) != len(basis):
        raise ArithmeticError('the supplied exact current basis is dependent')
    square = columns.extract(rows, range(len(basis)))
    inverse = square.to_DM(extension=True).to_field().inv().to_Matrix()
    target_columns = _columns(targets)
    result = _clean(inverse*target_columns.extract(rows, range(len(targets))))
    if not _zero(columns*result-target_columns):
        raise ArithmeticError('the supplied exact current span has the wrong action')
    return result


def _unit(i, j, dimension=6):
    result = s.zeros(dimension)
    result[i, j] = 1
    return result


def _su3(offset):
    result = []
    for a in range(3):
        for b in range(a+1, 3):
            e, f = _unit(offset+a, offset+b), _unit(offset+b, offset+a)
            result.extend((e-f, s.I*(e+f)))
    result += [s.I*(_unit(offset, offset)-_unit(offset+1, offset+1)),
               s.I*(_unit(offset+1, offset+1)-_unit(offset+2, offset+2))]
    return tuple(_clean(x) for x in result)


def _quadrupoles():
    return tuple(_clean(x) for x in (
        _unit(3, 3)-_unit(4, 4), _unit(4, 4)-_unit(5, 5),
        _unit(3, 4)+_unit(4, 3), _unit(3, 5)+_unit(5, 3),
        _unit(4, 5)+_unit(5, 4)))


def _commutant(basis):
    """Complex commutant, solved from all entries of the actual source span."""
    units = tuple(_unit(i, j) for i in range(6) for j in range(6))
    constraints = s.Matrix.vstack(*(s.Matrix.hstack(*(
        (g*u-u*g).reshape(36, 1) for u in units)) for g in basis))
    kernel = constraints.to_DM(extension=True).nullspace().to_Matrix()
    return tuple(_clean(s.Matrix(6, 6, kernel.row(i))) for i in range(kernel.rows))


@lru_cache(maxsize=1)
def _exact_data():
    native = response.preserving_responses(0, 1)
    coordinates = response.geometric_coordinates()
    skew, phase = response.lie_closure(native['basis'], coordinates['metric'])
    architecture = declared.exact_architecture()
    T = _clean(architecture['geometric_frame'])
    R = coordinates['frame']
    chart = _clean(coordinates['left_inverse']*T)
    inverse = _clean(T.T*R)
    if (not _zero(R*chart-T) or not _zero(inverse*chart-s.eye(6))
            or not _zero(chart.T*coordinates['metric']*chart-s.eye(6))):
        raise ArithmeticError('native and declared responses do not use the same geometric object')
    H = tuple(_clean(inverse*h*chart) for h in native['basis'])
    G = tuple(_clean(inverse*k*chart) for k in skew)+tuple(
        _clean(s.I*inverse*h*chart) for h in phase)
    D = tuple(_clean(g) for g in architecture['generators'])
    full, _, leakage = response.response_matrices(0, 1)
    # Repeat the full36 invariance identity in the source chart. This binds
    # every emitted Hermitian basis column to the actual78 port-pair rows.
    for coefficients, h in zip(native['source_coefficients'], native['basis'], strict=True):
        lift = sum((a*b for a, b in zip(coefficients, full, strict=True)), s.zeros(36))
        if leakage*coefficients.T != s.zeros(216, 1) or lift*R != R*h:
            raise ArithmeticError('a native Hermitian response lost its complete-output source lift')
    first, second = _su3(0), _su3(3)
    rotations = tuple(second[i] for i in (0, 2, 4))
    p, q = s.diag(1, 1, 1, 0, 0, 0), s.diag(0, 0, 0, 1, 1, 1)
    canonical_native = first+second+(_clean(s.I*s.eye(6)),)
    canonical_declared = first+rotations+(_clean(s.I*p),)
    intersection = first+rotations
    union = first+second+(_clean(s.I*p), _clean(s.I*q))
    native_from_canonical = _coordinates(canonical_native, G)
    declared_from_canonical = _coordinates(canonical_declared, D)
    canonical_from_native = _coordinates(G, canonical_native)
    canonical_from_declared = _coordinates(D, canonical_declared)
    _coordinates(G, intersection)
    _coordinates(D, intersection)
    if (len(skew), len(phase), _rank(_columns(G+D))) != (6, 11, 18):
        raise ArithmeticError('actual native/declared dimensions or union changed')
    if any(not _zero(g+g.conjugate().T) for g in G+D):
        raise ArithmeticError('the shared orthonormal current chart is not anti-Hermitian')
    # The union is already a Lie algebra: all153 brackets are checked in
    # the actual span, not counted as a mere vector-space union.
    brackets = tuple(_clean(a*b-b*a) for i, a in enumerate(union) for b in union[:i])
    _coordinates(union, G+D)
    _coordinates(union, brackets)
    commutant = _commutant(H)
    if len(commutant) != 2:
        raise ArithmeticError('native response has additional reducing subspaces')
    # A complex commutant basis is checked as a complex span separately;
    # its entries here are real, so the same real-coordinate solver applies.
    _coordinates((_clean(p), _clean(q)), commutant)
    return {
        'native_basis': G, 'declared_basis': D, 'native_hermitian': H,
        'native_source_coefficients': tuple(_clean(x) for x in native['source_coefficients']),
        'geometric_frame': T, 'native_chart': chart, 'native_chart_inverse': inverse,
        'canonical_native_basis': canonical_native, 'canonical_declared_basis': canonical_declared,
        'intersection_basis': intersection, 'union_basis': union,
        'native_from_canonical': native_from_canonical,
        'declared_from_canonical': declared_from_canonical,
        'canonical_from_native': canonical_from_native,
        'canonical_from_declared': canonical_from_declared,
        'commutant_basis': tuple(_clean(x) for x in (p, q)),
        'extra_hermitian_basis': _quadrupoles(), 'common_center': _clean(s.I*s.eye(6)),
        'relative_center': _clean(s.I*p),
        'source_constraint_rank': native['constraint_rank'],
        'source_coefficient_kernel_dimension': native['coefficient_kernel_dimension'],
    }


def exact_data():
    """Detached mapping; all matrix entries and nested basis tuples are immutable."""
    return dict(_exact_data())


def fixed_probes():
    """Eight coefficient-independent real inputs and complete Y measurements.

    The first five measure symmetric-traceless phase response on3prime.
    The last three compare matching coordinates across the two triplets.
    Each setting retains Y+, Y-, and the full four-dimensional complement.
    """
    eye = s.eye(6)
    settings = ((3, 4, (3, 4)), (4, 5, (4, 5)),
                (3, 4, (3,)), (3, 5, (3,)), (4, 5, (4,)),
                (0, 3, (0, 3)), (1, 4, (1, 4)), (2, 5, (2, 5)))
    result = []
    for a, b, support in settings:
        vector = sum((eye[:, i] for i in support), s.zeros(6, 1))/s.sqrt(len(support))
        plus, minus = (eye[:, a]+s.I*eye[:, b])/s.sqrt(2), (eye[:, a]-s.I*eye[:, b])/s.sqrt(2)
        effects = (plus*plus.conjugate().T, minus*minus.conjugate().T,
                   eye-_unit(a, a)-_unit(b, b))
        result.append({'input': _clean(vector), 'density': _clean(vector*vector.T),
                       'contrast': _clean(effects[0]-effects[1]),
                       'effects': tuple(_clean(x) for x in effects),
                       'coordinates': (a, b), 'support': support})
    return tuple(result)


def response_jacobian(generators):
    """Exact derivatives of complete fixed-record contrasts at zero parameter."""
    if type(generators) is not tuple or not generators:
        raise ValueError('a nonempty tuple of exact anti-Hermitian matrices is required')
    trusted = []
    concrete = (s.MutableDenseMatrix, s.ImmutableDenseMatrix, s.MutableSparseMatrix, s.ImmutableSparseMatrix)
    for original in generators:
        # Matrix subclasses can override iteration, shape or algebra. Reject
        # them before invoking any caller-controlled matrix method.
        if type(original) not in concrete:
            raise ValueError('a concrete trusted SymPy matrix representation is required')
        g = s.ImmutableMatrix(original)
        if (g.shape != (6, 6)
                or any(x.has(s.Float) or x.free_symbols or x.is_finite is not True for x in g)
                or not _zero(g+g.conjugate().T)):
            raise ValueError('finite exact numerical anti-Hermitian six-dimensional generators required')
        trusted.append(g)
    return _clean(s.Matrix([[s.trace(probe['contrast']*(g*probe['density']-probe['density']*g))
                            for g in trusted] for probe in fixed_probes()]))


def _entries(matrix):
    return [[str(x) for x in matrix.row(i)] for i in range(matrix.rows)]


@lru_cache(maxsize=1)
def _exact_certificate():
    data = _exact_data()
    G, D, H = data['native_basis'], data['declared_basis'], data['native_hermitian']
    native_visible = _rank(_columns(tuple(_clean(g-s.trace(g)*s.eye(6)/6) for g in G)))
    declared_visible = _rank(_columns(tuple(_clean(g-s.trace(g)*s.eye(6)/6) for g in D)))
    jacobian = response_jacobian(tuple(_clean(s.I*h) for h in H))
    declared_jacobian = response_jacobian(D)
    quadrupole_jacobian = response_jacobian(tuple(_clean(s.I*h) for h in data['extra_hermitian_basis']))[:5, :]
    central_jacobian = response_jacobian((data['relative_center'],))
    if (native_visible, declared_visible, _rank(jacobian[:5, :]), quadrupole_jacobian.det()) != (16, 12, 5, 24):
        raise ArithmeticError('fixed operational probes do not resolve the asserted response directions')
    if not _zero(declared_jacobian[:5, :]):
        raise ArithmeticError('the declared real-triplet grammar acquired an extra phase response')
    complete_native_jacobian = response_jacobian(G)
    if (not _zero(sum((complete_native_jacobian[i, :] for i in range(5, 8)), s.zeros(1, 17)))
            or sum(central_jacobian[5:, 0]) != -3):
        raise ArithmeticError('the operational relative-center separating functional failed')
    relative_phase = s.diag(s.I, s.I, s.I, 1, 1, 1)
    witness_probe = fixed_probes()[0]
    # A pi/4 relative phase gives maximal imaginary coherence and an exact
    # one-bit answer. It is an ideal group comparator, not an emitted word.
    z = (1+s.I)/s.sqrt(2)
    phase_witness = s.diag(1, 1, 1, z, s.conjugate(z), 1)
    probabilities = tuple(s.simplify(s.trace(effect*phase_witness*witness_probe['density']*
                                             phase_witness.conjugate().T))
                          for effect in witness_probe['effects'])
    if probabilities != (0, 1, 0):
        raise ArithmeticError('finite extra-phase discriminator has the wrong complete outcome law')
    return {
        'schema': 'oph.native_response_comparison.v1',
        'same_object': {'native_binary_load_four_dimension': 10626, 'full_color_dimension': 36,
            'geometric_dimension': 6, 'matching_colors': [0, 1],
            'coordinate_order': 'actual native geometric3, then geometric3prime',
            'source_port_pair_count': 78, 'source_constraint_rank': data['source_constraint_rank'],
            'source_coefficient_kernel_dimension': data['source_coefficient_kernel_dimension'],
            'full36_source_invariance_checked_before_compression': True,
            'native_to_orthonormal_chart': _entries(data['native_chart'])},
        'source_hermitian_basis': [_entries(h) for h in H],
        'source_coefficients': [_entries(c)[0] for c in data['native_source_coefficients']],
        'native_basis': [_entries(g) for g in G],
        'declared_basis': [_entries(g) for g in D],
        'canonical_native_basis': [_entries(g) for g in data['canonical_native_basis']],
        'canonical_declared_basis': [_entries(g) for g in data['canonical_declared_basis']],
        'span_certificates': {key: _entries(data[key]) for key in (
            'native_from_canonical', 'declared_from_canonical', 'canonical_from_native', 'canonical_from_declared')},
        'algebras': {'native_real_dimension': len(G), 'declared_real_dimension': len(D),
            'intersection_dimension': 11, 'union_Lie_closure_dimension': 18,
            'native': 'su(3)+su(3)+iR I6', 'declared': 'u(3)+so(3)',
            'intersection': 'su(3)+so(3)', 'union': 'u(3)+u(3)',
            'native_full_End6_observable_dimension': int(native_visible),
            'declared_full_End6_observable_dimension': int(declared_visible),
            'native_common_scalar_center_dimension': 1, 'native_relative_center_dimension': 0,
            'declared_relative_center_dimension': 1, 'declared_common_scalar_center_dimension': 0,
            'declared_algebra_is_a_subalgebra_of_native': False,
            'native_algebra_is_a_subalgebra_of_declared': False,
            'block_diagonal_tomography_native_dimension': 16,
            'block_diagonal_tomography_declared_dimension': 11},
        'invariant_restrictions': {'complex_commutant_dimension': 2,
            'commutant_basis': [_entries(x) for x in data['commutant_basis']],
            'reducing_subspaces': ['zero', 'geometric3', 'geometric3prime', 'all6'],
            'complex_dimensions': [0, 3, 3, 6], 'restricted_current_dimensions': [0, 9, 9, 17],
            'restricted_End_observable_dimensions': [0, 8, 8, 16],
            'faithful_declared12_obtained_by_invariant_state_subspace_restriction': False,
            'scope': 'fixed subspaces invariant under every native source response; no assertion about engineered controls or tensor encodings'},
        'fixed_records': {'settings': [{'input': _entries(row['input']),
            'effects': [_entries(effect) for effect in row['effects']]} for row in fixed_probes()],
            'outcomes_per_setting': 3, 'settings_count': 8,
            'native_source_Hermitian_jacobian': _entries(jacobian),
            'declared_jacobian': _entries(declared_jacobian),
            'five_phase_direction_jacobian': _entries(quadrupole_jacobian),
            'five_phase_generators_in_native_basis': _entries(_coordinates(
                G, tuple(_clean(s.I*h) for h in data['extra_hermitian_basis']))),
            'five_phase_direction_jacobian_determinant': '24',
            'extra_phase_rank_already_in_six_source_Hermitian_responses': 5,
            'declared_real_triplet_all_finite_Y_contrasts': '0',
            'finite_extra_phase_unitary': _entries(phase_witness),
            'finite_extra_phase_first_setting_probabilities': [str(x) for x in probabilities],
            'relative_center_jacobian': _entries(central_jacobian),
            'sum_last_three_native_contrast_derivatives': '0',
            'sum_last_three_relative_center_contrast_derivatives': '-3',
            'a_single_marker_is_a_complete_native_group_exclusion_test': False},
        'finite_group_invariant': {'native': 'det(U_three)/det(U_three_prime)=1',
            'common_phase_invariant': True,
            'uniform_declared_phase_ratio': 'exp(3 i t)',
            'uniform_declared_phase_native_iff': 't belongs to (2 pi/3) Z',
            'finite_central_witness': _entries(relative_phase),
            'finite_central_witness_determinant_ratio': str(relative_phase[:3, :3].det()/relative_phase[3:, 3:].det()),
            'identification_requires_complete_process_information': True},
        'scope': {'exact_current_and_fixed_record_comparison': True,
            'signed_native_port_Gram_schedules_source_executed_here': False,
            'fixed_pointer_source_words_executed_here': False,
            'native_C_E_repairs_select_either_current_grammar': False,
            'physical_gauge_group_or_particle_parameters_derived': False,
            'empirical_parameter_inputs': []}}


def exact_certificate():
    """JSON-safe owned certificate; returned records cannot poison the cache."""
    return deepcopy(_exact_certificate())


def source_bindings():
    paths = set(response.source_bindings())|set(declared.source_bindings())|{
        'model/common_source/native_response_comparison.py',
        'model/common_source/test_native_response_comparison.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
