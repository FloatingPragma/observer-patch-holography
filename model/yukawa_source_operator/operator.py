"""Exact duals for the three actual native Higgs invariant tensor lines.

An input is a Higgs derivative of an action kernel. This module never treats
an invariant tensor, a covariance, or a supplied mass as that derivative.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as s

from ..common_source import native_matter_attachment as matter


def clean(value):
    return s.ImmutableSparseMatrix(value).applyfunc(s.expand)


def zero(value):
    return all(s.simplify(x) == 0 for x in value)


@lru_cache(maxsize=1)
def tensor_lines():
    """Recompute the currents and solve their exact tensor equations."""
    representation = matter.representation()
    selected = list(representation['selected'])
    weak = representation['operators']['weak_C2'].extract(selected, selected)
    rows = []
    for row in matter.invariant_higgs_lines():
        support = {i for component in row['components']
                   for (i, _), value in component.todok().items() if value != 0}
        left = s.diag(*(int(i in support and weak[i, i] != 0) for i in range(15)))
        right = s.diag(*(int(i in support and weak[i, i] == 0) for i in range(15)))
        directed = tuple(clean(left * c * right) for c in row['components'])
        rank = s.trace(right)
        if any(not zero(d.H*d-right) for d in directed):
            raise ArithmeticError('native directed Higgs tensor normalization failed')
        if any(not zero(d+d.T-c) for d, c in zip(directed, row['components'], strict=True)):
            raise ArithmeticError('the oriented tensor lost its symmetric completion')
        rows.append({'origin': tuple(row['origin']), 'directed': directed,
                     'rank': int(rank), 'left': clean(left), 'right': clean(right),
                     'conjugate_higgs': row['conjugate']})
    if len(rows) != 3:
        raise ArithmeticError('exactly three native Higgs invariant lines required')
    return tuple(rows)


def assemble_derivatives(coefficients):
    """Two real-component Higgs directional derivatives, convention Y.

    For conventional Dirac Y the native_record_yukawa_field adapter inserts
    conjugate(Y). This API extracts the coefficient actually present in the
    all-left bilinear; complex conjugation is never inferred from eigenvalues.
    The directions H=v e_k use real v, so H and H* coincide. The retained
    conjugate_higgs flag distinguishes their imaginary-component derivatives.
    """
    lines = tensor_lines()
    if set(coefficients) != {row['origin'] for row in lines}:
        raise ValueError('all three exact charge origins are required')
    matrices = {origin: s.ImmutableMatrix(value) for origin, value in coefficients.items()}
    if any(value.shape != (3, 3) for value in matrices.values()):
        raise ValueError('three-by-three family coefficient matrices required')
    answer = []
    for component in range(2):
        derivative = s.zeros(45)
        for row in lines:
            directed = s.kronecker_product(row['directed'][component], matrices[row['origin']])
            derivative += directed + directed.T
        answer.append(clean(derivative))
    return tuple(answer)


def extract_coefficients(derivatives):
    """Contract an actual derivative with its exact Hilbert-Schmidt dual.

    With D the directed partial isometry, Y_ij is
    sum_ab conjugate(D_ab) V_(3a+i,3b+j) / rank(D).
    Exact reconstruction of both components rejects extra operators and
    inconsistent coefficients for the two components of one Higgs doublet.
    """
    if len(derivatives) != 2:
        raise ValueError('both Higgs component derivatives are required')
    derivatives = tuple(clean(value) for value in derivatives)
    if any(value.shape != (45, 45) for value in derivatives):
        raise ValueError('two complete 45-by-45 derivative matrices required')
    coefficients = {}
    for row in tensor_lines():
        d = row['directed'][1]
        coefficients[row['origin']] = clean(s.Matrix(3, 3, lambda i, j:
            sum(s.conjugate(weight)*derivatives[1][3*a+i, 3*b+j]
                for (a, b), weight in d.todok().items()) / row['rank']))
    replay = assemble_derivatives(coefficients)
    if any(not zero(actual-expected) for actual, expected in zip(derivatives, replay, strict=True)):
        raise ValueError('derivatives are outside the same three native Higgs tensor lines')
    return coefficients


def _inverse_root(metric):
    value = np.asarray(metric, complex)
    if (value.shape != (3, 3) or not np.isfinite(value).all()
            or not np.allclose(value, value.conj().T, atol=1e-13, rtol=0)):
        raise ValueError('finite Hermitian three-by-three kinetic metric required')
    eigenvalues, vectors = np.linalg.eigh(value)
    if min(eigenvalues) <= 0:
        raise ValueError('strictly positive kinetic metric required')
    return (vectors/np.sqrt(eigenvalues))@vectors.conj().T


def canonical_coefficient(coefficient, left_metric, right_metric, higgs_metric):
    """Whiten the all-left action; two-point metrics must be supplied.

    psi_L=E_L psi_L,c and psi_R=E_R psi_R,c give
    Y_c=E_L.T Y E_R/sqrt(Z_H). These are kinetic coefficients, rather than
    equal-time statistical covariance unless a transfer theorem proves it.
    """
    y = np.asarray(coefficient, complex)
    if y.shape != (3, 3) or not np.isfinite(y).all():
        raise ValueError('finite three-by-three action coefficient required')
    if (isinstance(higgs_metric, (bool, complex)) or not np.isscalar(higgs_metric)
            or not np.isfinite(higgs_metric) or higgs_metric <= 0):
        raise ValueError('finite strictly positive Higgs kinetic coefficient required')
    left, right = _inverse_root(left_metric), _inverse_root(right_metric)
    return left.T@y@right/np.sqrt(higgs_metric)


def composite_three_point():
    """Exact mixed CAR/composite overlap in the declared blank Fock state.

    This uses the physical Spin2 tensor M15 CAR attachment, which is a
    separate domain from the auxiliary five-mode representation. Occupation
    masks implement annihilation signs independently of the tensor formula.
    The state and pair operator are specified; no Higgs binding action is.
    """
    epsilon = s.Matrix([[0, 1], [-1, 0]])
    rows = []
    for line in tensor_lines():
        for component in range(2):
            directed = line['directed'][component]
            delta = clean(s.kronecker_product(epsilon, directed+directed.T))
            if delta.T != -delta:
                raise ArithmeticError('physical CAR pair kernel is not antisymmetric')
            amplitudes = {(1 << i) | (1 << j): delta[i, j]
                          for i in range(30) for j in range(i+1, 30) if delta[i, j] != 0}
            norm_squared = s.simplify(sum(s.conjugate(x)*x for x in amplitudes.values()))
            if norm_squared != 2*line['rank']:
                raise ArithmeticError('native scalar composite pair norm failed')
            for i in range(30):
                for j in range(i+1, 30):
                    mask = (1 << i) | (1 << j)
                    coefficient = amplitudes.get(mask, s.Integer(0))/s.sqrt(norm_squared)
                    # c_j c_i first annihilates i, then j, in sorted Fock order.
                    coefficient *= (-1)**((mask & ((1 << i)-1)).bit_count())
                    mask ^= 1 << i
                    coefficient *= (-1)**((mask & ((1 << j)-1)).bit_count())
                    mask ^= 1 << j
                    if mask != 0 or s.simplify(coefficient-delta[i, j]/s.sqrt(norm_squared)) != 0:
                        raise ArithmeticError('independent mixed CAR/composite contraction failed')
            rows.append({'charge_origin': list(map(str, line['origin'])),
                         'Higgs_tensor_component': component,
                         'pair_norm_squared_exact': str(norm_squared),
                         'normalized_nonzero_three_point_modulus_exact': str(1/s.sqrt(norm_squared)),
                         'nonzero_pair_amplitude_count': len(amplitudes),
                         'all_435_CAR_pair_contractions_checked': True})
    return {
        'domain': 'declared independent Spin2 tensor M15 CAR module and its blank Fock state',
        'physical_CAR_modes': 30,
        'native_matter_tensors_used': True,
        'composite_definition': 'B_dagger=(1/(2*sqrt(2*r))) sum_ij (epsilon tensor C)_ij c_i_dagger c_j_dagger',
        'mixed_three_point_exact': '<0|c_j c_i B_dagger|0>=(epsilon tensor C)_ij/sqrt(2*r)',
        'composite_vacuum_two_point_exact': '1',
        'fermion_vacuum_two_point_exact': 'I30',
        'rows': rows,
        'same_nonzero_overlap_exists_for_free_H_equals_omega_total_number': True,
        'independent_Higgs_field_in_free_action': False,
        'free_action_Higgs_derivative_exact': '0',
        'native_66_state_preparation_identified_with_this_CAR_blank': False,
        'native_source_record_execution_performed_here': False,
        'bound_Higgs_pole_or_interaction_coefficient_selected': False,
    }


def declared_action_coefficients():
    """Extract the derivative of the existing source-only trial Higgs action.

    The coefficient map a F, b F, J F and canonical metric are architectural
    choices of that producer. The exact source supplies a,b,J,F. This replay
    calculates their consequences without identifying them as physical Y.
    """
    action = matter.existing_action
    change, mapping = matter.action_intertwiner()
    change45 = s.kronecker_product(change, s.eye(3))
    derivatives = tuple(clean(change45.T*action.mass_kernel(1, component=k)*change45)
                        for k in range(2))
    for k in range(2):
        if (action.mass_kernel(0, component=k) != s.zeros(45)
                or not zero(action.mass_kernel(2, component=k)-2*action.mass_kernel(1, component=k))):
            raise ArithmeticError('the declared source Higgs dependence differs from its linear kernel')
    extracted = extract_coefficients(derivatives)
    expected = action.yukawa_matrices()
    rows = []
    for record in mapping:
        origin = tuple(map(s.Integer, record['charge_origin']))
        name = record['existing_action_channel']
        scale = s.sympify(record['same_current_tensor_scale_exact'])
        coefficient = extracted[origin]
        if not zero(coefficient-scale*expected[name]):
            raise ArithmeticError('source Higgs variation does not recover its declared coefficient')
        if not coefficient.is_diagonal():
            raise ArithmeticError('the source-only trial coefficient lost its diagonal form')
        singular = sorted((s.simplify(abs(x)) for x in coefficient.diagonal()), key=float)
        rows.append({'charge_origin': list(map(str, origin)), 'existing_action_channel': name,
                     'native_to_existing_tensor_scale_exact': str(scale),
                     'coefficient_matrix_exact': [[str(x) for x in row] for row in coefficient.tolist()],
                     'unit_Higgs_coordinate_singular_values_exact': list(map(str, singular)),
                     'unit_Higgs_coordinate_singular_values_diagnostic': list(map(float, singular)),
                     'largest_over_smallest_exact': str(s.simplify(singular[-1]/singular[0]))})
    return {
        'source_producer': 'model/common_source/chiral_candidate/action.py',
        'Higgs_derivative': 'real directional derivative H=v*e_k: P45.T K_existing(v=1,component=k) P45; K_existing(v)=v K_existing(1)',
        'all_three_native_duals_recover_existing_derivative_exactly': True,
        'source_log_odds_exact': {name: str(value) for name, value in action.source_coefficients().items()},
        'declared_coefficient_map': {'up': 'a*F', 'down': 'b*F', 'charged_lepton': 'J*F'},
        'canonical_fermion_kinetic_metrics_I_are_declared': True,
        'unit_Higgs_coordinate_is_declared': True,
        'dynamical_Higgs_kinetic_residue_supplied': False,
        'fully_canonical_Yukawa_requires_division_by_sqrt_Z_H': True,
        'rows': rows,
        'physical_Higgs_and_charged_particle_attachment_proved': False,
        'coefficient_map_forced_by_native_repair': False,
        'measured_mass_or_Yukawa_inputs': [],
    }


def certificate():
    """Exact recovery of all 27 independent complex coefficient coordinates."""
    lines = tensor_lines()
    symbolic = {row['origin']: s.Matrix(3, 3, s.symbols(f'y{k}_0:9'))
                for k, row in enumerate(lines)}
    derivatives = assemble_derivatives(symbolic)
    extracted = extract_coefficients(derivatives)
    if extracted != symbolic:
        raise ArithmeticError('native dual failed arbitrary complex family recovery')
    controls = []
    for scale in (1, 2):
        ys = {row['origin']: scale*s.eye(3) for row in lines}
        actual = extract_coefficients(assemble_derivatives(ys))
        if actual != ys:
            raise ArithmeticError('same-current coefficient control failed')
        controls.append({'coefficient': str(scale), 'family_singular_values': [str(scale)]*3})
    return {
        'native_current_count': len(matter.representation()['matter_currents']),
        'native_invariant_line_count': len(lines),
        'complex_coefficient_dimension': 27,
        'tensor_lines': [{'charge_origin': list(map(str, row['origin'])),
                          'directed_rank': row['rank'],
                          'unit_nonzero_singular_values_exact': True,
                          'dual_denominator_exact': str(row['rank']),
                          'uses_conjugate_Higgs': bool(row['conjugate_higgs'])} for row in lines],
        'all_27_symbolic_coordinates_recovered_exactly': True,
        'both_Higgs_components_reconstructed_exactly': True,
        'arbitrary_nonsymmetric_complex_coefficients_supported': True,
        'physical_source_selected_Higgs_derivative_supplied_by_this_certificate': False,
        'mixed_fermion_composite_three_point': composite_three_point(),
        'existing_source_only_declared_action': declared_action_coefficients(),
        'same_current_same_canonical_metric_same_unit_Higgs_controls': controls,
        'control_scope': 'representation and kinetic premises fixed; the Higgs derivative and interacting action differ',
        'canonical_normalization': 'Y_c=Z_L^(-T/2) Y Z_R^(-1/2)/sqrt(Z_H)',
        'missing_physical_source_datum': 'the projected Higgs derivative of the same source effective inverse fermion propagator, with fermion and Higgs kinetic residues in its declared clock',
    }
