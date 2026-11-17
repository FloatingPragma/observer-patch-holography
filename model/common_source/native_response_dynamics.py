"""Finite feedback from complete reference-test disturbance, without a menu cutoff.

The controller acts on eighteen real Hamiltonian coefficients, not on unknown
quantum states. Its declared update is x -> (I-eta G_p)x, where G_p is the
failure-probability quadratic form of the collective action on the biased Bell
reference. Neither a desired rank nor a twelve-generator list enters that
update. At fixed interior weight 0<p<1 its limiting observable range has dimension six.
The controlled joint limit p -> 0, m -> infinity, mp -> 0 instead retains
twelve observable currents. Reference cooling or unchanged C/E repair alone
does not implement this additional classical feedback law.

Convergence is uniform only on a fixed coefficient and actuation-time budget.
At every finite interior weight and update count the filter is invertible:
unbounded input amplification can restore every suppressed current. The
class-sum phase must also be available as a work actuator; a bias applied only
to the reference cannot create the missing center of the original source17.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import hashlib
from pathlib import Path

import sympy as s


def _exact(value):
    if type(value) is int:
        result = s.Integer(value)
    elif type(value) is Fraction:
        result = s.Rational(value.numerator, value.denominator)
    elif isinstance(value, s.Expr):
        result = value
    else:
        raise ValueError('an original finite exact real coefficient is required')
    if (result.has(s.Float) or result.free_symbols or result.is_real is not True
            or result.is_finite is not True):
        raise ValueError('an original finite exact real coefficient is required')
    return result


def _weight(p):
    value = _exact(p)
    if not bool(value >= 0) or not bool(value <= 1):
        raise ValueError('reference weight must lie in [0,1]')
    return value


def _step(value):
    value = _exact(value)
    if not bool(value > 0) or not bool(value <= s.Rational(3, 4)):
        raise ValueError('feedback step must lie in (0,3/4]')
    return value


def _steps(value):
    if type(value) is not int or value < 0:
        raise ValueError('feedback count must be an original nonnegative integer')
    return value


def _clean(value):
    return s.ImmutableMatrix(value.applyfunc(lambda x: s.simplify(s.expand(x))))


def _unit(i, j, dimension=3):
    value = s.zeros(dimension)
    value[i, j] = 1
    return value


def full_basis():
    """Identical complete u(3) bases on both geometric triplets, HS orthonormal.

    This frame contains all eighteen controls. Its nine entries per block are
    three real rotations, five imaginary quadrupoles and one scalar phase.
    It is a coordinate choice for the full source, not a retained-current menu.
    """
    rotations = tuple((_unit(i, j)-_unit(j, i))/s.sqrt(2)
                      for i in range(3) for j in range(i+1, 3))
    quadrupoles = tuple(s.I*(_unit(i, j)+_unit(j, i))/s.sqrt(2)
                        for i in range(3) for j in range(i+1, 3))
    quadrupoles += (s.I*s.diag(1, -1, 0)/s.sqrt(2),
                    s.I*s.diag(1, 1, -2)/s.sqrt(6))
    block = (*rotations, *quadrupoles, s.I*s.eye(3)/s.sqrt(3))
    return tuple(s.ImmutableMatrix(s.diag(x, s.zeros(3))) for x in block)+tuple(
        s.ImmutableMatrix(s.diag(s.zeros(3), x)) for x in block)


def _bell(offset):
    vector = s.zeros(36, 1)
    for j in range(offset, offset+3):
        vector[6*j+j] = 1/s.sqrt(3)
    return vector


@lru_cache(maxsize=1)
def _gram_polynomial():
    """Build the cost from all collective actions, before any spectral labels."""
    basis = full_basis()
    gram = _clean(s.Matrix([[s.re(s.trace(x.conjugate().T*y)) for y in basis] for x in basis]))
    if gram != s.eye(18):
        raise ArithmeticError('the complete native coefficient metric is not Hilbert-Schmidt orthonormal')
    images = []
    means = []
    for ray in (_bell(0), _bell(3)):
        moved = s.Matrix.hstack(*((s.kronecker_product(x, s.eye(6))
                                  +s.kronecker_product(s.eye(6), x))*ray for x in basis))
        images.append(moved)
        means.append(ray.conjugate().T*moved)
    if images[0].conjugate().T*images[1] != s.zeros(18):
        raise ArithmeticError('the complete triplet reference outputs are not orthogonal')
    p = s.Symbol('p', real=True)
    mean = p*means[0]+(1-p)*means[1]
    raw = p*images[0].conjugate().T*images[0]+(1-p)*images[1].conjugate().T*images[1]
    cost = _clean((raw-mean.conjugate().T*mean).applyfunc(s.re))
    return p, cost


def cost_gram(p):
    """Exact coefficient cost: 1-Pr(reference survives at time t)=t² x'Gx+O(t⁴).

    Both pointer outcomes are retained and the reference is reset between
    tests. This unweighted probability cost is not energy dissipation and
    is not a trace-decreasing dynamics on the retained work object.
    """
    p = _weight(p)
    symbol, polynomial = _gram_polynomial()
    return s.ImmutableMatrix(polynomial.subs(symbol, p))


@lru_cache(maxsize=1)
def _spectral_data():
    # Diagonalize a generic exact interior point. The eigenspaces are then
    # certified against the complete polynomial, not assumed from a target.
    point = cost_gram(Fraction(1, 4))
    projectors = {}
    for value, multiplicity, vectors in point.eigenvects():
        matrix = s.Matrix.hstack(*vectors)
        if matrix.cols != multiplicity:
            raise ArithmeticError('the complete reference cost eigenspaces differ')
        projectors[value] = _clean(matrix*(matrix.T*matrix).inv()*matrix.T)
    expected = {s.Integer(0): 7, s.Rational(1, 3): 5, s.Integer(1): 5, s.Rational(1, 2): 1}
    if set(projectors) != set(expected) or any(s.trace(projectors[v]) != n for v, n in expected.items()):
        raise ArithmeticError('the complete reference disturbance spectrum changed')
    zero, first, second, relative = (projectors[v] for v in (0, s.Rational(1, 3), 1, s.Rational(1, 2)))
    symbol, polynomial = _gram_polynomial()
    predicted = 4*symbol*first/3+4*(1-symbol)*second/3+8*symbol*(1-symbol)*relative/3
    if (_clean(polynomial-predicted) != s.zeros(18)
            or sum((zero, first, second, relative), s.zeros(18)) != s.eye(18)):
        raise ArithmeticError('the spectral decomposition fails on the full reference-weight polynomial')
    return zero, first, second, relative


def limit_projection():
    """Validation target derived as ker G_0; it does not enter feedback updates."""
    return s.ImmutableMatrix(s.eye(18)-s.Rational(3, 4)*cost_gram(0))


def feedback_matrix(p, steps, step=Fraction(1, 2)):
    """Exact m updates of I-eta G; all eighteen input coefficients remain accepted."""
    p, steps, step = _weight(p), _steps(steps), _step(step)
    zero, first, second, relative = _spectral_data()
    return s.ImmutableMatrix(zero+(1-step*4*p/3)**steps*first
        +(1-step*4*(1-p)/3)**steps*second
        +(1-step*8*p*(1-p)/3)**steps*relative)


def feedback_bound(p, steps, step=Fraction(1, 2)):
    """Exact coefficient-operator distance to the biased-limit Lie projection."""
    p, steps, step = _weight(p), _steps(steps), _step(step)
    if bool(p > s.Rational(1, 2)):
        raise ValueError('the second-triplet limit bound requires p<=1/2')
    low = 1-(1-step*8*p*(1-p)/3)**steps
    high = (1-step*4*(1-p)/3)**steps
    return {'low_current_loss': low, 'high_current_gain': high,
            'coefficient_operator_error': s.Max(low, high),
            'low_current_loss_upper_bound': 8*step*steps*p/3,
            'high_current_gain_upper_bound': (1-2*step/3)**steps}


def continuous_bound(p, time):
    """Exact bound for the declared coefficient flow dx/dtau=-G_p x."""
    p, time = _weight(p), _exact(time)
    if bool(p > s.Rational(1, 2)) or bool(time < 0):
        raise ValueError('the second-triplet limit requires p<=1/2 and nonnegative time')
    low, high = 1-s.exp(-8*time*p*(1-p)/3), s.exp(-4*time*(1-p)/3)
    return {'low_current_loss': low, 'high_current_gain': high,
            'coefficient_operator_error': s.Max(low, high),
            'low_current_loss_upper_bound': 8*time*p/3,
            'high_current_gain_upper_bound': s.exp(-2*time/3)}


def limiting_group_failure_bound(p):
    """Sharp survival-failure bound for the entire limiting current group.

    The collective overlap is p*a+(1-p)*z, where |a|<=1 is the first
    triplet Bell overlap and |z|=1 is the second triplet's phase. For
    p<=1/2 its modulus is at least 1-2p, attained by a=-z. For p>=1/2,
    take U_A U_A^T=z*diag(-1,exp(i theta),exp(-i theta)), with
    cos(theta)=(1-3(1-p)/p)/2 in [-1,1/2]. This gives
    a=-(1-p)*z/p and cancels the overlap.
    The bound therefore cannot remain 4p(1-p) on the upper half interval.
    """
    p = _weight(p)
    return s.simplify(4*p*(1-p)) if bool(p <= s.Rational(1, 2)) else s.Integer(1)


def _coordinates(matrix):
    return _clean(s.Matrix([s.re(s.trace(x.conjugate().T*matrix)) for x in full_basis()]))


def canonical_certificate():
    """Exact spectrum, limiting Lie closure, and non-quotient control."""
    zero, first, second, relative = _spectral_data()
    projection = limit_projection()
    if projection != zero+first+relative or projection*projection != projection:
        raise ArithmeticError('the biased limiting cost kernel differs')
    basis = full_basis()
    columns = projection.columnspace()
    retained = tuple(_clean(sum((a*x for a, x in zip(column, basis, strict=True)), s.zeros(6)))
                     for column in columns)
    if len(retained) != 13:
        raise ArithmeticError('the limiting feedback range has the wrong dimension')
    for i, x in enumerate(retained):
        for y in retained[:i]:
            coordinates = _coordinates(x*y-y*x)
            if projection*coordinates != coordinates:
                raise ArithmeticError('the limiting feedback range is not closed under current brackets')
    excluded = tuple(_clean(sum((a*x for a, x in zip(column, basis, strict=True)), s.zeros(6)))
                     for column in second.columnspace())
    for i, x in enumerate(excluded):
        for y in excluded[:i]:
            coordinates = _coordinates(x*y-y*x)
            if projection*coordinates != coordinates:
                raise ArithmeticError('the high-high current bracket leaves the low range')
        for y in retained:
            coordinates = _coordinates(x*y-y*x)
            if second*coordinates != coordinates:
                raise ArithmeticError('the low-high current bracket leaves the high range')
    # Erased right quadrupoles have a retained rotation as their bracket.
    a, b = basis[12], basis[15]
    bracket = _coordinates(a*b-b*a)
    if (projection*_coordinates(a) != s.zeros(18, 1)
            or projection*_coordinates(b) != s.zeros(18, 1)
            or bracket == s.zeros(18, 1) or projection*bracket != bracket):
        raise ArithmeticError('the non-homomorphic-projection control failed')
    common = _coordinates(s.I*s.eye(6))
    if projection*common != common or zero*common != common:
        raise ArithmeticError('the common scalar phase is not retained')
    native17 = s.Matrix.hstack(*(s.eye(18)[:, j] for j in range(18) if j not in (8, 17)), common)
    if native17.rank() != 17 or (projection*native17).rank() != 12:
        raise ArithmeticError('feedback supplied a relative center absent from the original source')
    return {'schema': 'oph.native_response_dynamics.v1',
        'cost_definition': 'G_p is the unweighted collective-reference failure-probability quadratic coefficient (half the Hessian) in the complete HS-orthonormal native18 frame',
        'spectrum': [['0', 7], ['4p/3', 5], ['4(1-p)/3', 5], ['8p(1-p)/3', 1]],
        'exact_polynomial_spectral_identity_checked': True,
        'finite_interior_kernel_raw_dimension': 7,
        'finite_interior_kernel_observable_dimension': 6,
        'feedback_law': 'x_(m+1)=(I-eta G_p)x_m; 0<eta<=3/4',
        'finite_interior_coefficients_removed_exactly': False,
        'finite_interior_filter_is_invertible': True,
        'finite_reachable_algebra_with_unbounded_input_coefficients': 'full18; no algebra reduction',
        'uniform_limit_requires_fixed_input_and_time_budgets': True,
        'original17_controls_limit_without_work_class_sum_phase': {'raw_dimension': 12, 'observable_dimension': 11},
        'limiting_projection_definition': 'P=I-(3/4)G_0; validation only, never used in the feedback recurrence',
        'joint_limit': 'p->0, m->infinity, mp->0 at fixed positive eta',
        'Euler_exact_operator_error': 'max(1-(1-8 eta p(1-p)/3)^m, (1-4 eta(1-p)/3)^m)',
        'Euler_error_upper_bound': 'max((8 eta/3)mp,(1-2 eta/3)^m), for0<=p<=1/2',
        'continuous_exact_operator_error': 'max(1-exp(-8 tau p(1-p)/3),exp(-4 tau(1-p)/3))',
        'continuous_joint_limit': 'p->0, tau->infinity, tau*p->0',
        'limit_raw_dimension': 13, 'limit_observable_dimension': 12,
        'limit_range_bracket_closed': True,
        'low_high_and_high_high_bracket_inclusions_checked': True,
        'projection_is_a_Lie_homomorphism': False,
        'fixed_p_long_time_observable_dimension': 6,
        'fixed_p_long_time_domain': '0<p<1; endpoint kernels instead have13 raw and12 observable directions',
        'finite_bracket_leakage_bound': '||(I-P)[F_m X,F_m Y]||_HS <=2 high_current_gain ||X||_HS ||Y||_HS',
        'complete_channel_bound': 'diamond distance <=min(2,2 epsilon B), B=sum_j sup_h |t_j(h)| ||x_j(h)||_2, for matched adaptive schedules with uniform per-stage bounds and all records retained',
        'arbitrary_untouched_reference_domain': True,
        'low_group_reference_failure_bound': '4p(1-p) for0<=p<=1/2; 1 for1/2<=p<=1, sharply and uniformly for every finite product in the limiting current group',
        'relative_center_reference_failure': '4p(1-p) sin(t)^2 for X=i P_first',
        'finite_cost_is_energy_loss': False,
        'unchanged_native_C_E_repairs_implement_feedback': False,
        'ordinary_reference_cooling_alone_implements_feedback': False,
        'additional_architecture': ['reference test and reset with both outcomes retained',
            'complete eighteen-coefficient classical controller and its stated metric',
            'promoted class-sum phase available as a work actuator, not only as a reference bias',
            'declared gradient feedback update and current actuator',
            'bias and controller count joint scaling with bounded input coefficients and actuation time'],
        'physical_gauge_redundancy_or_parameters_derived': False}


def source_certificate(bias=Fraction(16), steps=16, step=Fraction(1, 2)):
    """Join the complete native phase augmentation and its biased ground state."""
    from . import native_response_biased as biased, native_response_comparison as comparison
    from . import native_persistent_response as response
    data = comparison.exact_data()
    actual_phase = s.Matrix(response.rotation_endogeneity()['near_fivefold_class_sum_exact'])
    phase = comparison._clean(s.I*data['native_chart_inverse']*actual_phase*data['native_chart'])
    augmented = (*data['native_basis'], phase)
    basis = full_basis()
    comparison._coordinates(augmented, basis)
    comparison._coordinates(basis, augmented)
    # Identify the limit only after computing the full-source feedback law.
    # The declared current basis is a comparison target, never a filter input.
    retained = tuple(_clean(sum((a*x for a, x in zip(column, basis, strict=True)), s.zeros(6)))
                     for column in limit_projection().columnspace())
    sections = tuple(_clean(x-s.trace(x[3:, 3:])*s.eye(6)/3) for x in retained)
    pivots = comparison._columns(sections).to_DM(extension=True).rref()[1]
    section = tuple(sections[j] for j in pivots)
    comparison._coordinates(data['declared_basis'], section)
    comparison._coordinates(section, data['declared_basis'])
    ground = biased.exact_data(bias)
    p = ground['first_triplet_weight']
    # Independently recover the unweighted Gram from the actual full36 ground.
    ray = s.Matrix(ground['ground_geometric'])
    if ray.shape != (36, 1) or s.simplify((ray.conjugate().T*ray)[0]) != 1:
        raise ArithmeticError('the actual biased reference is not a normalized full36 ray')
    moved = s.Matrix.hstack(*((s.kronecker_product(x, s.eye(6))
                              +s.kronecker_product(s.eye(6), x))*ray for x in basis))
    direct = _clean((moved.conjugate().T*(s.eye(36)-ray*ray.conjugate().T)*moved).applyfunc(s.re))
    if _clean(direct-cost_gram(p)) != s.zeros(18):
        raise ArithmeticError('the actual biased reference has a different complete disturbance Gram')
    bounds = feedback_bound(p, steps, step)
    return {'canonical': canonical_certificate(), 'bias': str(bias), 'steps': _steps(steps),
        'step': str(_step(step)), 'source_weight_p_exact': str(p),
        'complete_native18_bidirectional_span_checked': True,
        'limit_common_phase_quotient_matches_declared12_as_validation_only': True,
        'actual_full36_ground_disturbance_Gram_checked': True,
        'finite_bounds_exact': {key: str(value) for key, value in bounds.items()}}


def source_bindings():
    from . import native_response_biased as biased, native_response_comparison as comparison
    paths = set(biased.source_bindings())|set(comparison.source_bindings())|{
        'model/common_source/native_response_dynamics.py',
        'model/common_source/test_native_response_dynamics.py'}
    root = Path(__file__).resolve().parents[2]
    return {path: hashlib.sha256((root/path).read_bytes()).hexdigest() for path in sorted(paths)}
