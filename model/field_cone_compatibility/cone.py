"""Exact clock, principal-cone and frame compatibility over Q(sqrt(5)).

The source supplies the twelve directions and their tensor decoder. Spatial
coefficients, a coordinate Jacobian and clock pullback remain supplied data.
No evolution, radial action or source chart is changed by these checks.
"""
from fractions import Fraction as F
from itertools import permutations

from ..native_lorentz_transition import transition as t

SOURCE_FILES = (
    'model/field_cone_compatibility/__init__.py',
    'model/field_cone_compatibility/__main__.py',
    'model/field_cone_compatibility/cone.py',
    'model/field_cone_compatibility/test_cone.py',
    'model/field_cone_compatibility/README.md',
    'model/field_cone_compatibility/PROOF.md',
    'model/native_lorentz_transition/transition.py',
    'model/common_source/null_tomography.py',
    'model/common_source/slow_geometry.py',
)


def m3(rows):
    if (not isinstance(rows, (tuple, list)) or len(rows) != 3
            or any(not isinstance(row, (tuple, list)) or len(row) != 3
                   for row in rows)):
        raise ValueError('a three by three exact matrix is required')
    return tuple(tuple(t.q(x) for x in row) for row in rows)


def _transpose(rows):
    return tuple(zip(*rows))


def _multiply(a, b):
    return tuple(tuple(t._dot(row, column) for column in zip(*b)) for row in a)


def _det(rows):
    c = t.source()
    total = t.q(0)
    for order in permutations(range(len(rows))):
        inversions = sum(order[i] > order[j] for i in range(len(rows))
                         for j in range(i+1, len(rows)))
        term = t.q((-1)**inversions)
        for i, j in enumerate(order):
            term = c.q5_mul(term, rows[i][j])
        total = c.q5_add(total, term)
    return total


def _positive(value, name):
    value = t.q(value)
    if t.source().q5_sign(value) <= 0:
        raise ValueError(name + ' must be strictly positive')
    return value


def inverse(rows):
    """Exact three-dimensional inverse, including off-diagonal pivots."""
    a, c = m3(rows), t.source()
    work = [list(row)+[t.q(int(i == j)) for j in range(3)]
            for i, row in enumerate(a)]
    for j in range(3):
        pivot = next((i for i in range(j, 3) if work[i][j] != t.q(0)), None)
        if pivot is None:
            raise ValueError('the matrix must be invertible')
        work[j], work[pivot] = work[pivot], work[j]
        inv = c.q5_inv(work[j][j])
        work[j] = [c.q5_mul(inv, x) for x in work[j]]
        for i in range(3):
            if i != j:
                coefficient = work[i][j]
                work[i] = [c.q5_sub(x, c.q5_mul(coefficient, y))
                           for x, y in zip(work[i], work[j], strict=True)]
    return tuple(tuple(row[3:]) for row in work)


def positive_tensor(rows):
    """Sylvester's criterion uses exact ordered-real Q(sqrt(5)) minors."""
    a = m3(rows)
    if a != _transpose(a):
        raise ValueError('the spatial tensor must be symmetric')
    for n in (1, 2, 3):
        if t.source().q5_sign(_det(tuple(row[:n] for row in a[:n]))) <= 0:
            raise ValueError('the spatial tensor must be positive definite')
    return a


def field_form(spatial_tensor, jacobian, clock_squared):
    """Pull back -dt^2+dx^T a^-1 dx by dx=J dz, dt=sqrt(s) d tau.

Only the square of the positive clock conversion is needed. The returned
form acts on tangent displacements, not on frequency covectors.
"""
    a, j = positive_tensor(spatial_tensor), m3(jacobian)
    s = _positive(clock_squared, 'clock squared')
    if _det(j) == t.q(0):
        raise ValueError('the spatial Jacobian must be invertible')
    h = _multiply(_transpose(j), _multiply(inverse(a), j))
    return t.matrix([[t.source().q5_sub(t.q(0), s), 0, 0, 0],
                     *[[0, *row] for row in h]])


def match(spatial_tensor, jacobian, clock_squared):
    """Accept precisely when all twelve supplied native null rays are null.

The imported exact rank-nine decoder makes this a sufficient test for
equality of the entire principal cone. It does not identify source charts.
"""
    form = field_form(spatial_tensor, jacobian, clock_squared)
    factor = t.conformal_form(form)
    r2 = t.source().q5_sub(t.q(0), t.metric()[0][0])
    expected = t.source().q5_mul(t.q(clock_squared), t.source().q5_inv(r2))
    if factor != expected or t.source().q5_sign(factor) <= 0:
        raise AssertionError('the positive field clock did not fix the cone factor')
    return {'form': form, 'factor': factor, 'null_charges': t.form_charges(form)}


def congruence(form, change):
    form, change = t.matrix(form), t.matrix(change)
    if form != _transpose(form):
        raise ValueError('the form must be symmetric')
    if t.determinant(change) == t.q(0):
        raise ValueError('a frame change must be invertible')
    return t.multiply(_transpose(change), t.multiply(form, change))


def frame_covariance(spatial_tensor, jacobian, clock_squared, change):
    """Change both charts by the same arbitrary invertible linear map."""
    matched = match(spatial_tensor, jacobian, clock_squared)
    field = congruence(matched['form'], change)
    observer = congruence(t.metric(), change)
    if field != t.scaled(matched['factor'], observer):
        raise AssertionError('a common change of frame broke cone compatibility')
    return {'field_form': field, 'observer_form': observer,
            'factor': matched['factor']}


def radial_product_form(radius, angular_coefficient, clock_squared, direction=(0, 0, 1)):
    """Round product control under the FIXED map x=r*n and one time scaling.

The spatial metric in x coordinates is nn^T+(kappa*r^2)^-1(I-nn^T).
This is a control for the product geometry, not a claimed formula for the
source's nonconstant angular operator or for every coordinate dictionary.
"""
    r = _positive(radius, 'radius')
    k = _positive(angular_coefficient, 'angular coefficient')
    s = _positive(clock_squared, 'clock squared')
    if not isinstance(direction, (tuple, list)) or len(direction) != 3:
        raise ValueError('a unit three-direction is required')
    n, c = tuple(t.q(x) for x in direction), t.source()
    if t._dot(n, n) != t.q(1):
        raise ValueError('the radial direction must be a unit vector')
    tangential = c.q5_inv(c.q5_mul(k, c.q5_mul(r, r)))
    h = tuple(tuple(c.q5_add(c.q5_mul(n[i], n[j]),
                       c.q5_mul(tangential, c.q5_sub(t.q(int(i == j)),
                                                    c.q5_mul(n[i], n[j]))))
                    for j in range(3)) for i in range(3))
    return t.matrix([[c.q5_sub(t.q(0), s), 0, 0, 0], *[[0, *row] for row in h]])


def certificate():
    # This replays the native Gram, rank-nine and exact decoder identities.
    native = t.source_contract()
    c, r2 = t.source(), t.source().q5_sub(t.q(0), t.metric()[0][0])
    j = m3([[1, F(1, 2), 0], [0, 2, F(1, 3)], [0, 0, 3]])
    a = _multiply(j, _transpose(j))
    positive = match(a, j, r2)
    if positive['form'] != t.metric() or positive['factor'] != t.q(1):
        raise AssertionError('the constructive non-diagonal match failed')
    boost = t.boost()
    transformed = frame_covariance(a, j, r2, boost)
    if transformed['observer_form'] != t.metric():
        raise AssertionError('the native boost did not preserve the observer form')
    dilation = t.scaled(2, boost)
    dilated = congruence(positive['form'], dilation)
    if t.conformal_form(dilated) != t.q(4):
        raise AssertionError('frame conformal factors did not multiply')
    first = radial_product_form(1, 1, r2)
    second = radial_product_form(2, 1, r2)
    if t.conformal_form(first) != t.q(1):
        raise AssertionError('the matched radius failed')
    charges = t.form_charges(second)
    if all(x == t.q(0) for x in charges):
        raise AssertionError('the second-radius mismatch disappeared')
    encode = lambda value: c.q5_str(value)
    em = lambda matrix: [[encode(x) for x in row] for row in matrix]
    return {'schema': 'oph.field_cone_compatibility.v1',
            'native_source': native['registered_binding'],
            'native_null_design_rank': native['exact_design_rank'],
            'constructive_match': {'spatial_tensor': em(a), 'jacobian': em(j),
                'clock_squared': encode(r2), 'form': em(positive['form']),
                'factor': encode(positive['factor']), 'all_twelve_charges_zero': True},
            'frame_covariance': {'boost_factor': '1', 'dilated_boost_factor': '4',
                'same_change_must_act_on_both_forms': True},
            'radial_control': {'radii': ['1', '2'], 'angular_coefficient': '1',
                'clock_squared': encode(r2), 'first_form': em(first),
                'second_form': em(second), 'second_null_charges': [encode(x) for x in charges],
                'fixed_dictionary': 'x=r*n; t=sqrt(clock_squared)*tau',
                'one_clock_cannot_match_both_radii': True},
            'scope': {'principal_cone_criterion': True,
                'coordinate_clock_and_spatial_tensor_are_supplied': True,
                'source_selected_common_causal_field_clock': False,
                'full_field_law_or_interacting_continuum_derived': False,
                'curved_Lorentzian_or_other_coordinate_descriptions_excluded': False,
                'source_radial_action_changed': False}}
