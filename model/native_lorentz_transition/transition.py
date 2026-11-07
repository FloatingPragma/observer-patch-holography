"""Exact source-direction tests for conditional conformal Lorentz maps.

The native source supplies the twelve spatial directions. The real time
coordinate, null interpretation and proposed linear transition are inputs.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import permutations
import json

from ..common_source.null_tomography import (
    TENSOR_COMPONENTS, _max_row_sum, _null_tomography_decoder,
    decode_null_charges, null_tomography_certificate, null_tomography_design,
)
from ..common_source.slow_geometry import ARCHIVE, SPEC, _dot, raw_frame, source


def q(value=0):
    """Accept exact Q(sqrt(5)) data; floating approximations are rejected."""
    if type(value) in (int, Fraction):
        return (Fraction(value), Fraction(0))
    if (isinstance(value, (tuple, list)) and len(value) == 2
            and all(type(x) in (int, Fraction) for x in value)):
        return tuple(Fraction(x) for x in value)
    raise ValueError('an exact rational or rational Q(sqrt5) pair is required')


def matrix(rows):
    if (not isinstance(rows, (tuple, list)) or len(rows) != 4
            or any(not isinstance(row, (tuple, list)) or len(row) != 4
                   for row in rows)):
        raise ValueError('a four by four exact matrix is required')
    return tuple(tuple(q(value) for value in row) for row in rows)


def identity():
    return matrix([[int(i == j) for j in range(4)] for i in range(4)])


def metric():
    c = source()
    r2 = c.q5_add(c.GOLDEN, q(2))
    return matrix([[c.q5_sub(q(0), r2) if i == j == 0
                    else q(int(i == j)) for j in range(4)] for i in range(4)])


def multiply(left, right):
    left, right = matrix(left), matrix(right)
    return tuple(tuple(_dot(row, column) for column in zip(*right)) for row in left)


def scaled(value, rows):
    value, rows = q(value), matrix(rows)
    return tuple(tuple(source().q5_mul(value, x) for x in row) for row in rows)


def determinant(rows):
    rows, c = matrix(rows), source()
    total = q(0)
    for permutation in permutations(range(4)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(4) for j in range(i + 1, 4))
        term = q((-1) ** inversions)
        for i, j in enumerate(permutation):
            term = c.q5_mul(term, rows[i][j])
        total = c.q5_add(total, term)
    return total


def pullback(rows):
    rows = matrix(rows)
    return multiply(tuple(zip(*rows)), multiply(metric(), rows))


def null_charges(rows):
    """The twelve exact g(L k_p,L k_p), without fitting or projection."""
    return form_charges(pullback(rows))


def form_charges(tensor):
    """Evaluate a supplied symmetric bilinear form on the twelve source lifts."""
    tensor = matrix(tensor)
    if tensor != tuple(zip(*tensor)):
        raise ValueError('the form must be symmetric')
    components = tuple(tensor[i][j] for i, j in TENSOR_COMPONENTS)
    return tuple(_dot(row, components) for row in null_tomography_design())


def conformal_form(tensor):
    """Return lambda when the twelve zero tests imply tensor=lambda*G.

    This form-level function permits zero and negative lambda. Positivity
    follows only for an invertible congruence with the Lorentz metric.
    """
    tensor, c = matrix(tensor), source()
    charges = form_charges(tensor)
    if any(value != q(0) for value in charges):
        raise ValueError('the form fails a native null-direction test')
    if any(value != q(0) for value in decode_null_charges(charges)):
        raise AssertionError('zero charges did not decode to zero')
    g = metric()
    factor = c.q5_mul(tensor[0][0], c.q5_inv(g[0][0]))
    if tensor != scaled(factor, g):
        raise AssertionError('the native null kernel is not the metric span')
    return factor


def classify(rows):
    """Certify an invertible conformal Lorentz map from all twelve tests.

    Time reversal and spatial parity are classified separately; neither is
    silently accepted as a proper future-oriented transition.
    """
    rows, c = matrix(rows), source()
    det = determinant(rows)
    if det == q(0):
        raise ValueError('a chart transition must be invertible')
    # Reuse the existing rank-nine decoder rather than imposing the metric
    # identity as an additional independent acceptance requirement.
    factor = conformal_form(pullback(rows))
    if c.q5_sign(factor) <= 0:
        raise AssertionError('invertible Lorentz congruence must have positive scale')
    time_sign = c.q5_sign(rows[0][0])
    if time_sign == 0:
        raise AssertionError('a conformal Lorentz map cannot send time to zero time')
    return {'matrix': rows, 'metric_factor': factor, 'determinant': det,
            'future_preserving': time_sign > 0, 'proper': c.q5_sign(det) > 0}


def proper_future(rows):
    result = classify(rows)
    if not result['future_preserving']:
        raise ValueError('the transition reverses the chosen future cone')
    if not result['proper']:
        raise ValueError('the transition reverses spacetime orientation')
    return result


def compose(first, second):
    """First apply first, then second; verify multiplicative metric scale."""
    a, b = classify(first), classify(second)
    result = classify(multiply(b['matrix'], a['matrix']))
    if result['metric_factor'] != source().q5_mul(a['metric_factor'], b['metric_factor']):
        raise AssertionError('conformal factors did not compose multiplicatively')
    if result['future_preserving'] != (a['future_preserving'] == b['future_preserving']):
        raise AssertionError('time orientations did not compose')
    if result['proper'] != (a['proper'] == b['proper']):
        raise AssertionError('spacetime orientations did not compose')
    return result


def metric_residual(rows):
    """An exact entrywise bound from approximate null tests.

    This bounds the metric residual, not distance to a Lorentz matrix. No
    invertibility, positivity of the chosen factor or conditioning is inferred.
    """
    rows, c = matrix(rows), source()
    tensor, g = pullback(rows), metric()
    factor = c.q5_mul(tensor[0][0], c.q5_inv(g[0][0]))
    residual = tuple(tuple(c.q5_sub(tensor[i][j], c.q5_mul(factor, g[i][j]))
                           for j in range(4)) for i in range(4))
    charges = null_charges(rows)
    decoded = decode_null_charges(charges)
    if decoded != tuple(residual[i][j] for i, j in TENSOR_COMPONENTS):
        raise AssertionError('the decoder did not recover the gauge-fixed metric residual')
    def max_abs(values):
        maximum = q(0)
        for value in values:
            absolute = value if c.q5_sign(value) >= 0 else c.q5_sub(q(0), value)
            if c.q5_sign(c.q5_sub(absolute, maximum)) > 0:
                maximum = absolute
        return maximum
    actual = max_abs(x for row in residual for x in row)
    maximum_charge = max_abs(charges)
    constant = _max_row_sum(_null_tomography_decoder())
    bound = c.q5_mul(constant, maximum_charge)
    if c.q5_sign(c.q5_sub(bound, actual)) < 0:
        raise AssertionError('the exact decoder norm failed to bound the residual')
    return {'metric_factor': factor, 'residual': residual, 'maximum_charge': maximum_charge,
            'maximum_metric_residual': actual, 'decoder_norm': constant, 'bound': bound}


def boost(parameter=Fraction(1, 4)):
    """A rationally parametrized exact boost in the declared raw chart."""
    c, u = source(), q(parameter)
    r2 = c.q5_sub(q(0), metric()[0][0])
    ru2 = c.q5_mul(r2, c.q5_mul(u, u))
    denominator = c.q5_sub(q(1), ru2)
    if c.q5_sign(denominator) <= 0:
        raise ValueError('the future boost chart requires r2*u^2 < 1')
    inverse = c.q5_inv(denominator)
    diagonal = c.q5_mul(c.q5_add(q(1), ru2), inverse)
    off_diagonal = c.q5_mul(c.q5_mul(q(2), u), inverse)
    rows = [list(row) for row in identity()]
    rows[0][0] = rows[1][1] = diagonal
    rows[0][1] = off_diagonal
    rows[1][0] = c.q5_mul(r2, off_diagonal)
    return matrix(rows)


def controls():
    """Exact examples, not transition laws selected by the source scheduler."""
    rotation = matrix([[1, 0, 0, 0], [0, Fraction(3, 5), Fraction(-4, 5), 0],
                       [0, Fraction(4, 5), Fraction(3, 5), 0], [0, 0, 0, 1]])
    cyclic = matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0]])
    parity = matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    time_reversal = matrix([[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    null = (q(1),) + raw_frame()[0]
    collapse = matrix([[null[i] if j == 0 else q(0) for j in range(4)] for i in range(4)])
    return {'identity': identity(), 'boost': boost(), 'rotation': rotation,
            'native_cyclic_rotation': cyclic, 'dilation': scaled(Fraction(3, 2), identity()),
            'parity': parity, 'time_reversal': time_reversal, 'null_collapse': collapse}


def source_contract():
    """Replay the existing source/Gram/decoder identities used by the theorem."""
    # The underlying source loader is cached. Recheck its retained bytes even
    # if a caller populated that cache before this certificate was requested.
    pins = json.loads(SPEC.read_text())['source_sha256']
    for name in ('oph_exact/carrier.py', 'oph_fpe/dynamics/self_readback_repair_closure.py'):
        if hashlib.sha256((ARCHIVE / name).read_bytes()).hexdigest() != pins[name]:
            raise ValueError('the pinned native geometry source has changed')
    return null_tomography_certificate()
