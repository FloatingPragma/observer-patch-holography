"""Exact identification of the native definite-current instrument.

Nine scalar quadratures identify the nine real coordinates of the existing
record-law family. This identifies its common-current CP branches, not an
arbitrary fine record basis or a coherent dilation of those branches.

All public exact scalars are detached. No caller-owned mutable Fraction or
family.S object is cached or retained in a probe descriptor.
"""
from fractions import Fraction as F
from math import gcd
from typing import NamedTuple

from ..record_laws import family


PARAMETERS = ('k', 'g', 'h', 'm', 'n')
COORDINATES = ('k_re', 'k_im', 'g_re', 'g_im', 'h_re', 'h_im',
               'm_re', 'n_re', 'n_im')
SOURCE_PATHS = ('model/repair_identification/core.py',
                'model/record_laws/family.py', 'model/observer/native.py',
                'model/predictions/observer_instrument.py')


class Probe(NamedTuple):
    name: str
    parameter: str
    phase: int
    current: int
    keys: tuple
    input_pair: tuple
    output_pair: tuple
    baseline: tuple


def _probes():
    specifications = (
        ('k', ((0, 0), (1, 1)), 0),
        ('g', ((0, 0), (0, 1)), 0),
        ('h', ((1, 1), (0, 1)), 0),
        ('m', ((1, 0), (0, 1)), 0),
        ('n', ((0, 1), (0, 2)), -1),
    )
    result = []
    for parameter, keys, current in specifications:
        inputs = tuple(family.representative(key) for key in keys)
        outputs = tuple(family.INDEX[family.move(family.STATES[i], current)] for i in inputs)
        baseline = sum((family.PROBABILITY[key, current] for key in keys), F(0))/4
        if len(set(inputs)) != 2 or len(set(outputs)) != 2:
            raise ArithmeticError('each tomography probe needs two distinct states')
        for phase in ((0,) if parameter == 'm' else (0, 1)):
            name = parameter+('_re' if phase == 0 else '_im')
            result.append(Probe(name, parameter, phase, current, keys, inputs, outputs,
                                (baseline.numerator, baseline.denominator)))
    if tuple(probe.name for probe in result) != COORDINATES:
        raise ArithmeticError('the complete nine-coordinate design changed')
    return tuple(result)


PROBES = _probes()


def _fraction(value):
    if type(value) is int:
        return F(value)
    if type(value) is not F:
        raise ValueError('original integers or canonical Fractions are required')
    numerator, denominator = value.numerator, value.denominator
    if (type(numerator) is not int or type(denominator) is not int or denominator <= 0
            or gcd(numerator, denominator) != 1):
        raise ValueError('corrupt or noncanonical rational slots are not exact input')
    return F(numerator, denominator)


def _scalar(value):
    if type(value) in (int, F):
        return family.S(_fraction(value))
    if type(value) is not family.S:
        raise ValueError('an original exact Q(sqrt(2))(i) scalar is required')
    return family.S(*(_fraction(part) for part in value.parts()))


def _real(value):
    result = _scalar(value)
    if not result.is_real():
        raise ValueError('an exact real scalar is required')
    return result


def _unit_interval(value):
    result = _real(value)
    if result.sign() < 0 or (result-1).sign() > 0:
        raise ValueError('an exact probability in [0,1] is required')
    return result


def validate_table(table, require_cp=True):
    """Return a fresh exact table, checking the full existing CP criterion."""
    if type(require_cp) is not bool:
        raise ValueError('an original Boolean CP-validation selector is required')
    if (type(table) is not dict or any(type(key) is not str for key in table)
            or set(table) != set(PARAMETERS)):
        raise ValueError('the five original named table coordinates are required')
    result = {parameter: _scalar(table[parameter]) for parameter in PARAMETERS}
    if not result['m'].is_real():
        raise ValueError('endpoint covariance requires real m')
    if require_cp and not family.positive_semidefinite(result):
        raise ValueError('the reconstructed table is not a completely positive instrument')
    return result


def _probe(value):
    def same_original(left, right):
        if type(left) is not type(right):
            return False
        if type(right) in (Probe, tuple):
            return len(left) == len(right) and all(same_original(a, b) for a, b in zip(left, right))
        return left == right

    if type(value) is str:
        for item in PROBES:
            if item.name == value:
                return item
    elif type(value) is Probe:
        for item in PROBES:
            if same_original(value, item):
                return item
    raise ValueError('one of the nine original probe descriptors is required')


def _prediction(table, probe):
    value = table[probe.parameter]
    coordinate = family.S(value.a, value.b) if probe.phase == 0 else family.S(value.c, value.d)
    return family.S(F(*probe.baseline))+coordinate*F(1, 2)


def probability(table, probe):
    """Unconditional P(current,plus); phase 1 prepares relative amplitude i.

    The output projector uses the two current-translated input states.
    Every other current/readout outcome remains part of the instrument.
    """
    return _prediction(validate_table(table), _probe(probe))


def predictions(table):
    table = validate_table(table)
    return tuple(_prediction(table, probe) for probe in PROBES)


def reconstruct(probabilities):
    """Invert all nine probabilities and reject an inadmissible exact table.

    This is an exact inverse, not a CP projection of noisy empirical means.
    Error intervals must retain their feasible table set separately.
    """
    if type(probabilities) is not tuple or len(probabilities) != len(PROBES):
        raise ValueError('an original tuple of all nine exact probabilities is required')
    coordinates = {probe.name: (_unit_interval(value)-F(*probe.baseline))*2
                   for probe, value in zip(PROBES, probabilities, strict=True)}
    table = {}
    for parameter in PARAMETERS:
        real = coordinates[parameter+'_re']
        imaginary = family.S() if parameter == 'm' else coordinates[parameter+'_im']
        table[parameter] = family.S(real.a, real.b, imaginary.a, imaginary.b)
    return validate_table(table)


def apply_current(table, current, operator):
    """Apply a complete common-current branch to an exact sparse operator.

    Matrix-unit and non-Hermitian inputs are admitted for complex-linear
    identity checks. This function does not assert positivity of its input.
    """
    table = validate_table(table)
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('an original native current -1,0,1 is required')
    if type(operator) is not dict:
        raise ValueError('an original sparse operator dictionary is required')
    saved = {}
    for unit, value in operator.items():
        if (type(unit) is not tuple or len(unit) != 2
                or any(type(i) is not int or not 0 <= i < len(family.STATES) for i in unit)):
            raise ValueError('original carrier matrix-unit addresses are required')
        saved[unit] = _scalar(value)
    coefficients = family.multipliers(table)
    result = {}
    for (left, right), value in saved.items():
        x, y = family.STATES[left], family.STATES[right]
        coefficient = coefficients.get((current, family.key(x), family.key(y)))
        if coefficient:
            target = (family.INDEX[family.move(x, current)],
                      family.INDEX[family.move(y, current)])
            result[target] = result.get(target, family.S())+coefficient*value
    return {unit: value for unit, value in result.items() if value}


def _absolute_upper(value):
    # |a+b sqrt2+i(c+d sqrt2)| <= |a|+|c|+3(|b|+|d|)/2.
    # This exact rational bound avoids rounding or choosing a numeric root.
    return abs(value.a)+abs(value.c)+F(3, 2)*(abs(value.b)+abs(value.d))


def instrument_bound(first, second):
    """Reference-uniform half-diamond upper bound for common-current outputs.

    Write P_a for the six orthogonal seam-key projectors, with arbitrary
    spectator/reference factors. A positive unit-trace input obeys
    ||P_a rho P_b||_1 <= sqrt(p_a p_b). A same-current translation is an
    isometry on its supported input keys. The CQ output difference is thus
    at most sum_(j,a<b) |Delta G_j(a,b)| sqrt(p_a p_b), bounded by half the
    largest absolute row sum. This also bounds the channel diamond distance
    since the difference is Hermiticity preserving. Fine record refinements
    and coherent access to Kraus labels are outside this output alphabet.
    """
    first, second = validate_table(first), validate_table(second)
    k, g, h, m, n = (_absolute_upper(first[p]-second[p]) for p in PARAMETERS)
    return min(F(1), max(k+2*g, k+2*h, g+h+m+n, n)/2)


def coordinate_error_bound(delta):
    """Two admissible tables with nine probability differences <=delta.

    Complex coordinate differences have real-plus-imaginary bound 4delta;
    m has bound 2delta. The maximum row sum is therefore at most 14delta.
    A non-CP inverse of noisy point estimates is not an admissible comparator.
    """
    delta = _fraction(delta)
    if not 0 <= delta <= 1:
        raise ValueError('a probability error in [0,1] is required')
    return min(F(1), 7*delta)


def history_bound(call_count, per_call_error, prep_error=0, terminal=0):
    """Complete-output bound with a certified worst-branch native call cap.

    Both histories use identical intervening quantum instruments and
    classical feedback on the common current records. They retain every
    stopping/timeout branch and arbitrary untouched references. Local call
    errors telescope; preparation and final-service errors are separate.
    This does not infer an expected stopping cost or a physical calibration.
    """
    if type(call_count) is not int or call_count < 0:
        raise ValueError('an original nonnegative integer call cap is required')
    errors = tuple(_fraction(value) for value in (per_call_error, prep_error, terminal))
    if any(not 0 <= value <= 1 for value in errors):
        raise ValueError('complete-output errors must lie in [0,1]')
    return min(F(1), call_count*errors[0]+errors[1]+errors[2])


def controls():
    """Admissible blind directions for the old k/g tests and binary domain."""
    result = {name: family.extract_table(family.INSTRUMENTS[name]())
              for name in ('C0', 'C', 'E', 'T', 'C_signed')}
    result.update({
        'C_phase_quarter': family.table(1, family.S(0, 0, -F(1, 2)),
                                       family.S(0, 0, -F(1, 2)), 0, F(1, 2)),
        'T_h_quarter': family.table(0, 0, F(1, 4), 0, 0),
        'T_h_imag': family.table(0, 0, family.S(0, 0, F(1, 4)), 0, 0),
        'T_m_quarter': family.table(0, 0, 0, F(1, 4), 0),
        'T_n_quarter': family.table(0, 0, 0, 0, F(1, 4)),
        'E_n_quarter': family.table(1, 0, 0, 0, F(1, 4)),
    })
    return {name: validate_table(table) for name, table in result.items()}


def rank_certificate():
    """Rank nine and minimality among scalar affine probability statistics.

    All nine real coordinates vary independently in a cube of radius 1/16
    about the zero table. Each G block is strictly diagonally dominant there
    with margin >=3/16. Thus positivity does not reduce the local dimension.
    This is not a lower bound for joint multi-outcome settings or a finite
    preselected list of candidates.
    """
    zero = family.table(0, 0, 0, 0, 0)
    base = predictions(zero)
    radius = F(1, 16)
    columns = []
    for probe in PROBES:
        table = validate_table(zero)
        table[probe.parameter] = family.S(radius) if probe.phase == 0 else family.S(0, 0, radius)
        column = []
        for value, initial in zip(predictions(table), base, strict=True):
            change = (value-initial)*16
            if change.b or change.c or change.d:
                raise ArithmeticError('the rational probe Jacobian changed field')
            column.append(F(change.a))
        columns.append(tuple(column))
    response = tuple(tuple(columns[j][i] for j in range(9)) for i in range(9))
    if any(response[i][j] != F(int(i == j), 2) for i in range(9) for j in range(9)):
        raise ArithmeticError('the complete nine-coordinate inverse lost rank')
    return {'coordinates': COORDINATES, 'jacobian': response, 'rank': 9,
            'determinant': F(1, 512), 'interior_coordinate_radius': radius,
            'interior_diagonal_dominance_margin': F(3, 16),
            'minimum_scalar_affine_statistics': 9}
