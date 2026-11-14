"""Independently bind literal port controls to exact full six-level response.

Actual integrated_controls outputs are interpreted as their exact binary64
rationals. Generator entries use the retained Q(sqrt(5)) source law, not
the floating expm matrices. Every complement export is kept in the channel.
"""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path

from ..observer.source_scalar import compiler
from .interval import Interval, sqrt_fraction

ROOT = Path(__file__).resolve().parents[2]
ZERO = (F(), F())


def plus(a, b):
    return a[0]+b[0], a[1]+b[1]


def multiply(a, b):
    return a[0]*b[0]+5*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def scale(a, b):
    return a[0]*b, a[1]*b


def field_interval(value):
    return Interval.point(value[0])+value[1]*sqrt_fraction(F(5))


def _target(name):
    result = [[(F(), F()) for _ in range(6)] for _ in range(6)]
    if name.startswith('phase_'):
        i = int(name[-1]); result[i][i] = F(), F(1)
    elif name == 'relative_phase':
        for i in range(3): result[i][i] = F(), F(1)
    else:
        i, j = map(int, name[-2:])
        if name.startswith('real_'):
            result[i][j], result[j][i] = (F(1), F()), (F(-1), F())
        elif name.startswith('imaginary_'):
            result[i][j] = result[j][i] = F(), F(1)
        else:
            raise ValueError('unknown exact base source control')
    return result


@lru_cache(maxsize=1)
def specification():
    packet = json.loads((ROOT/compiler.RESPONSE).read_text())
    full = []
    if len(packet['generators']) != 12:
        raise ArithmeticError('the source must retain all twelve response ports')
    for blocks in packet['generators']:
        matrix = [[(ZERO, ZERO) for _ in range(6)] for _ in range(6)]
        for b, block in enumerate(blocks):
            for i, row in enumerate(block):
                for j, value in enumerate(row):
                    if len(value) != 4 or any(type(x) is not str for x in value):
                        raise ValueError('exact Q(sqrt(5)) source entry required')
                    a, c, d, e = map(F, value)
                    matrix[3*b+i][3*b+j] = ((a, c), (d, e))
        for i in range(6):
            for j in range(6):
                a, b = matrix[i][j]; c, d = matrix[j][i]
                if plus(a, c) != ZERO or b != d:
                    raise ArithmeticError('the exact source is not anti-Hermitian')
        full.append(matrix)
    rows = json.loads((ROOT/compiler.CONTROL_RECEIPT).read_text())['evidence']['controls']['rows']
    expected = ['phase_0', 'phase_1', 'phase_2', 'real_01', 'imaginary_01',
                'real_02', 'imaginary_02', 'real_12', 'imaginary_12', 'relative_phase']
    if [r['name'] for r in rows] != expected:
        raise ArithmeticError('all ten exact control identities must be present')
    for row in rows:
        coefficients = [tuple(map(F, c)) for c in row['coefficients']]
        if len(coefficients) != 12:
            raise ArithmeticError('all twelve control coefficients required')
        target = _target(row['name'])
        for i in range(6):
            for j in range(6):
                for part in (0, 1):
                    result = ZERO
                    for c, g in zip(coefficients, full, strict=True):
                        result = plus(result, multiply(c, g[i][j][part]))
                    if result != (target[i][j][part], F()):
                        raise ArithmeticError('retained exact control identity failed')
    return tuple(tuple(tuple(row) for row in matrix) for matrix in full)


def literal_controls(kind, angle):
    values = compiler.integrated_controls(kind, angle)
    if len(values) != 12 or not all(float(x) == x for x in values):
        raise ValueError('twelve finite binary64 integrated port controls required')
    return tuple(F(float(x)) for x in values)


@lru_cache(maxsize=32768)
def generator_error(kind, angle, controls):
    """Max-row bound on G_actual-angle*T_kind, hence on unitary error.

    Both matrices are anti-Hermitian. Absolute row and column bounds agree,
    so the row bound controls spectral norm. Duhamel requires no e^||G||.
    """
    angle = F(angle)
    if kind not in ('X', 'Y', 'Z', 'I', 'J', 'H') or len(controls) != 12:
        raise ValueError('a native kind with all twelve exact controls required')
    if any(type(x) is not F for x in controls):
        raise ValueError('literal controls must be exact rational values')
    spec = specification()
    target = [[(Interval.point(0), Interval.point(0)) for _ in range(6)] for _ in range(6)]
    z, a = Interval.point(0), Interval.point(angle)
    if kind == 'X': target[0][1] = target[1][0] = z, a
    if kind == 'Y': target[0][1], target[1][0] = (a, z), (-a, z)
    if kind in ('Z', 'I'):
        target[0][0] = z, a; target[1][1] = z, a if kind == 'I' else -a
    if kind == 'J':
        for i in range(3): target[i][i] = z, a
    if kind == 'H':
        b = angle*sqrt_fraction(F(1, 2))
        target[0][0] = target[0][1] = target[1][0] = z, b
        target[1][1] = z, -b
    maximum = F()
    for i in range(6):
        row_bound = F()
        for j in range(6):
            for part in (0, 1):
                result = ZERO
                for c, g in zip(controls, spec, strict=True):
                    result = plus(result, scale(g[i][j][part], c))
                row_bound += (field_interval(result)-target[i][j][part]).absolute_upper()
        maximum = max(maximum, row_bound)
    return maximum


def pulse_error(event):
    controls = literal_controls(event['kind'], event['angle'])
    return generator_error(event['kind'], F(event['angle']), controls)


def complete_export_certificate():
    """Literal coherent/dephasing exports, including all outside-code rows."""
    from fractions import Fraction
    for dimension in (2, 4):
        for coherent in (False, True):
            kraus = compiler.transfer_kraus(dimension, coherent)
            expected_count = (1 if coherent else dimension)+6-dimension
            if len(kraus) != expected_count:
                raise ArithmeticError('complete export Kraus branch was dropped')
            gram = [[Fraction() for _ in range(6)] for _ in range(6)]
            for k in kraus:
                if k.shape != (dimension+1, 6):
                    raise ArithmeticError('export carrier changed')
                for row in k:
                    if any(x not in (0, 1) for x in row):
                        raise ArithmeticError('declared exact export contains another scalar')
                    for i in range(6):
                        for j in range(6): gram[i][j] += F(int(row[i].real))*F(int(row[j].real))
            if gram != [[F(i == j) for j in range(6)] for i in range(6)]:
                raise ArithmeticError('full source export is not trace preserving')
    return {'proper_code_dimensions': [2, 4], 'coherent_and_dephasing_export_checked': True,
            'complement_Kraus_outputs_retained': True,
            'complete_export_isometry_preserves_operator_error': True}
