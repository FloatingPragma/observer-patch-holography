"""Rational error certificate for a public two-reading event on the emitted cell."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial, isqrt

from model.executable_model.carrier import CELL_SEAMS

SCALE = 1 << 48
DEGREE = 16


def pi_interval():
    """Machin's identity, with alternating-series bounds on both arctangents."""
    def atan_bounds(x):
        partial = sum((-1)**j * x**(2*j+1)/(2*j+1) for j in range(8))
        return partial, partial + x**17/17

    a, b = atan_bounds(F(1, 5))
    c, d = atan_bounds(F(1, 239))
    return 16*a - 4*d, 16*b - 4*c


def generators():
    """Twice the cell Hamiltonian in the emitter energy frame, on fuel + 66 binary pairs."""
    pairs = list(combinations(range(12), 2))
    index = {pair: i+1 for i, pair in enumerate(pairs)}
    adjacent = [index[edge] for edge in CELL_SEAMS]
    port_pair = index[(0, 1)]
    # H = K + 8 (12 - 2 adjacency). At the emitting seam Ef = 4 + 8*10.
    rows = [{} for _ in range(1+len(pairs))]
    for pair, i in index.items():
        neighbors = []
        for a, b in CELL_SEAMS:
            if (a in pair) != (b in pair):
                changed = set(pair) ^ {a, b}
                neighbors.append(index[tuple(sorted(changed))])
        rows[i] = {j: -1 for j in neighbors}
        rows[i][i] = len(neighbors) + 16*(12 - 2*(i in adjacent)) - 2*84
    pulse = [row.copy() for row in rows]
    pulse[0][port_pair] = pulse[port_pair][0] = 128
    return rows, pulse, adjacent, port_pair


def evolve(rows, duration, seed, error=F(0)):
    """Fixed-point Taylor actions with a rational Euclidean error radius.

    Input and output components are complex integer pairs in units of 2^-48.
    The reference exact unitary starts from a unit vector within the input radius.
    """
    twice_norm = max(sum(map(abs, row.values())) for row in rows)
    steps = 1
    while abs(duration)*twice_norm > steps:
        steps *= 2
    step = duration/(2*steps)
    num, den = step.numerator, step.denominator
    edges = [tuple(row.items()) for row in rows]
    real, imag = (list(part) for part in seed)
    rounding = F(isqrt(2*len(rows))+1, SCALE)
    tail = F(2, 2**(DEGREE+1)*factorial(DEGREE+1))
    # Each computed starting vector has norm <= 2 by induction using this bound.
    total_error = error + steps*(2*tail + 2*DEGREE*rounding)
    if total_error >= 1:
        raise ArithmeticError('Taylor norm bootstrap failed')
    for _ in range(steps):
        term_r, term_i = real.copy(), imag.copy()
        for j in range(1, DEGREE+1):
            divisor = den*j
            next_r = [num*sum(v*term_i[k] for k, v in row)//divisor for row in edges]
            next_i = [-num*sum(v*term_r[k] for k, v in row)//divisor for row in edges]
            term_r, term_i = next_r, next_i
            real = [a+b for a, b in zip(real, term_r)]
            imag = [a+b for a, b in zip(imag, term_i)]
    return (real, imag), total_error


def event_certificate():
    """Intervals for first presence followed by a positive record on the emitter seam."""
    free, pulse, seams, emitting = generators()
    size = len(free)

    def basis(index):
        return [SCALE*int(i == index) for i in range(size)], [0]*size

    low_pi, high_pi = pi_interval()
    pi_center = F(355, 113)
    twice_norm = max(sum(map(abs, row.values())) for row in pulse)
    pulse_error = max(abs(pi_center-low_pi), abs(pi_center-high_pi))*twice_norm/256
    emitted, error = evolve(pulse, pi_center/128, basis(0), pulse_error)
    (ar, ai), error_a = evolve(free, F(4), emitted, error)
    (br, bi), error_b = evolve(free, F(-4), basis(emitting))
    real = [br[i]*ar[i] + bi[i]*ai[i] for i in seams]
    imag = [br[i]*ai[i] - bi[i]*ar[i] for i in seams]
    coarse = F(sum(real)**2 + sum(imag)**2, SCALE**4)
    resolved = F(sum(r*r + j*j for r, j in zip(real, imag)), SCALE**4)
    amplitude_error = error_a + error_b + error_a*error_b
    probability_error = amplitude_error*(2+amplitude_error)
    lower = abs(resolved-coarse) - 2*probability_error - F(48, 16384)
    return {'aggregate_event': coarse, 'resolved_event': resolved,
            'each_probability_error': probability_error, 'always_on_tv_lower': lower}


if __name__ == '__main__':
    certificate = event_certificate()
    for name, value in certificate.items():
        print(f'{name} = {value} (approximately {float(value):.10g})')
