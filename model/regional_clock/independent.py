"""Independent sector target and rational literal-angle certificates.

This module never calls the producer's diagonal, phase, log or Walsh helpers.
The declared native encoding is an interface input, not a source of weights.
Binary64 matrix residuals are separate diagnostics. Analytic certificates use
exact Q arithmetic and the source's checked exact generator identities.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import math
import json
from pathlib import Path

from ..observer.source_discrimination import circuit


STATES = tuple(sorted(tuple(int(p in pair) for p in range(12))
                      for pair in combinations(range(12), 2)))
CODES = tuple(circuit.encode(x) for x in STATES)


def units(ports):
    """All matrix units of direct_sum_n M_choose(|U|,n) tensor I_ext,n."""
    outside = tuple(p for p in range(12) if p not in ports)
    index = {tuple(i for i, bit in enumerate(x) if bit): j for j, x in enumerate(STATES)}
    for n in range(3):
        for a in combinations(ports, n):
            for b in combinations(ports, n):
                support = tuple((index[tuple(sorted((*a, *e)))],
                                 index[tuple(sorted((*b, *e)))])
                                for e in combinations(outside, 2-n))
                if support:
                    yield n, a, b, support


def target(values, ports):
    w = tuple(F(x)/F(values[7]) for x in values)
    return tuple(math.prod(w[p] for p in ports if x[p]) for x in STATES)


def reduced_reference(values, ports):
    """Actual normalized reduction, with all local-number sectors retained."""
    products = tuple(math.prod(F(w)**n for w, n in zip(values, x, strict=True)) for x in STATES)
    z = sum(products)
    reduced = {}
    for x, probability in zip(STATES, products, strict=True):
        local = tuple(p for p in ports if x[p])
        reduced[local] = reduced.get(local, F())+probability/z
    return reduced


def algebra_certificate(values, ports, ratios=None):
    expected = target(values, ports)
    ratios = expected if ratios is None else tuple(ratios)
    if len(ratios) != 66 or any(x <= 0 for x in ratios):
        raise ValueError('positive ratio on every H66 state required')
    reduced = reduced_reference(values, ports)
    counts = [0, 0, 0]
    for n, a, b, support in units(ports):
        counts[n] += 1
        ratio = reduced[a]/reduced[b]
        if any(ratios[i]/ratios[j] != ratio for i, j in support):
            raise ArithmeticError('wrong reference, phase sign, or omitted regional sector')
    exterior = tuple(p for p in range(12) if p not in ports)
    exterior_units = 0
    for _, _, _, support in units(exterior):
        exterior_units += 1
        if any(ratios[i] != ratios[j] for i, j in support):
            raise ArithmeticError('actuator disturbs the exterior algebra')
    return {'regional_matrix_units_by_number': counts,
            'exterior_matrix_units': exterior_units,
            'full_H66_and_arbitrary_reference': True}


def expectation_diagonal(values, ports, observable):
    """Exact weighted conditional expectation on diagonal input operators."""
    probabilities = [math.prod(F(w)**n for w, n in zip(values, x, strict=True)) for x in STATES]
    groups = {}
    for x, p, f in zip(STATES, probabilities, observable, strict=True):
        key = tuple(x[j] for j in ports)
        a, b = groups.get(key, (F(), F()))
        groups[key] = a+p*f, b+p
    return tuple(groups[tuple(x[j] for j in ports)][0]/groups[tuple(x[j] for j in ports)][1]
                 for x in STATES)


@lru_cache(maxsize=256)
def log_interval(value):
    """Independent atanh series, exact outward remainder after 48 terms."""
    value = F(value)
    if value <= 0:
        raise ValueError('positive logarithm argument required')
    k = value.numerator.bit_length()-value.denominator.bit_length()
    r = value/F(2)**k
    if r < 1:
        r *= 2
        k -= 1
    def local(q):
        x = (q-1)/(q+1)
        lo = 2*sum((x**(2*j+1)/F(2*j+1) for j in range(48)), F())
        return lo, lo+2*x**97/(97*(1-x*x))
    a, b = local(r)
    c, d = sorted(k*t for t in local(F(2)))
    return a+c, b+d


@lru_cache(maxsize=1)
def pi_interval():
    # Machin: pi = 16 atan(1/5) - 4 atan(1/239). Alternating remainders.
    def atan(q):
        lo = sum(((-1)**j*q**(2*j+1)/F(2*j+1) for j in range(96)), F())
        return lo, lo+q**193/193
    a, b = atan(F(1, 5))
    c, d = atan(F(1, 239))
    return 16*a-4*d, 16*b-4*c


def native_control_error(tape):
    """Rigorous finite-control rounding cost against exact native D_p.

    The source executor emits finite binary64 integrated port coefficients.
    Reconstruct sum_p c_p D_p exactly in Q(sqrt(5)), then compare it with the
    intended ideal code generator using rational sqrt enclosures. The entry
    L1 majorant bounds the full M6 operator norm, including complement rows.
    Duhamel applies because both generators are antihermitian. This is not an
    error estimate for scipy.expm or an assumed physical control precision.
    """
    from ..observer.source_scalar import compiler
    root = Path(__file__).resolve().parents[2]
    raw = json.loads((root/compiler.RESPONSE).read_text())['generators']
    generators = []
    for blocks in raw:
        entries = {}
        for block_index, block in enumerate(blocks):
            for i, row in enumerate(block):
                for j, cell in enumerate(row):
                    entries[3*block_index+i, 3*block_index+j] = tuple(F(x) for x in cell)
        generators.append(entries)
    scale = 1 << 192
    def sqrt_box(n):
        lo = math.isqrt(n*scale*scale)
        return F(lo, scale), F(lo+1, scale)
    sqrt5, sqrt2 = sqrt_box(5), sqrt_box(2)
    def field_box(a, b):
        low, high = sorted(b*x for x in sqrt5)
        return a+low, a+high
    total, cached = F(), {}
    for pulse in tape:
        kind, angle = pulse['kind'], F(pulse['angle'])
        key = kind, angle
        if key not in cached:
            c = tuple(F(float(x)) for x in compiler.integrated_controls(kind, float(angle)))
            error = F()
            for i in range(6):
                for j in range(6):
                    actual = tuple(sum((c[k]*generators[k].get((i, j), (F(),)*4)[part]
                                        for k in range(12)), F()) for part in range(4))
                    real, imag = F(), F()
                    if kind == 'X' and (i, j) in ((0, 1), (1, 0)):
                        imag = angle
                    elif kind == 'Y' and (i, j) in ((0, 1), (1, 0)):
                        real = angle*(1 if i == 0 else -1)
                    elif kind == 'Z' and i == j and i < 2:
                        imag = angle*(1 if i == 0 else -1)
                    elif kind == 'I' and i == j and i < 2:
                        imag = angle
                    elif kind == 'J' and i == j and i < 3:
                        imag = angle
                    re_target, im_target = (real, real), (imag, imag)
                    if kind == 'H' and i < 2 and j < 2:
                        sign = -1 if i == j == 1 else 1
                        im_target = tuple(sorted(sign*angle/x for x in sqrt2))
                    for observed, target_box in ((field_box(actual[0], actual[1]), re_target),
                                                 (field_box(actual[2], actual[3]), im_target)):
                        error += max(abs(observed[0]-target_box[1]), abs(observed[1]-target_box[0]))
            cached[key] = error
        total += cached[key]
    return total


def angle_certificate(tape, values, ports, parameter=F(1)):
    """Prove the emitted literal word's operator error on ALL512 columns.

    Replace each recognized H/J/H triple by its exact CNOT, charging each
    angle displacement by unit-norm-generator Duhamel and telescoping.
    Symbolically propagate every address through the resulting permutation
    circuit and sum its exact dyadic I/Z phase. Compare with independently
    reconstructed -s log(product_U w). No numerical exponential is trusted.
    """
    pi_lo, pi_hi = pi_interval()
    positions = list(range(512))
    phases = [F() for _ in positions]
    cnot_error, cnots, k = F(), 0, 0
    while k < len(tape):
        event = tape[k]
        kind, angle, targets = event['kind'], F(event['angle']), event['targets']
        if kind in ('I', 'Z') and len(targets) == 1 and 0 <= targets[0] < 9:
            for j, address in enumerate(positions):
                phases[j] += angle*(-1 if kind == 'Z' and address >> targets[0] & 1 else 1)
            k += 1
        else:
            triple = tape[k:k+3]
            if (len(triple) != 3 or [x['kind'] for x in triple] != ['H', 'J', 'H']
                    or len(triple[1]['targets']) != 2
                    or triple[0]['targets'] != triple[2]['targets']
                    or triple[0]['targets'] != triple[1]['targets'][1:]):
                raise ArithmeticError('unrecognized source parity circuit')
            control, dest = triple[1]['targets']
            if not 0 <= control < 9 or not 0 <= dest < 9 or control == dest:
                raise ArithmeticError('invalid source parity targets')
            for pulse, a, b in zip(triple, (pi_lo/2, -pi_hi, pi_lo/2),
                                  (pi_hi/2, -pi_lo, pi_hi/2), strict=True):
                cnot_error += max(abs(F(pulse['angle'])-a), abs(F(pulse['angle'])-b))
            positions = [j ^ ((1 << dest) if j >> control & 1 else 0) for j in positions]
            cnots += 1
            k += 3
    if positions != list(range(512)):
        raise ArithmeticError('parity circuit does not close on all addresses')
    ratios = dict(zip(CODES, target(values, ports), strict=True))
    phase_error = F()
    for address, phase in enumerate(phases):
        a, b = sorted(-parameter*x for x in log_interval(ratios.get(address, F(1))))
        phase_error = max(phase_error, abs(phase-a), abs(phase-b))
    controls_error = native_control_error(tape)
    total = phase_error+cnot_error+controls_error
    if total > F(1, 10**9):
        raise ArithmeticError('literal regional word exceeds the fixed analytic cap')
    return {'all_address_columns': 512, 'ideal_CNOT_count': cnots,
            'binary64_diagonal_angle_error_upper': str(phase_error),
            'binary64_CNOT_angle_error_upper': str(cnot_error),
            'binary64_integrated_controls_error_upper': str(controls_error),
            'analytic_operator_error_upper': str(total),
            'analytic_halfdiamond_error_upper': str(min(F(1), total)),
            'ideal_control_processor_complement_probability': '0',
            'finite_control_first_processor_failure_probability_upper': str(min(F(1), controls_error**2)),
            'exact_code_generator_error': '0',
            'full_address_padding_included': True}
