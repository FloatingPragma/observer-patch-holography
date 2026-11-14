"""Additive ten-bit diagonal compiler with outward fast-Walsh phase targets."""
from fractions import Fraction as F
import math
import sys

sys.dont_write_bytecode = True
import providers

nf = providers.native()


def transform(values):
    """Unnormalized Walsh transform, using only exact interval additions."""
    result = list(values); n = len(result)
    if n < 2 or n & (n-1) or n > 1024:
        raise ValueError('a complete power-of-two phase target through ten bits required')
    step = 1
    while step < n:
        for start in range(0, n, 2*step):
            for offset in range(step):
                i, j = start+offset, start+offset+step
                left, right = result[i], result[j]
                result[i], result[j] = left+right, left-right
        step *= 2
    return result


def compile_intervals(phases):
    n = len(phases); bits = n.bit_length()-1
    coefficients = [x/n for x in transform(phases)]
    word = []
    for mask, value in enumerate(coefficients):
        angle = float(F(value.lo+value.hi, 2*nf.iv.SCALE))
        label = ''.join('Z' if mask >> j & 1 else 'I' for j in range(bits))
        word.extend(nf.compiler.pauli_rotation(label, angle))
    assert len(word) == (3*bits-5)*n+6
    return word, certify(word, phases)


def certify(word, phases):
    n = len(phases); bits = n.bit_length()-1
    targets = [x/n for x in transform(phases)]
    offset = 0; error = F()
    for mask, target in enumerate(targets):
        label = ''.join('Z' if mask >> j & 1 else 'I' for j in range(bits))
        offset, fixed, variable = nf.cert.consume(word, offset, nf.cert.pauli(label))
        assert len(variable) == 1
        error += fixed+(nf.iv.Interval.point(variable[0])-target).absolute_upper()
    assert offset == len(word)
    return error+nf.cert.controls_bound(word, {})


def controls():
    # Exact integer butterflies versus an independent direct character sum.
    for bits in (1, 2, 4, 8, 10):
        n = 1 << bits; values = [((j*j+17*j+3) % 31)-15 for j in range(n)]
        actual = transform(values)
        for mask in range(n):
            expected = sum(x*(-1)**((j & mask).bit_count()) for j, x in enumerate(values))
            assert actual[mask] == expected
        assert transform(actual) == [n*x for x in values]
    for n in (1, 3, 2048):
        try: transform([0]*n)
        except ValueError: pass
        else: raise AssertionError('unsupported Walsh domain admitted')
    return {'integer_character_sums_through10bits':True, 'transform_square_n_identity':True,
            'all_scalar_and_zero_coefficients_retained':True, 'unsupported_domains_rejected':True}
