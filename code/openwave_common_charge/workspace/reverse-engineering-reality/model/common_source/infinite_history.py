"""Certified interacting history on two infinite oscillators.

Finite polynomials approximate states, not the Hamiltonian or final detector.
Each time step applies the full infinite Hamiltonian to expanding support,
charges its spectral Taylor remainder, then charges all discarded and rounded
amplitude. Coefficients use unnormalized Fock monomials with factorial norms.
No quantum branch is normalized, selected, or replaced by a target probability.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
import math


Key = tuple[int, int, int]
GaussianInteger = tuple[int, int]
ROOT_SCALE = 10**30


@dataclass
class _Polynomial:
    coefficients: dict[Key, GaussianInteger]
    denominator: int = 1


def _add(target, key, real, imaginary):
    old = target.get(key, (0, 0))
    value = old[0]+real, old[1]+imaginary
    if value == (0, 0):
        target.pop(key, None)
    else:
        target[key] = value


def _scaled_add(target, key, value, factor):
    _add(target, key, value[0]*factor, value[1]*factor)


def _hamiltonian_numerator(vector, coupling=F(1, 4)):
    """Apply 128 H exactly, with no support projection between powers.

    H=X0+X1+E^2+(73/64)(N0+N1+1)
      +(9/128) sum_i(a_i^2+a_i^{*2})
      -(9/64)(a0+a0*)(a1+a1*)+(1/4)(a0+a0*) Z0.
    Creation shifts a monomial; annihilation also multiplies by occupation.
    """
    if type(coupling) is not F or coupling not in (F(0), F(1, 4)):
        raise ValueError('exact declared coupling 1/4 or disconnected control 0 required')
    out = {}
    for (n0, n1, gauge), value in vector.items():
        b0, b1 = gauge & 1, gauge >> 1
        _scaled_add(out, (n0, n1, gauge), value,
                    146*(n0+n1+1)+128*(b1-b0)**2)
        for bit in (1, 2):
            _scaled_add(out, (n0, n1, gauge ^ bit), value, 128)
        _scaled_add(out, (n0+2, n1, gauge), value, 9)
        _scaled_add(out, (n0, n1+2, gauge), value, 9)
        if n0 >= 2:
            _scaled_add(out, (n0-2, n1, gauge), value, 9*n0*(n0-1))
        if n1 >= 2:
            _scaled_add(out, (n0, n1-2, gauge), value, 9*n1*(n1-1))
        for delta0, factor0 in ((1, 1), (-1, n0)):
            for delta1, factor1 in ((1, 1), (-1, n1)):
                if factor0 and factor1:
                    _scaled_add(out, (n0+delta0, n1+delta1, gauge), value,
                                -18*factor0*factor1)
        if coupling:
            _scaled_add(out, (n0+1, n1, gauge), value, 32*(1-2*b0))
            if n0:
                _scaled_add(out, (n0-1, n1, gauge), value, 32*n0*(1-2*b0))
    return out


def _norm_numerator(coefficients):
    return sum(math.factorial(n0)*math.factorial(n1)*(real*real+imag*imag)
               for (n0, n1, _), (real, imag) in coefficients.items())


def _root_upper(value):
    if type(value) is not F or value < 0:
        raise ValueError('a nonnegative exact rational squared norm required')
    n = math.isqrt(value.numerator*ROOT_SCALE**2//value.denominator)
    if F(n, ROOT_SCALE)**2 < value:
        n += 1
    return F(n, ROOT_SCALE)


def _integer(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f'{name} must be an integer in [{low}, {high}]')
    return value


def _initial(maximum):
    """Projection of exp[-(a0* - a1*)^2/36]|vacuum>, before normalization."""
    coefficients = {}
    denominator = 36**(maximum//2)*math.factorial(maximum//2)
    for pair in range(maximum//2+1):
        factor = denominator//(36**pair*math.factorial(pair))
        for n0 in range(2*pair+1):
            n1 = 2*pair-n0
            coefficients[n0, n1, 0] = (
                factor*(-1)**(pair+n1)*math.comb(2*pair, n0), 0)
    return _Polynomial(coefficients, denominator)


def _taylor(vector, time, degree, coupling=F(1, 4)):
    """Exact polynomial and exact next Taylor term, both on full support."""
    term = dict(vector.coefficients)
    total = dict(term)
    denominator = vector.denominator
    for power in range(1, degree+2):
        term = {key: (time.numerator*imag, -time.numerator*real)
                for key, (real, imag) in _hamiltonian_numerator(term, coupling).items()}
        factor = 128*time.denominator*power
        denominator *= factor
        if power == degree+1:
            return _Polynomial(total, denominator//factor), _Polynomial(term, denominator)
        total = {key: (real*factor, imag*factor) for key, (real, imag) in total.items()}
        for key, (real, imag) in term.items():
            _add(total, key, real, imag)
    raise AssertionError('unreachable')


def _cut_and_round(vector, cutoff, bits, normalization):
    keep, removed = {}, {}
    for key, value in vector.coefficients.items():
        (keep if key[0]+key[1] <= cutoff else removed)[key] = value
    discarded = _root_upper(F(_norm_numerator(removed), vector.denominator**2)/normalization)
    grid = 1 << bits
    rounded, differences = {}, {}
    for key, (real, imag) in keep.items():
        r, i = ((2*part*grid+vector.denominator)//(2*vector.denominator)
                for part in (real, imag))
        _add(rounded, key, r, i)
        _add(differences, key, real*grid-r*vector.denominator,
             imag*grid-i*vector.denominator)
    rounding = _root_upper(F(_norm_numerator(differences),
                            (vector.denominator*grid)**2)/normalization)
    return _Polynomial(rounded, grid), discarded, rounding


def _advance(vector, time, steps, degree, cutoff, bits, normalization, coupling=F(1, 4)):
    errors = {'Taylor': F(0), 'discarded_support': F(0), 'coefficient_rounding': F(0)}
    maximum_support = 0
    for _ in range(steps):
        polynomial, remainder = _taylor(vector, time/steps, degree, coupling)
        maximum_support = max(maximum_support, max(
            (n0+n1 for n0, n1, _ in remainder.coefficients), default=0))
        errors['Taylor'] += _root_upper(F(_norm_numerator(remainder.coefficients),
                                        remainder.denominator**2)/normalization)
        vector, discarded, rounding = _cut_and_round(polynomial, cutoff, bits, normalization)
        errors['discarded_support'] += discarded
        errors['coefficient_rounding'] += rounding
    return vector, errors, maximum_support


def _first_record(vector, bit, enabled):
    selected = {key: value for key, value in vector.coefficients.items() if key[2] & 1 == bit}
    if not enabled or bit == 0:
        return _Polynomial(selected, vector.denominator)
    # The SAME recorded-bit feedback as the finite experiment: (3I-4iX0)/5.
    output = {}
    for (n0, n1, gauge), (real, imag) in selected.items():
        _add(output, (n0, n1, gauge), 3*real, 3*imag)
        _add(output, (n0, n1, gauge ^ 1), 4*imag, -4*real)
    return _Polynomial(output, 5*vector.denominator)


def _rational_vector(vector):
    return {key: (F(real, vector.denominator), F(imag, vector.denominator))
            for key, (real, imag) in vector.coefficients.items()}


def _interval(low, high, places=12):
    scale = 10**places
    return [str(F(low.numerator*scale//low.denominator, scale)),
            str(F(-((-high.numerator*scale)//high.denominator), scale))]


def certificate(cutoff=16, preparation_cutoff=12, degree=16, steps_per_unit=16,
                bits=112, coupling=F(1, 4), retain_polynomials=False):
    """Recompute an unconditional infinite-oscillator feedback contrast.

    Exact rational residuals establish the errors; float simulation and the
    finite d=2 signal are not inputs. The fixed native preparation current is
    zero. This certifies the declared infinite two-site action, not a spatial
    continuum, physical calibration, or derivation of that action from axioms.
    """
    from .halfline import positive_weight_interval, norm_squared
    cutoff = _integer(cutoff, 2, 64, 'retained total occupation')
    preparation_cutoff = _integer(preparation_cutoff, 0, cutoff, 'preparation cutoff')
    degree = _integer(degree, 1, 48, 'Taylor degree')
    steps_per_unit = _integer(steps_per_unit, 2, 128, 'steps per unit')
    bits = _integer(bits, 32, 256, 'coefficient precision')
    if type(coupling) is not F or coupling not in (F(0), F(1, 4)):
        raise ValueError('exact declared coupling 1/4 or disconnected control 0 required')
    if type(retain_polynomials) is not bool:
        raise ValueError('retained-polynomial output switch must be Boolean')
    if preparation_cutoff % 2 or steps_per_unit % 2:
        raise ValueError('even preparation cutoff and steps per unit required')
    initial = _initial(preparation_cutoff)
    normalization = F(_norm_numerator(initial.coefficients), initial.denominator**2)
    sqrt5_lower = F(math.isqrt(5*ROOT_SCALE**2), ROOT_SCALE)
    tail = 1-F(4, 9)*sqrt5_lower*normalization
    preparation_error = _root_upper(tail)
    first, first_errors, maximum = _advance(initial, F(1), steps_per_unit,
                                           degree, cutoff, bits, normalization, coupling)
    first_error = sum(first_errors.values(), F(0))
    configurations = []
    for enabled in (True, False):
        branches, totals, second_errors = [], [], []
        for bit in (0, 1):
            selected = _first_record(first, bit, enabled)
            state, errors, support = _advance(selected, F(1, 2), steps_per_unit//2,
                                              degree, cutoff, bits, normalization, coupling)
            maximum = max(maximum, support)
            second_errors.append(sum(errors.values(), F(0)))
            vector = _rational_vector(state)
            positive = positive_weight_interval(vector, normalization)
            mass = norm_squared(vector)/normalization
            totals.append(positive)
            rows = [[*key, str(real), str(imag)] for key, (real, imag) in sorted(vector.items())]
            digest = hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()
            branches.append({'first_record': bit,
                             'pre_readout_polynomial_sha256': digest,
                             'pre_readout_nonzero_coefficients': len(rows),
                             'approximate_pre_readout_norm_squared': str(mass),
                             'approximate_positive_weight': _interval(*positive),
                             'approximate_negative_weight': _interval(mass-positive[1], mass-positive[0]),
                             'second_segment_vector_errors': {k: str(v) for k, v in errors.items()}})
            if retain_polynomials:
                branches[-1]['pre_readout_polynomial'] = rows
        # Orthogonal retained-record branches combine in Euclidean norm.
        error = first_error+_root_upper(sum((x*x for x in second_errors), F(0)))
        probability_error = preparation_error+2*error+error*error
        lo, hi = (sum((x[i] for x in totals), F(0)) for i in (0, 1))
        configurations.append({
            'feedback_enabled': enabled, 'retained_first_records': branches,
            'normalized_preparation_to_polynomial_vector_error_upper': str(error),
            'complete_output_half_trace_norm_error_upper': str(preparation_error+error+error*error/2),
            'event_probability_error_upper': str(probability_error),
            'infinite_positive_probability_interval': _interval(max(F(0), lo-probability_error),
                                                               min(F(1), hi+probability_error)),
        })
    on, off = (tuple(map(F, row['infinite_positive_probability_interval'])) for row in configurations)
    contrast = _interval(max(F(-1), on[0]-off[1]), min(F(1), on[1]-off[0]))
    return {'schema': 'oph.common_source.infinite_history.v1',
            'parameters': {'retained_total_occupation': cutoff, 'initial_total_occupation': preparation_cutoff,
                           'Taylor_degree': degree, 'steps_per_unit': steps_per_unit, 'coefficient_bits': bits,
                           'coupling': str(coupling)},
            'initial_polynomial_norm_squared': str(normalization),
            'initial_infinite_vacuum_tail_probability_upper': str(tail),
            'initial_trace_distance_upper': str(preparation_error),
            'first_segment_vector_errors': {k: str(v) for k, v in first_errors.items()},
            'maximum_expanded_total_occupation': maximum,
            'configurations': configurations,
            'infinite_feedback_contrast_interval': contrast,
            'nonzero_response_certified': F(contrast[1]) < 0 or F(contrast[0]) > 0,
            'scope': {'action': 'same declared interacting two-site Hamiltonian, infinite onsite Fock spaces',
                      'initial_state': 'infinite free coupled vacuum and blank gauge code; native preparation current zero',
                      'readout': 'actual half-line coordinate projectors, not sign of a compressed coordinate matrix',
                      'terminal_states': 'half-line polynomial wavefunctions with all first records and both final outcomes retained',
                      'errors': 'initial projection, full-H Taylor remainder, support discard and exact dyadic rounding included',
                      'spatial_limit_or_physical_services': False,
                      'independent_empirical_validation': False}}


def run(retain_polynomials=False):
    """Fixed exact certificate, refinement control and disconnected control."""
    from .halfline import inverse_sqrt_two_pi_interval
    from .occupation import certificate_table
    fine = certificate(retain_polynomials=retain_polynomials)
    coarse = certificate(cutoff=12, preparation_cutoff=8, bits=96)
    disconnected = certificate(cutoff=12, preparation_cutoff=8, bits=96, coupling=F(0))
    interval = lambda result: tuple(map(F, result['infinite_feedback_contrast_interval']))
    f, c0, zero = map(interval, (fine, coarse, disconnected))
    if not (f[1] < 0 and c0[1] < 0 and f[1]-f[0] < c0[1]-c0[0]
            and max(f[0], c0[0]) <= min(f[1], c0[1])
            and f[1]-f[0] < F(1, 100000) and zero[0] <= 0 <= zero[1]):
        raise ArithmeticError('fixed infinite signal/refinement/control certificate failed')
    for configuration in disconnected['configurations']:
        low, high = map(F, configuration['infinite_positive_probability_interval'])
        if not low <= F(1, 2) <= high:
            raise ArithmeticError('disconnected even scalar vacuum must have half-line weight 1/2')
    low, high = inverse_sqrt_two_pi_interval()
    return {'schema': 'oph.common_source.infinite_occupation_audit.v1',
            'fine': fine, 'coarse_control': coarse, 'disconnected_control': disconnected,
            'rectangular_cutoff_comparisons': certificate_table(),
            'two_level_readout_effect_error_interval': _interval(F(1, 2)-high, F(1, 2)-low),
            'all_pass': True,
            'conclusion': 'nonzero declared interacting feedback survives infinite occupation and the actual coordinate-sign readout',
            'physical_clock_calibration_or_independent_data': False}


def canonical(value):
    """Encode the exact receipt tree without JSON's implicit type coercions."""
    def validate(node):
        if type(node) is dict:
            if any(type(key) is not str for key in node):
                raise ValueError('receipt object keys must be strings')
            for item in node.values():
                validate(item)
        elif type(node) is list:
            for item in node:
                validate(item)
        elif type(node) not in (str, int, bool, type(None)):
            raise ValueError('exact receipt values must be JSON objects, lists, strings, integers, Booleans or null')
    validate(value)
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def compare_receipt(expected, actual):
    """Every exact field is binding, including types, errors and state hashes."""
    try:
        identical = canonical(expected) == canonical(actual)
    except (TypeError, ValueError) as error:
        raise ValueError('an exact JSON receipt required') from error
    if not identical:
        raise ValueError('infinite-occupation receipt differs from exact recomputation')


def main():
    import argparse
    from pathlib import Path
    from tools.strict_json import load
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path, help='create exact receipt; refuse overwrite')
    group.add_argument('--check', type=Path, help='compare every field by exact recomputation')
    parser.add_argument('--polynomials', action='store_true', help='include regenerable terminal polynomial coefficients')
    args = parser.parse_args()
    expected = load(args.check) if args.check else None
    result = run(retain_polynomials=args.polynomials)
    if args.check:
        compare_receipt(expected, result)
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    print(canonical({'all_pass': result['all_pass'],
                     'contrast_interval': result['fine']['infinite_feedback_contrast_interval'],
                     'physical_clock_calibration_or_independent_data': False}))


if __name__ == '__main__':
    main()
