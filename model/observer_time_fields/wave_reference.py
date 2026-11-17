"""A bounded native-return protocol for the coherent wave modular reference.

Positive random-length closure reconstructs exp[-arccos(K)]/Tr, without
preparing that coherent density on the working source. Length sampling and
full endpoint reads are declared services; all failures and caps are outputs.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from types import MappingProxyType

import sympy as sp

from . import spatial_dynamics as dynamics
from tools import strict_json

ROOT = dynamics.ROOT
SOURCE_PATHS = tuple(sorted(set((*dynamics.SOURCE_PATHS,
    'model/observer_time_fields/wave_reference.py'))))


def _cap(value, maximum=512):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError(f'an original length cap from zero through {maximum} required')
    return value


def _atan_inverse(integer, terms=48):
    x = F(1, integer)
    total = sum(((-1)**n*x**(2*n+1)/(2*n+1) for n in range(terms)), F())
    following = (-1)**terms*x**(2*terms+1)/(2*terms+1)
    return min(total, total+following), max(total, total+following)


@lru_cache(maxsize=1)
def _normalization():
    """Rational enclosure of c0=exp(-pi/2), using Machin and positive series."""
    a, b = _atan_inverse(5), _atan_inverse(239)
    pi = dynamics._outward((16*a[0]-4*b[1], 16*a[1]-4*b[0]), 1 << 96)
    lo, hi = pi[0]/2, pi[1]/2
    degree = 64
    low_exp = sum((lo**n/math.factorial(n) for n in range(degree+1)), F())
    high_exp = sum((hi**n/math.factorial(n) for n in range(degree+1)), F())
    high_exp += hi**(degree+1)/math.factorial(degree+1)/(1-hi/(degree+2))
    lower, upper = dynamics._outward((1/high_exp, 1/low_exp), 1 << 96)
    grid = 1 << 64
    implemented = F((lower*grid).numerator//(lower*grid).denominator, grid)
    if not F(1, 5) < implemented <= lower < upper < F(21, 100):
        raise ArithmeticError('the finite sampler normalization enclosure changed')
    return tuple((value.numerator, value.denominator) for value in (*pi, lower, upper, implemented))


def normalization():
    """Detached exact enclosures; callers cannot mutate the cached certificate."""
    a, b, lower, upper, implemented = (F(*pair) for pair in _normalization())
    return MappingProxyType({'pi_interval': (a, b), 'c0_interval': (lower, upper),
                            'implemented_c0': implemented})


def coefficients(cap):
    """Rational coefficients b_n of exp(arcsin(z)); no target spectrum input."""
    cap = _cap(cap)
    values = [F(1), F(1)]
    for n in range(cap-1):
        values.append(F(n*n+1, (n+2)*(n+1))*values[n])
    return tuple(values[:cap+1])


def sampler(cap):
    """Exact finite categorical law, with one explicit unexecuted-length cap."""
    values = coefficients(cap)
    c0 = normalization()['implemented_c0']
    weights = tuple(c0*value for value in values)
    cap_mass = 1-sum(weights, F())
    if any(value <= 0 for value in weights) or not 0 < cap_mass < 1:
        raise ArithmeticError('bounded length sampler lost normalization or positivity')
    common = math.lcm(*(value.denominator for value in (*weights, cap_mass)))
    counts = tuple(value.numerator*(common//value.denominator) for value in (*weights, cap_mass))
    if sum(counts) != common:
        raise ArithmeticError('finite categorical randomizer lost its complete partition')
    return MappingProxyType({'weights': weights, 'cap': cap_mass,
        'uniform_denominator': common, 'uniform_integer_counts': counts,
        'expected_native_attempts': sum((n*weight for n, weight in enumerate(weights)), F())})


def _negative_log_upper(value, terms=96):
    """Exact outward log bound via the positive atanh series."""
    value = F(value)
    if not 0 < value <= 1:
        raise ValueError('a positive probability required')
    x = (1-value)/(1+value)
    partial = 2*sum((x**(2*n+1)/(2*n+1) for n in range(terms)), F())
    return partial+2*x**(2*terms+1)/((2*terms+1)*(1-x*x))


def cap_bounds(cap):
    cap = _cap(cap)
    mass = sum(coefficients(cap), F())
    lo, hi = normalization()['c0_interval']
    lower_tail, upper_tail = 1-hi*mass, 1-lo*mass
    implemented = sampler(cap)
    if not 0 < lower_tail <= upper_tail <= implemented['cap'] < 1:
        raise ArithmeticError('analytic and implemented length-tail bounds disagree')
    log_tail = dynamics._outward((upper_tail, upper_tail), 1 << 64)[1]
    return {'cap': cap,
        'ideal_length_tail_interval': dynamics._outward((lower_tail, upper_tail)),
        'complete_implemented_cap_probability': dynamics._outward((implemented['cap'], implemented['cap']))[1],
        'conditional_unforced_return_history_TV_upper': dynamics._outward((upper_tail, upper_tail))[1],
        'reconstructed_reference_trace_distance_upper': dynamics._outward((upper_tail, upper_tail))[1],
        'wave_modular_generator_norm_error_upper': dynamics._outward((F(0), _negative_log_upper(1-log_tail)))[1],
        'uniform_return_probability_lower': normalization()['implemented_c0'],
        'expected_executed_native_attempts_interval': dynamics._outward((implemented['expected_native_attempts'], implemented['expected_native_attempts'])),
        'maximum_executed_native_attempts': cap,
        'uniform_categorical_denominator_bits': implemented['uniform_denominator'].bit_length(),
        'ideal_uncapped_expected_length_is_infinite': True,
        'ideal_unforced_return_conditioned_expected_length_is_infinite': True}


def refinement_budget(level):
    """Sufficient joint wave-scale/cap theorem, not an executed long sampler.

    Positive PGF coefficients give delta_h<=4/sqrt(h+1). With s=4^level
    and h=64s^4, delta<=1/(2s^2), logarithm error<=1/s^2, and multiplying
    by sqrt(630)s bounds the scaled wave-generator error by26/s.
    """
    level = dynamics._level(level)
    scale = 4**level
    return {'level': level, 'sufficient_native_length_cap': 64*scale**4,
        'reference_trace_distance_upper': F(1, 2*scale**2),
        'unscaled_modular_generator_error_upper': F(1, scale**2),
        'scaled_wave_generator_error_upper': F(26, scale),
        'cap_sampler_at_this_budget_executed': False,
        'bound_is_sufficient_not_minimal': True}


@lru_cache(maxsize=4)
def _powers(cap):
    counts = dynamics.native_density()['counts']
    power = sp.eye(66)
    result = []
    for n in range(cap+1):
        result.append(sp.ImmutableMatrix(power))
        power = counts*power
    return tuple(result)


def complete_law(cap=2):
    """All initial/final configurations and sampled lengths, plus the cap.

    Matrix entries are joint probabilities, not quantum output-state entries.
    Both unequal endpoints and equal endpoints are kept without conditioning.
    """
    cap = _cap(cap, 8)
    weights = sampler(cap)
    tables = tuple(sp.ImmutableMatrix(sp.Rational(weights['weights'][n])*power/(66*60**n))
                   for n, power in enumerate(_powers(cap)))
    mass = sum((sum(table) for table in tables), sp.Rational(weights['cap']))
    if mass != 1 or any(entry < 0 for table in tables for entry in table):
        raise ArithmeticError('the complete length/endpoint law is not normalized')
    return MappingProxyType({'endpoint_tables': tables, 'cap_probability': weights['cap']})


def endpoint_oracle(cap=2):
    """Direct endpoint/coin recursion; no kernel, resolvent or matrix powers."""
    cap = _cap(cap, 3)
    states = tuple(dynamics.native.states(binary=True))
    index = {state: i for i, state in enumerate(states)}
    seams, _ = dynamics.source_green.passage.recovered_carrier()
    counts = [sp.zeros(66) for _ in range(cap+1)]
    for initial, before in enumerate(states):
        frontier = {before: 1}
        for n in range(cap+1):
            for state, multiplicity in frontier.items():
                counts[n][index[state], initial] += multiplicity
            following = defaultdict(int)
            for state, multiplicity in frontier.items():
                for edge in seams:
                    for coin in (0, 1):
                        after, _ = dynamics.native.native(state, edge, coin)
                        following[after] += multiplicity
            frontier = following
    weights = sampler(cap)
    return MappingProxyType({'endpoint_tables': tuple(sp.ImmutableMatrix(
        sp.Rational(weights['weights'][n])*table/(66*60**n)) for n, table in enumerate(counts)),
        'cap_probability': weights['cap']})


def reference(cap=4):
    """Reconstructed positive density; not the source's postselected density."""
    cap = _cap(cap, 8)
    b = coefficients(cap)
    polynomial = sum((sp.Rational(value)*power/60**n
        for n, (value, power) in enumerate(zip(b, _powers(cap), strict=True))), sp.zeros(66))
    trace = sp.trace(polynomial)
    density = polynomial/trace
    transfer = polynomial/sp.Rational(sum(b, F()))
    law = complete_law(cap)
    combined = sum(law['endpoint_tables'], sp.zeros(66))
    return_mass = sp.trace(combined)
    if combined/return_mass != density:
        raise ArithmeticError('retained endpoint law does not reconstruct the reference')
    native_return_state = sp.diag(*(combined[i, i]/return_mass for i in range(66)))
    if cap and native_return_state == density:
        raise ArithmeticError('coherent reconstruction was confused with native closure preparation')
    if transfer*sp.ones(66, 1) != sp.ones(66, 1):
        raise ArithmeticError('finite wave-reference transfer lost its ground eigenvalue')
    return MappingProxyType({'density': sp.ImmutableMatrix(density),
        'transfer': sp.ImmutableMatrix(transfer), 'return_probability': dynamics._fraction(return_mass),
        'different_probability': dynamics._fraction(sum(combined)-return_mass),
        'native_return_state': sp.ImmutableMatrix(native_return_state),
        'cap_probability': law['cap_probability']})


def run():
    normalization()
    rows = []
    for cap in (0, 1, 2, 4, 8):
        value = reference(cap)
        law = complete_law(cap)
        if cap <= 2 and law != endpoint_oracle(cap):
            raise ArithmeticError('native length/endpoint oracle disagrees')
        rows.append({'cap': cap, 'reconstructed_density_sha256': dynamics._digest(value['density']),
            'complete_endpoint_law_sha256': tuple(dynamics._digest(x) for x in law['endpoint_tables']),
            'return_probability': value['return_probability'],
            'different_probability': value['different_probability'],
            'cap_probability': value['cap_probability'],
            'native_and_reconstructed_density_matrices_equal': cap == 0})
    return dynamics._json({'schema': 'oph.observer_time_fields.wave_reference.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'length_law': {'PGF': 'exp(arcsin(z)-pi/2)', 'b0': 1, 'b1': 1,
            'positive_rational_recurrence': 'b[n+2]=(n*n+1)*b[n]/((n+2)*(n+1))',
            'implemented_c0': normalization()['implemented_c0'],
            'infinite_ideal_mean_length': True},
        'complete_native_protocols': rows,
        'all_source_cap_bounds': [cap_bounds(cap) for cap in (0, 2, 8, 16, 64, 256, 512)],
        'joint_wave_refinement_budgets': [refinement_budget(level) for level in (0,1,2,4,10)],
        'modular_reference': {'density': 'exp(-arccos K)/Tr exp(-arccos K)',
            'positive_generator': 'Theta=arccos K',
            'finite_cap_generator': 'Theta_h=-log(P_h(K)/P_h(1))',
            'finite_cap_operator_bound': '0<=Theta-Theta_h<=-log(1-delta_h) I',
            'dual_modular_unitary': 'exp(+i u Theta); coherent wave phases use u=n',
            'executed_arc_walk_compression': 'A*W^nA=cos(n Theta), containing both signs',
            'scaled_coherent_wave_parameter': 't=n/sqrt(M*4^level)',
            'finite_cap_unitary_operator_error': 'min(2,abs(u)*[-log(1-delta_h)])'},
        'scope': {'same_positive_native_configuration_transfer': True,
            'independent_finite_categorical_length_sampler_admitted': True,
            'complete_initial_final_configuration_reads_admitted': True,
            'every_return_failure_and_length_cap_retained': True,
            'no_renormalization_of_the_complete_native_output': True,
            'conditional_reference_bound_is_for_unforced_return': True,
            'coherent_record_and_reflection_phase_services_required': True,
            'arc_walk_executes_vertex_modular_unitary_or_all_modular_interventions': False,
            'variable_length_reference_forced_by_committed_repair_law': False,
            'reconstructed_density_physically_prepared_on_native_working_source': False,
            'existing_three_class_probe_equals_density_wave_generator': False,
            'physical_QFT_vacuum_CFT_or_Lorentz_time_selected': False}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        from model.observer.scalar_response.experiment import compare_receipt
        compare_receipt(strict_json.load(args.check), result)
        print('Native random-length closure, wave modular reference and complete cap bounds replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, indent=2, sort_keys=True); stream.write('\n')
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
