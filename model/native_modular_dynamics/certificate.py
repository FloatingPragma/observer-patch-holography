"""Exact native-word polynomials for modular correlations and sector references."""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from math import comb
from pathlib import Path
import json

import numpy as np

from tools import strict_json
from ..native_quantum_reconstruction import certificate as source
from ..native_thermal_return import certificate as observer

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
D = source.DIMENSION
ZERO = (Q(), Q())
ONE = (Q(1), Q())
require = source.require


def rational(value):
    require(type(value) in (int, Q), 'exact rational required; no float or Boolean')
    return Q(value)


def integer(value, lo, hi, name):
    require(type(value) is int and lo <= value <= hi, f'{name}: integer in [{lo},{hi}] required')
    return value


def add(z, w):return (z[0]+w[0], z[1]+w[1])
def mul(z, w):return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def scale(a, z):return (a*z[0], a*z[1])
def conjugate(z):return (z[0], -z[1])
def magnitude_upper(z):return abs(z[0])+abs(z[1])


def x_coefficients(t, degree):
    t = rational(t); integer(degree, 0, 64, 'polynomial degree')
    out = [ONE]
    for k in range(degree):
        out.append(scale(Q(1, k+1), mul(out[-1], (Q(k), -t))))
    return tuple(out)


def coefficients(t, degree):
    """Expand sum c_k(t)(I-P)^k into exact Gaussian-rational P powers."""
    integer(degree, 0, 64, 'polynomial degree')
    out = [ZERO]*(degree+1)
    for k, c in enumerate(x_coefficients(t, degree)):
        for j in range(k+1):
            out[j] = add(out[j], scale((-1)**j*comb(k, j), c))
    return tuple(out)


def tail_bound(t, degree):
    t = abs(rational(t)); integer(degree, 0, 64, 'polynomial degree')
    if not t:return Q()
    n = (t.numerator+t.denominator-1)//t.denominator
    require(n <= 8, 'certified time window requires absolute time at most eight')
    r = Q(1, 3)
    value = (1-r)**(-n)-sum((Q(comb(n+k-1, k))*r**k for k in range(degree+1)), Q())
    require(value > 0, 'binomial tail bound is not positive')
    return value


@lru_cache(maxsize=130)
def power_counts(k):
    """Exact integer native branch counts, including the zero-attempt identity."""
    integer(k, 0, 128, 'native power')
    if k == 0:return tuple(tuple(int(i == j) for j in range(D)) for i in range(D))
    counts, _ = source.source_counts()
    result = np.array(power_counts(k-1), dtype=object) @ np.array(counts, dtype=object)
    return tuple(tuple(int(v) for v in row) for row in result)


def matrix_entry(k, a, b):
    source.index(a); source.index(b)
    _, denominator = source.source_counts()
    return Q(power_counts(k)[a][b], denominator**k)


def normalization(steps):
    integer(steps, 1, 16, 'reference loop length')
    return sum(power_counts(steps)[i][i] for i in range(D))


def reference_entry(steps, power, a, b):
    return Q(power_counts(steps+power)[a][b], normalization(steps)*60**power)


def word_moment(steps, indices, powers):
    require(type(indices) is tuple and len(indices) >= 1 and type(powers) is tuple
            and len(powers) == len(indices), 'matching nonempty coordinate and power tuples required')
    for i in indices:source.index(i)
    for k in powers:integer(k, 0, 64, 'word power')
    value = reference_entry(steps, powers[0], indices[-1], indices[0])
    for j in range(1, len(indices)):
        value *= matrix_entry(powers[j], indices[j-1], indices[j])
    return value


def relative_times(times):
    require(type(times) is tuple and 1 <= len(times) <= 4, 'one through four exact times required')
    times = tuple(rational(t) for t in times)
    return (times[0]-times[-1],)+tuple(times[j]-times[j-1] for j in range(1, len(times)))


def correlation(steps, indices, times, degree=16):
    """Finite linear combination of native closed-word moments."""
    integer(steps, 1, 16, 'reference loop length')
    require(type(indices) is tuple and len(indices) == len(times), 'one coordinate for each time required')
    for i in indices:source.index(i)
    result = ONE
    for j, t in enumerate(relative_times(times)):
        value = ZERO
        for k, c in enumerate(coefficients(t, degree)):
            entry = (reference_entry(steps, k, indices[-1], indices[0]) if j == 0
                     else matrix_entry(k, indices[j-1], indices[j]))
            value = add(value, scale(entry, c))
        result = mul(result, value)
    return result


def approximation_bound(times, degree):
    value = Q(1)
    for t in relative_times(times):value *= 1+tail_bound(t, degree)
    return value-1


def record_error_bound(times, degree, delta):
    """Numerator-probability and independent-return-normalizer errors."""
    delta = rational(delta)
    require(0 <= delta < Q(1, D), 'record error must be below the return lower bound')
    weight, norm = Q(1), Q(1)
    for t in relative_times(times):
        weight *= sum(magnitude_upper(c) for c in coefficients(t, degree))
        norm *= 1+tail_bound(t, degree)
    return delta*(weight+norm)/(Q(1, D)-delta)


def entry_error_bound(steps, times, degree, entry_error):
    """Requires Hermitian contraction P_hat; not an empirical PSD assertion."""
    integer(steps, 1, 16, 'reference loop length')
    e = rational(entry_error)
    require(e >= 0, 'entry error must be nonnegative')
    eta = D*e
    require(D*steps*eta < 1, 'matrix reference normalization error is too large')
    true_norm, estimated_norm = Q(1), Q(1)
    for t in relative_times(times):
        n = 1+tail_bound(t, degree)
        derivative = sum(k*magnitude_upper(c) for k, c in enumerate(coefficients(t, degree)))
        true_norm *= n
        estimated_norm *= n+eta*derivative
    state_error = 2*D*steps*eta/(1-D*steps*eta)
    return state_error*estimated_norm+(estimated_norm-true_norm)


def rounded_upper(value, denominator=10**42):
    value = rational(value)
    return Q((value.numerator*denominator+value.denominator-1)//value.denominator, denominator)


def observer_budget(steps, times, degree):
    body = steps+len(times)*degree
    budget = observer.bounds(clock='attempt', length=body, seed_scans=16, retries=128,
                             warmup=2048, absorb=True)
    delta = rounded_upper(budget['complete_law_error_upper'])
    return {'maximum_body_attempts': body, 'total_attempts_max': budget['native_attempts_max'],
            'complete_event_probability_error_upper': delta,
            'polynomial_bias_upper': approximation_bound(times, degree),
            'record_amplification_error_upper': rounded_upper(record_error_bound(times, degree, delta)),
            'finite_sampling_error_included': False}


def corners(steps=2):
    """Exact common-reference conditioning and its modular compatibility test."""
    density = tuple(tuple(Q(v, normalization(steps)) for v in row) for row in power_counts(steps))
    cells = tuple(tuple(i for i, state in enumerate(source.STATES) if state[0] == bit) for bit in (0, 1))
    weights = tuple(sum(density[i][i] for i in cell) for cell in cells)
    require(sum(weights) == 1 and all(w > 0 for w in weights), 'sector references lost faithfulness')
    counts, denominator = source.source_counts()
    cross = next((i, j) for i in cells[0] for j in cells[1] if counts[i][j])
    require(density[cross[0]][cross[1]] > 0, 'proper-corner nonintertwining witness lost')
    within = next((i, j) for i in cells[0] for j in cells[0] if i != j and counts[i][j])
    require(density[within[0]][within[1]] > 0, 'completed sector action became tracial')
    nested = tuple(i for i in cells[0] if source.STATES[i][1] == 0)
    direct_weight = sum(density[i][i] for i in nested)
    conditional_weight = sum(density[i][i]/weights[0] for i in nested)
    require(weights[0]*conditional_weight == direct_weight, 'nested conditioning probability failed')
    require(all((density[i][j]/weights[0])/conditional_weight == density[i][j]/direct_weight
                for i in nested for j in nested), 'nested reference restrictions do not compose')
    a, _ = cross
    leak = sum((Q(counts[a][j]*counts[j][a], denominator**2) for j in cells[1]), Q())
    require(leak > 0, 'compression falsely commutes with transfer powers')
    # All cross-sector entries vanish in the completed reference, while each
    # normalized block is exactly the corresponding conditional reference.
    completed = tuple(tuple(density[i][j] if source.STATES[i][0] == source.STATES[j][0] else Q()
                            for j in range(D)) for i in range(D))
    require(all(completed[i][j] == weights[k]*(density[i][j]/weights[k])
                for k, cell in enumerate(cells) for i in cell for j in cell), 'block reference identity failed')
    return {'sector_sizes': tuple(map(len, cells)), 'sector_weights': weights,
            'nested_sector_size': len(nested), 'nested_weight': direct_weight,
            'cross_sector_indices': cross, 'cross_sector_reference_entry': density[cross[0]][cross[1]],
            'within_sector_indices': within, 'within_sector_reference_entry': density[within[0]][within[1]],
            'two_step_compression_leakage_diagonal': leak,
            'completed_reference_differs_from_original': completed != density,
            'completed_sector_modular_intertwining': True,
            'original_coordinate_corner_modular_intertwining': False,
            'general_tensor_region_no_go_claimed': False}


def numerical_control(steps, indices, times, degree, diagnostics=None):
    """Independent floating spectral evaluation, diagnostic not proof input."""
    counts, denominator = source.source_counts()
    P = np.array(counts, dtype=float)/denominator
    eigenvalues, vectors = np.linalg.eigh(P)
    rho = (vectors*(eigenvalues**steps/sum(eigenvalues**steps))) @ vectors.T
    result = 1+0j
    for j, t in enumerate(relative_times(times)):
        U = (vectors*np.exp(1j*float(t)*np.log(eigenvalues))) @ vectors.T
        result *= (rho@U)[indices[-1], indices[0]] if j == 0 else U[indices[j-1], indices[j]]
    exact = correlation(steps, indices, times, degree)
    residual = abs(complex(float(exact[0]), float(exact[1]))-result)
    bound = float(approximation_bound(times, degree))
    require(residual <= bound+1e-13, 'independent modular spectral diagnostic disagrees')
    if diagnostics is not None:
        diagnostics.append(format(residual, '.8g'))
    return {'within_analytic_bound_plus_roundoff': True,
            'analytic_bias_bound': str(approximation_bound(times, degree)),
            'floating_evaluation_is_proof_premise': False}


def encode(value):
    if isinstance(value, Q):return str(value)
    if isinstance(value, dict):return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):return [encode(v) for v in value]
    return value


def inherited_pins():
    checked = {}
    for receipt in ('model/native_quantum_reconstruction/results.json', 'model/native_thermal_return/results.json'):
        stored = strict_json.load(ROOT/receipt)
        for path, expected in stored['source_sha256'].items():
            require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected,
                    f'inherited source binding changed: {path}')
        checked[receipt] = len(stored['source_sha256'])
    return checked


def report(diagnostics=None):
    pins = inherited_pins()
    counts, _ = source.source_counts()
    a, b = next((i, j) for i in range(D) for j in range(D) if i != j and counts[i][j])
    c = next(j for j in range(D) if j not in (a, b) and counts[b][j])
    cases = []
    for indices, times in (((a, b), (Q(), Q(1, 2))), ((a, b, c), (Q(-1, 2), Q(), Q(1, 2)))):
        value = correlation(2, indices, times, 20)
        require(value == correlation(2, indices, tuple(t+7 for t in times), 20),
                'stationary common time shift changed polynomial histories')
        cases.append({'indices': indices, 'times_in_common_generator_units': times,
                      'exact_correlation': value, 'budget': observer_budget(2, times, 20),
                      'spectral_diagnostic': numerical_control(2, indices, times, 20, diagnostics)})
    require(coefficients(2*Q(1, 4), 20) == coefficients(4*Q(1, 8), 20), 'cross-reference clock scaling failed')
    require(reference_entry(1, 0, a, a) != reference_entry(2, 0, a, a)
            or reference_entry(1, 0, a, b) != reference_entry(2, 0, a, b),
            'distinct loop references falsely identified')
    from . import interaction
    interaction_result = interaction.certificate()
    paths = ['model/native_modular_dynamics/certificate.py', 'model/native_modular_dynamics/interaction.py',
             'model/native_effective_fields/construction.py', 'model/native_modular_dynamics/PROOF.md',
             'model/native_modular_dynamics/test_certificate.py',
             'model/native_quantum_reconstruction/certificate.py', 'model/native_quantum_reconstruction/results.json',
             'model/native_thermal_return/certificate.py', 'model/native_thermal_return/results.json']
    return encode({'schema': 'oph.native_modular_dynamics.v1', 'all_pass': True,
        'inherited_source_pins': pins, 'native_branch_counts': D*60,
        'source_spectral_bounds': source.source_spectral_bounds(),
        'ordered_modular_histories': cases,
        'same_generator_clock_relation': 't=L*s for U(t)=P^(it); Hamiltonian Heisenberg time tau=-L*s',
        'sector_references': corners(2),
        'native_interaction_diagnosis': interaction_result,
        'matrix_stability_example': {'entry_error': '1/100000000000000000000',
            'correlation_error_upper': rounded_upper(entry_error_bound(2, (Q(), Q(1, 2)), 12, Q(1, 10**20))),
            'Hermitian_contraction_hypothesis_required': True},
        'scope': {'native_weighted_record_statistics_constructed': True,
                  'modular_flow_executed_as_native_quantum_channel': False,
                  'time_shift_stationarity_in_reconstructed_reference': True,
                  'same_state_as_unconditioned_native_working_carrier': False,
                  'coordinate_corner_obstruction_is_all_regional_algebras_no_go': False,
                  'physical_clock_or_energy_units_identified': False},
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true');parser.add_argument('--check', action='store_true')
    args = parser.parse_args(); diagnostics = []; result = report(diagnostics); path = HERE/'results.json'
    if args.write:path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check:require(strict_json.load(path) == result, 'native modular dynamics receipt changed')
    print(json.dumps({'all_pass': True, 'modular_history_cases': len(result['ordered_modular_histories']),
                      'coordinate_sectors': result['sector_references']['sector_sizes'],
                      'spectral_absolute_differences_diagnostic_only': diagnostics}))


if __name__ == '__main__':main()
