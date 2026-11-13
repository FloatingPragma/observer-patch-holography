"""Full-source correlation enclosures with explicit representation-tail errors."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import numpy as np
from model.openwave_coupled_haar_spectrum.certificate import calculate, display, ROUND_DENOMINATOR
from model.openwave_su2_spectral.spectrum import sqrt_upper, exp_bounds


def dyadic_outer(interval, bits=80):
    """Outward rounding is an exact rational inequality, never a tolerance."""
    lo, hi = map(F, interval)
    assert lo <= hi
    denominator = 2**bits
    lower = F((lo*denominator).__floor__(), denominator)
    upper = F((hi*denominator).__ceil__(), denominator)
    assert lower <= lo <= hi <= upper
    return lower, upper


def strings(interval):
    return [str(x) for x in dyadic_outer(interval)]


def upper_text(value):
    return str(dyadic_outer((value, value))[1])


@lru_cache(None)
def exponential_decay_interval(x):
    x = F(x)
    if x > 64:
        # Taylor's rational lower bound proves exp(64)>2**92.
        assert exp_bounds(F(64))[0] > 2**92
        return F(0), F(1, 2**92)
    lo, hi = dyadic_outer((x, x), bits=60)
    result = exp_bounds(-hi)[0], exp_bounds(-lo)[1]
    return dyadic_outer(result)


@lru_cache(None)
def spectral_data(key):
    record, data = calculate(key)
    packet, u, norm = data['packet'], data['u'], data['norm']
    centered = [F(int(x), ROUND_DENOMINATOR)-data['m1']*int(y) for x, y in zip(data['w1']@u, u)]
    centered_norm2 = sum(x*x for x in centered)/norm
    assert centered_norm2 < 2
    frame = packet['complete_finite_frame']
    D = frame['dyadic_denominator']
    O = np.array([[int(x) for x in row] for row in frame['dyadic_eigenframe']], dtype=object)
    weights = [sum(F(int(O[i, j]), D)*centered[i] for i in range(len(u)))**2/norm for j in range(len(u))]
    energies = [F(x)/D for x in frame['dyadic_eigenvalues']]
    assert min(energies) > 0 and min(weights) >= 0
    return record, data, centered, centered_norm2, weights, energies


def euclidean_record(key, time):
    time = F(time)
    if time < 0:
        raise ValueError('Euclidean time is nonnegative')
    record, data, centered, norm2, weights, energies = spectral_data(key)
    cap, beta = data['packet']['cap'], data['beta']
    mu, norm = data['mu'], data['norm']
    exact_lo = exact_hi = F(0)
    for weight, energy in zip(weights, energies):
        lo, hi = exponential_decay_interval(time*(energy-mu))
        exact_lo += weight*lo
        exact_hi += weight*hi

    frame_eta = F(record['complete_finite_spectrum']['exact_Gram_error_upper'])
    generator_error = F(record['complete_finite_spectrum']['all_eigenvalue_error_upper'])
    norm_bound = sqrt_upper(norm2)
    # Every path to the finite spin boundary changes a loop label by one per
    # magnetic insertion. The minimum is chosen among valid certified bounds.
    choices = []
    for inner in range(cap):
        tail2 = sum(centered[i]**2 for i, (a, b, _) in enumerate(data['basis']) if max(a, b) > inner)/norm
        tail = sqrt_upper(tail2)
        distance = cap-inner
        x = 2*beta*time
        reach = exponential_decay_interval(-x)[1]*x**distance/factorial(distance)
        finite_boundary_error = norm_bound*2*beta*time*(tail+reach*norm_bound)
        choices.append((finite_boundary_error, inner, tail, reach))
    boundary_error, inner, tail, reach = min(choices)
    outside_initial = data['boundary_norm']
    coefficient_error = data['error1']
    # The exact centered trial has norm<=1. Its represented counterpart differs
    # by at most the missing W component plus two coefficient-rounding errors.
    initial_vector_error = outside_initial+2*coefficient_error
    initial_quadratic_error = (2+initial_vector_error)*initial_vector_error
    ground_error = 12*data['distance']
    E0lo, E0hi = data['ground_interval']
    normalization_error = time*max(abs(E0lo-mu), abs(E0hi-mu))*norm2
    amplification = exponential_decay_interval(-time*max(E0hi, mu))[1]
    finite_matrix_error = norm2*(time*generator_error+3*frame_eta)
    error = ground_error+amplification*(boundary_error+initial_quadratic_error+finite_matrix_error+normalization_error)
    variance = tuple(map(F, record['Wilson_variance_interval_exact']))
    full = (max(F(0), exact_lo-error), min(F(1), variance[1], exact_hi+error))
    assert full[0] <= full[1]
    if time == 0:
        full = (max(full[0], variance[0]), min(full[1], variance[1]))
        assert full[0] > 0
    return {'Euclidean_time_exact': str(time), 'full_source_connected_correlation_interval_exact': strings(full),
            'full_source_connected_correlation_interval': display(full),
            'absolute_error_upper_exact': upper_text(error),
            'boundary_buffer_inner_twice_spin': inner,
            'initial_high_spin_centered_norm_upper_exact': upper_text(tail),
            'Dyson_boundary_reach_upper_exact': upper_text(reach),
            'full_ground_comparison_error_upper_exact': upper_text(ground_error),
            'positive_finite_spectral_sum_interval_exact': strings((exact_lo, exact_hi))}


def reflection_form(weights, energies, times, coefficients):
    """A sum of nonnegative squares is an exact positive reflection form."""
    if len(times) != len(coefficients) or any(F(t) < 0 for t in times):
        raise ValueError('paired coefficients and nonnegative reflection times required')
    if any(F(w) < 0 for w in weights) or any(F(e) < 0 for e in energies):
        raise ValueError('a positive spectral measure and generator are required')
    if len(weights) != len(energies):
        raise ValueError('the full finite spectral measure is paired')
    lower = upper = F(0)
    for weight, energy in zip(weights, energies):
        lo = hi = F(0)
        for time, coefficient in zip(times, coefficients):
            a, b = exponential_decay_interval(F(time)*F(energy))
            coefficient = F(coefficient)
            lo += coefficient*(a if coefficient >= 0 else b)
            hi += coefficient*(b if coefficient >= 0 else a)
        square = (F(0) if lo <= 0 <= hi else min(lo*lo, hi*hi), max(lo*lo, hi*hi))
        lower += weight*square[0]
        upper += weight*square[1]
    return lower, upper


def clock_scaling(record, spacings, physical_time=F(1, 4)):
    gap = tuple(map(F, record['full_Haar_gap_interval_exact']))
    variance = tuple(map(F, record['Wilson_variance_interval_exact']))
    assert gap[0] > 0 and variance[0] > 0
    rows = []
    for a in spacings:
        a = F(a)
        if a <= 0:
            raise ValueError('positive spatial spacing is required')
        decay = exponential_decay_interval(gap[0]*F(physical_time)/a)[1]
        rows.append({'spacing_exact': str(a), 'source_clock_changed': False,
                     'physical_clock_gap_interval_exact': strings((gap[0]/a, gap[1]/a)),
                     'fixed_positive_time_exact': str(physical_time),
                     'connected_correlation_upper_exact': upper_text(variance[1]*decay),
                     'zero_time_variance_lower_exact': str(dyadic_outer(variance)[0])})
    return rows


def run():
    keys = ['shared_g2_2_cap4', 'shared_g2_2_cap6', 'shared_g2_2_cap8']
    times = (F(0), F(1, 8), F(1, 4), F(1, 2))
    records = []
    for key in keys:
        record, data, centered, norm2, weights, energies = spectral_data(key)
        rows = [euclidean_record(key, t) for t in times]
        reflection = reflection_form(weights, energies, (0,F(1,8),F(1,4)), (1,-2,1))
        assert reflection[0] >= 0
        records.append({'proof_coordinate_key': key,
                        'source_Haar_representations_truncated': False,
                        'all_finite_Galerkin_spectral_weights_exact': [str(w) for w in weights],
                        'all_finite_Galerkin_energy_coordinates_exact': [str(e) for e in energies],
                        'vacuum_weight_retained': True,
                        'correlations': rows,
                        'finite_positive_reflection_form_interval_exact': strings(reflection)})
    for i in range(len(times)):
        intervals = [tuple(map(F, r['correlations'][i]['full_source_connected_correlation_interval_exact'])) for r in records]
        assert max(x[0] for x in intervals) <= min(x[1] for x in intervals)
    primary = calculate(keys[-1])[0]
    positive_time = records[-1]['correlations'][2]
    assert F(positive_time['full_source_connected_correlation_interval_exact'][0]) > 0
    trajectory = clock_scaling(primary, (1,F(1,2),F(1,4),F(1,8),F(1,16),F(1,32)))
    assert all(F(a['connected_correlation_upper_exact']) > F(b['connected_correlation_upper_exact']) for a,b in zip(trajectory,trajectory[1:]))
    return {'schema': 'oph.coupled_Haar_Euclidean_refinement.v1',
            'verdict': 'PASS_FULL_SOURCE_EUCLIDEAN_RECORDS_AND_CLOCK_OBSTRUCTION',
            'source': 'Same full-Haar coupled two-plaquette Hamiltonian, fixed g²2, declared source clock',
            'representation_refinement_records': records,
            'fixed_coupling_spacing_trajectory': trajectory,
            'proved': {'full_source_Euclidean_correlations_enclosed': True,
                       'nonzero_positive_time_Wilson_response': True,
                       'temporal_reflection_positivity_of_same_Hamiltonian': True,
                       'finite_nonvacuum_centered_energy_witness': True,
                       'fixed_coupling_source_clock_gap_diverges_as_inverse_spacing': True,
                       'fixed_positive_time_response_collapses_in_that_sequence': True},
            'continuum_acceptance': {'fixed_graph_fixed_coupling_trajectory_admitted': False,
                'reason': 'A positive zero-time variance and a source-clock gap proportional to1/a give a discontinuous limiting correlator; this is not a nontrivial strongly continuous field limit',
                'required_physical_gap_scaling': 'Along a volume-preserving changing-coupling source family, the dimensionless physical gap must scale with a, together with uniform local Euclidean field limits',
                'physical_clock_rescaling_is_a_different_attachment': True,
                'missing_source_results': ['volume-preserving spatial refinement', 'four-dimensional Euclidean invariance and regularity',
                    'noncollapsed common local field limits', 'pure Yang-Mills short-distance identification',
                    'uniform physical spectral-gap transfer']},
            'scope': {'imaginary_time_filter_is_real_time_observer': False,
                      'positive_gap_exact_conformal_vacuum_theory': False,
                      'four_dimensional_interacting_continuum_constructed': False,
                      'Clay_problem_solved': False}}
