"""Execute certified powers of the ideal source word and infer eight weights.

This is an accelerated, error-certified simulation of the specified long
native instrument. It does not individually execute billions of raw pulses
or evaluate their separate first-failure probabilities.
"""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from ..certified_reference import interval as iv
from . import factorization as f, certified_power as cp, budgets, inversion


TRIALS = 9*2**29
SEED = 20261102
SYSTEMATIC = F(1, 65536)
SAMPLING = F(1, 32768)


def expanded(center, error):
    return iv.Interval(iv.Interval.point(center-error).lo, iv.Interval.point(center+error).hi)


@lru_cache(maxsize=1)
def execution():
    base = cp.from_intervals(*f.step_intervals())
    bound = budgets.certificate()
    q = f.probe_vector()
    random = np.random.Generator(np.random.PCG64(SEED))
    settings = []; powers = []; zero_counts = []; quadrature_counts = []
    for index, time in enumerate(budgets.TIMES):
        power = cp.power(base, index*budgets.STEPS_PER_INTERVAL)
        response = cp.expectation(power, q)
        real = expanded(F(response['real_center']), power.error)
        imag = expanded(F(response['imag_center']), power.error)
        # H_N6=T(H_N2+2), whereas R has H eigenvalue3. Hence +time.
        phase_real, phase_imag = inversion.trig_interval(F(-time), F(-time))
        cr = F(1, 16)+F(15, 16)*(phase_real*real-phase_imag*imag)
        ci = F(15, 16)*(phase_real*imag+phase_imag*real)
        probabilities = ((1+cr)/2, (1-ci)/2)
        rows = []
        for phase, probability in enumerate(probabilities):
            midpoint = F(probability.lo+probability.hi, 2*iv.SCALE)
            center = min(F(1), max(F(), midpoint))
            radius = max(abs(center-F(probability.lo, iv.SCALE)), abs(center-F(probability.hi, iv.SCALE)))
            float_probability = float(center)
            binary_rounding = abs(F(float_probability)-center)
            native = F(bound['schedule'][index]['native_systematic_probability_bounds'][phase])
            systematic = native+radius+binary_rounding
            if systematic >= SYSTEMATIC: raise ArithmeticError('long-time source/power budget exceeds the inversion premise')
            count = int(random.binomial(TRIALS, float_probability))
            if abs(F(count, TRIALS)-F(float_probability)) > SAMPLING:
                raise ArithmeticError('seeded simulated counts miss the declared simultaneous sampling event')
            (zero_counts if phase == 0 else quadrature_counts).append(count)
            rows.append({'phase_over_pi': '0' if phase == 0 else '1/2',
                         'ideal_source_word_probability_interval': probability.json(),
                         'sampling_probability_binary64': float_probability,
                         'sampling_probability_exact_binary_rational': str(F(float_probability)),
                         'certified_power_and_probability_rounding_error': str(radius+binary_rounding),
                         'total_systematic_probability_bound': str(systematic),
                         'trials': TRIALS, 'pointer_zero_count': count,
                         'all_other_complete_outcomes_count': TRIALS-count,
                         'simulated_frequency': str(F(count, TRIALS))})
        settings.append({'time': time, 'exact_word_power': power.certificate(),
                         'parent_response': response,
                         'reference_relative_response_interval': {'real': cr.json(), 'imag': ci.json()},
                         'pointer_settings': rows})
        powers.append({'time': time, 'precision_bits': power.bits,
                       'real_integer_numerators': power.real.tolist(),
                       'imag_integer_numerators': power.imag.tolist()})
    reconstructed = inversion.reconstruct(zero_counts, quadrature_counts, TRIALS, SYSTEMATIC, SAMPLING)
    if not reconstructed['all_eight_weights_strictly_positive']:
        raise ArithmeticError('the simulated complete instrument did not resolve all eight supplied lines')
    comparisons = []
    # This independent comparison is deliberately AFTER inference. The
    # known source residues are never supplied to the inverse or counts.
    for observed, source in zip(reconstructed['weights'], inversion.source_rows(), strict=True):
        low, high = map(F, observed['weight_interval'])
        exact_low, exact_high = map(F, source['positive_normalized_probe_weight_interval'])
        if not low < exact_low <= exact_high < high:
            raise ArithmeticError('inferred weight interval excludes the independent source residue')
        comparisons.append({'line': observed['line'], 'frequency_interval': source['omega_interval'],
                            'inferred_weight_interval': observed['weight_interval'],
                            'independent_source_weight_interval': source['positive_normalized_probe_weight_interval'],
                            'independent_exact_weight_contained': True})
    return {'settings': settings, 'retained_dyadic_power_matrices': powers,
            'weight_reconstruction': reconstructed, 'independent_residue_comparison': comparisons,
            'total_simulated_reset_trials': 16*TRIALS,
            'sampling_seed': SEED, 'sampling_generator': 'NumPy PCG64; binary zero-versus-all-other complete outcomes',
            'simultaneous_sampling_failure_bound': str(32*F(3, 8)**9),
            'sampling_error_per_setting': str(SAMPLING),
            'systematic_error_per_setting': str(SYSTEMATIC),
            'sampled_probability_convention': 'binary64 proxy to the certified powered ideal source word; its distance from the complete native instrument is explicitly bounded',
            'actual_raw_native_long_word_executed': False,
            'complete_native_instrument_probability_error_certified': True,
            'all_eight_source_support_weights_resolved': True,
            'unknown_frequency_support_inferred': False,
            'physical_particle_species_assigned': 0}
