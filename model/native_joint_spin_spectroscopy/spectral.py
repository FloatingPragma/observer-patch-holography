"""Certified positive joint bins from the retained full finite eigensystem.

All tabulated intervals use exact rational arithmetic. The eigenvectors are
untrusted dyadic proposals admitted by the original complete certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
import copy
from math import factorial
import numpy as np
from model.native_ordered_ground import source, certificate, experiment
from model.native_spin_observer.precision import pi_interval
from model.centered_particle_spectrum.certificate import decimal_bounds
from . import algebra


def cosine_interval(midpoint, pi_coefficient=F()):
    if type(midpoint) is not F or type(pi_coefficient) is not F:
        raise ValueError('exact rational angle coefficients required')
    return _cosine_interval(midpoint, pi_coefficient)


@lru_cache(maxsize=100000)
def _cosine_interval(midpoint, pi_coefficient):
    plo, phi = pi_interval(110)
    scale = 1 << 110
    plo = F((plo * scale).numerator // (plo * scale).denominator, scale)
    upper = phi * scale
    phi = F(-((-upper.numerator) // upper.denominator), scale)
    pm = (plo + phi) / 2
    center = midpoint + pi_coefficient * pm
    turns = (center / (2 * pm) + F(1, 2)).numerator // (center / (2 * pm) + F(1, 2)).denominator
    coefficient = pi_coefficient - 2 * turns
    x = midpoint + coefficient * pm
    uncertainty = abs(coefficient) * (phi - plo) / 2
    scaled = x * (1 << 64) + F(1, 2)
    rounded = F(scaled.numerator // scaled.denominator, 1 << 64)
    uncertainty += abs(x - rounded)
    x = rounded
    if abs(x) > 4:
        raise ArithmeticError('certified range reduction failed')
    value = sum(((-1) ** k * x ** (2 * k) / factorial(2 * k)
                 for k in range(17)), F())
    remainder = abs(x) ** 34 / factorial(34) + uncertainty
    return max(F(-1), value - remainder), min(F(1), value + remainder)


def fejer_interval(energy, bin_index, size=4, tau=F(1)):
    if type(size) is not int or size < 2 or type(bin_index) is not int or not 0 <= bin_index < size:
        raise ValueError('finite complete clock and actual bin required')
    if type(energy) is not F or type(tau) is not F or tau <= 0:
        raise ValueError('exact energy and positive time required')
    lo = hi = F(1, size)
    for n in range(1, size):
        a, b = cosine_interval(n * tau * energy, F(-2 * n * bin_index, size))
        c = F(2 * (size - n), size * size)
        lo += c * a
        hi += c * b
    return max(F(), lo), min(F(1), hi)


def operator_error(size, tau, radius, origin_error, gram_defect, state_error):
    if type(size) is not int or size < 2:
        raise ValueError('complete clock size required')
    if any(type(x) is not F or x < 0 for x in (tau, radius, origin_error, gram_defect, state_error)) or tau == 0:
        raise ValueError('nonnegative exact error inputs and positive tau required')
    if gram_defect >= 1:
        raise ValueError('polar frame must be invertible')
    return tau * F(size * size - 1, 3 * size) * (radius + origin_error) + gram_defect * (2 + gram_defect) + 2 * state_error


def cube_readout(size=4, tau=F(1, 2)):
    if type(size) is not int or size < 2 or type(tau) is not F or tau <= 0:
        raise ValueError('exact clock size and positive rational time required')
    return copy.deepcopy(_cube_readout(size, tau))


@lru_cache(maxsize=8)
def _cube_readout(size, tau):
    w = experiment.witness()
    checked = certificate.verify(w)
    h, scale = source.matrix()
    group = algebra.audit(h, int(scale))
    integers = w['vectors'].astype(object)
    ground = integers[:, 0]
    norm = sum(int(x) ** 2 for x in ground)
    d = 1 << certificate.BITS
    energies = [F(int(x), d) for x in w['energies']]
    origin = energies[0]
    e0lo, e0hi = map(F, checked['ground_energy_interval'])
    origin_error = max(abs(origin - e0lo), abs(origin - e0hi))
    error = operator_error(size, tau, F(checked['complete_polar_Weyl_radius']),
                           origin_error, F(checked['Gram_defect']),
                           F(checked['normalized_ground_vector_error']))
    probe = ground * np.array([1 - 2 * (m & 1) for m in range(256)], dtype=object)
    kernels = [[fejer_interval(e - origin, j, size, tau) for e in energies]
               for j in range(size)]
    rows = []
    for k in range(8):
        q = algebra.character_projection(probe, k)
        coefficients = integers.T @ q
        weights = [F(int(x) ** 2, d * d * 64 * norm) for x in coefficients]
        momentum = F(sum(int(x) ** 2 for x in q), 64 * norm)
        output = []
        for j in range(size):
            lo = sum((weight * interval[0] for weight, interval in zip(weights, kernels[j], strict=True)), F())
            hi = sum((weight * interval[1] for weight, interval in zip(weights, kernels[j], strict=True)), F())
            output.append(decimal_bounds(max(F(), lo - error), min(F(1), hi + error), 11))
        state_error = F(checked['normalized_ground_vector_error'])
        rows.append({'character_bits': format(k, '03b'),
                     'momentum_probability_interval': decimal_bounds(max(F(), momentum - 2 * state_error), min(F(1), momentum + 2 * state_error), 11),
                     'energy_bin_probability_intervals': output})
    excitation_upper = energies[-1] + F(checked['complete_polar_Weyl_radius']) - e0lo
    alias_free = tau * excitation_upper < 2 * pi_interval(110)[0]
    return {'source': 'FA98 selected open two-cube: J=1, B=1/4; local encoded Z probe',
            'group_audit': group, 'clock_size': size, 'tau': str(tau),
            'complete_excitation_energy_upper': str(excitation_upper),
            'all_source_excitation_energies_in_one_phase_period': alias_free,
            'per_bin_certified_operator_and_state_error': str(error),
            'positive_complete_joint_instrument': True,
            'open_cube_is_periodic_dispersion_sample': False,
            'relativistic_mass_identified': False,
            'rows': rows}
