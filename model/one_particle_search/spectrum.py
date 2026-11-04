"""Exact complete search on the unchanged side-eight one-particle source.

All shell decisions use Q(sqrt(2)) and integer arithmetic. Floating arrays
implement the independently checked ideal detector, never the certificate.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb

import numpy as np
import sympy as s

from ..particle_species_detector import detector, native_fft

SIDE, VOLUME, MODES = 8, 512, 1024
# Twice the sine at pi*k/4, as a+b*sqrt(2).
TWICE_SINE = ((0, 0), (0, 1), (2, 0), (0, 1),
              (0, 0), (0, -1), (-2, 0), (0, -1))


def coordinates(momentum):
    if (type(momentum) not in (tuple, list) or len(momentum) != 3 or
            any(type(k) is not int or not 0 <= k < SIDE for k in momentum)):
        raise ValueError('three integer momentum labels in 0..7 required')
    return tuple(momentum)


def direction(momentum):
    return tuple((-a, -b) for a, b in (TWICE_SINE[k] for k in coordinates(momentum)))


def shell(momentum):
    return sum(a*a+2*b*b for a, b in direction(momentum))//2


def symbolic_block(momentum):
    x, y, z = (a+b*s.sqrt(2) for a, b in direction(momentum))
    return s.Matrix([[z, x-s.I*y], [x+s.I*y, -z]])


def exact_certificate():
    counts = Counter()
    sums = {n: [[0, 0] for _ in range(3)] for n in range(7)}
    checked = set()
    for k in product(range(SIDE), repeat=3):
        n = shell(k)
        counts[n] += 1
        d = direction(k)
        for j, (a, b) in enumerate(d):
            sums[n][j][0] += a
            sums[n][j][1] += b
        if d not in checked:
            h = symbolic_block(k)
            if h != h.conjugate().T or h.trace() != 0 or s.expand(h*h) != 2*n*s.eye(2):
                raise ArithmeticError('exact source block identity failed')
            checked.add(d)
    if counts != {n: 8*comb(6, n) for n in range(7)} or any(
            a or b for rows in sums.values() for a, b in rows):
        raise ArithmeticError('full spectrum multiplicity or local spin sum failed')
    rows = []
    for n in range(7):
        multiplicity = 2*counts[n] if n == 0 else counts[n]
        diagonal = F(multiplicity, MODES)
        rows.append({'shell_n': n, 'energy': '0' if n == 0 else '+/-sqrt(2*n)/16',
                     'energy_squared': str(F(n, 128)),
                     'signed_branches': 1 if n == 0 else 2,
                     'multiplicity_per_energy': multiplicity,
                     'one_site_occupancy_upper': str(diagonal),
                     'necessary_sites_for_exact_support': (MODES+multiplicity-1)//multiplicity})
    return {'all_momenta': VOLUME, 'all_single_particle_outcomes': MODES,
            'distinct_exact_blocks_checked': len(checked),
            'shells': rows, 'zero_energy_modes': 16,
            'complete_count': sum(r['signed_branches']*r['multiplicity_per_energy'] for r in rows),
            'spectrum_polynomial': 'd_n = [z^n](2+4z+2z^2)^3 = 8*binomial(6,n)',
            'scaled_block_identity': '(16*h(k))^2=2*n*I2; trace(h(k))=0',
            'projectors_for_n_positive': '(I2 +/- 16*h(k)/sqrt(2*n))/2',
            'zero_block_effects': 'two CAR component projectors, both assigned energy zero',
            'contact_on_complete_N1_sector': 'zero for every U; original H is unchanged',
            'single_particle_ground_energy': '-sqrt(12)/16',
            'original_vacuum_is_positive_energy_reference': False}


def rotation(momentum):
    """SU(2) rows are the +energy and -energy eigenbras, with polar charts."""
    n = shell(momentum)
    if not n:
        return np.eye(2, dtype=complex)
    x, y, z = (a+b*np.sqrt(2) for a, b in direction(momentum))
    if x == y == 0:
        return np.eye(2, dtype=complex) if z > 0 else np.array([[0, 1], [-1, 0]], complex)
    e = np.sqrt(2*n)
    norm = np.sqrt(2*e*(e+z))
    return np.array([[e+z, x-1j*y], [-x-1j*y, e+z]], complex)/norm


def energy_labels():
    values = []
    for k in product(range(SIDE), repeat=3):
        e = np.sqrt(2*shell(k))/16
        values.extend((e, -e))
    return np.asarray(values)


def apply(value, *, inverse=False):
    """All 1,024 ideal N=1 rows, including both components of every zero mode."""
    if type(inverse) is not bool or np.ma.isMaskedArray(value):
        raise ValueError('unmasked input and Boolean inverse required')
    raw = np.asarray(value)
    if (raw.ndim not in (1, 2) or raw.shape[0] != MODES or
            (raw.ndim == 2 and not 1 <= raw.shape[1] <= MODES) or
            raw.dtype.kind not in 'ficu' or not np.isfinite(raw).all()):
        raise ValueError('all 1024 finite one-particle rows required')
    out = np.array(raw, dtype=complex, copy=True)
    if not inverse:
        out = native_fft.apply(out)
    for i, k in enumerate(product(range(SIDE), repeat=3)):
        r = rotation(k)
        out[2*i:2*i+2] = (r.conjugate().T if inverse else r)@out[2*i:2*i+2]
    return native_fft.apply(out, inverse=True) if inverse else out


def source_replay():
    """Independent sparse real-space H on every column of the complete detector."""
    vectors = apply(np.eye(MODES), inverse=True)
    energies = energy_labels()
    residual = np.max(abs(detector.original_hamiltonian()@vectors-vectors*energies))
    inverse_error = np.max(abs(apply(vectors)-np.eye(MODES)))
    if residual > 2e-13 or inverse_error > 2e-13:
        raise ArithmeticError('complete source spectral detector failed')
    # These ceilings are diagnostics, not the proof or outward error bounds.
    return {'columns_executed': MODES,
            'original_H_eigenvector_residual_below': '2e-13',
            'complete_detector_inverse_residual_below': '2e-13',
            'probability_outputs_retained': MODES,
            'probability_renormalization': False,
            'floating_values_used_for_exact_multiplicities': False}


def localization_certificate():
    # A stationary density commutes with H and is a convex mixture of density
    # matrices on its distinct energy eigenspaces. Translation and exact spin
    # cancellation make P_E(x,x)=d_E/1024 I2 on each of those eigenspaces.
    peak = max(F(row['one_site_occupancy_upper']) for row in exact_certificate()['shells'])
    if peak != F(5, 32):
        raise ArithmeticError('stationary one-site bound changed')
    return {'stationary_one_site_occupancy_upper': str(peak),
            'any_m_site_region_occupancy_upper': 'min(1,5*m/32)',
            'six_site_exterior_occupancy_lower': '1/16',
            'one_site_bound_attained': 'normalized P_E |x,spin>, n=3, either sign',
            'finite_reference': 'any normalized density commuting with the unchanged N1 Hamiltonian',
            'finite_torus_point_spectrum_present': True,
            'infinite_c_positive_compact_stationary_eigenvectors': 0,
            'infinite_c_positive_square_summable_stationary_eigenvectors': 0,
            'infinite_c_positive_nonzero_positive_trace_class_stationary_densities': 0,
            'c_zero_countercontrol': 'H=0 on N1; every one-site spin state is stationary',
            'infinite_proof': 'nonzero Laurent determinant and Fourier/Fubini null level sets',
            'moving_packets_or_multibody_localization_excluded': False}


def laurent_certificate():
    """All energies at once: a nonzero z_x^2 coefficient survives at c != 0."""
    z = s.symbols('z0:3', nonzero=True)
    c, energy = s.symbols('c energy', real=True)
    sine = [(q-1/q)/(2*s.I) for q in z]
    x, y, zz = [-c*q for q in sine]
    h = s.Matrix([[zz, x-s.I*y], [x+s.I*y, -zz]])
    matrix = h-energy*s.eye(2)
    determinant = s.expand(matrix.det())
    expected = s.expand(energy**2-c**2*sum(q*q for q in sine))
    if determinant != expected or s.expand(matrix.adjugate()*matrix) != determinant*s.eye(2):
        raise ArithmeticError('Laurent adjugate identity failed')
    if s.expand(determinant).coeff(z[0], 2) != c**2/4:
        raise ArithmeticError('nonzero determinant coefficient missing')
    return {'coefficient_of_z0_squared': 'c^2/4',
            'domain': 'infinite Z^3, c>0, any real E, original N1 sector',
            'determinant': 'E^2-c^2*sum_j sin(k_j)^2',
            'compact_support_argument': 'Laurent polynomial ring is an integral domain',
            'L2_argument': 'for fixed ky,kz the determinant has at most four kx zeros; Fubini and Plancherel',
            'mixed_stationary_argument': 'a positive trace-class density has a nonzero finite-dimensional eigenspace preserved by bounded H, hence an H eigenvector',
            'zero_energy_included': True,
            'c_zero_excluded_from_no_go': True}
