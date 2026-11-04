"""Finite native source words for every momentum spin rotation.

The independent dyadic certificate checks literal angles against radical
matrix entries. High-precision acos/atan2 propose words; they prove no bound.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import mpmath as mp
import numpy as np

from ..certified_reference import certificate as cert, continuation as templates, interval as iv
from ..observer.source_scalar import compiler
from ..particle_species_detector import native_fft
from . import spectrum


def reciprocal(value):
    if value.lo <= 0:
        raise ValueError('strictly positive interval required')
    return iv.Interval(
        iv.Interval.point(F(iv.SCALE, value.hi)).lo,
        iv.Interval.point(F(iv.SCALE, value.lo)).hi)


def square_root(value):
    if value.lo < 0:
        raise ValueError('nonnegative interval required')
    return iv.Interval(iv.sqrt_fraction(F(value.lo, iv.SCALE)).lo,
                       iv.sqrt_fraction(F(value.hi, iv.SCALE)).hi)


def radical(value):
    a, b = value
    return iv.Interval.point(a)+b*iv.sqrt_fraction(F(2))


def target_entries(momentum):
    """Enclosures of cosine and real/imaginary upper off-diagonal of R(k)."""
    n = spectrum.shell(momentum)
    dx, dy, dz = spectrum.direction(momentum)
    if not n or (dx == dy == (0, 0) and (dz[0] > 0 or dz[1] > 0)):
        return tuple(iv.Interval.point(q) for q in (1, 0, 0))
    if dx == dy == (0, 0):
        return tuple(iv.Interval.point(q) for q in (0, 1, 0))
    x, y, z = map(radical, (dx, dy, dz))
    e = iv.sqrt_fraction(F(2*n))
    inverse_norm = reciprocal(square_root(2*e*(e+z)))
    return (e+z)*inverse_norm, x*inverse_norm, -y*inverse_norm


def rotation_word(momentum):
    """Chronological D* exp(-i t sigma_x) D in the complete local Fock block."""
    n = spectrum.shell(momentum)
    if not n:
        return (), F(), F()
    # Fixed high precision avoids platform-libm differences in proposed
    # literal words. The interval check below establishes their accuracy.
    with mp.workdps(80):
        x, y, z = (mp.mpf(a)+b*mp.sqrt(2) for a, b in spectrum.direction(momentum))
        if x == y == 0 and z > 0:
            return (), F(), F()
        if x == y == 0:
            alpha, time = float(mp.pi/2), float(mp.pi/2)
        else:
            alpha = float(mp.atan2(x, y))
            time = float(mp.acos(z/mp.sqrt(2*n))/2)
    def phase(angle):
        return [compiler.pulse('I', angle/2, [0]), compiler.pulse('Z', -angle/2, [0])]
    word = (phase(-alpha)+compiler.pauli_rotation('XX', -time/2)+
            compiler.pauli_rotation('YY', -time/2)+phase(alpha))
    return tuple(word), F(alpha), F(time)


@lru_cache(maxsize=512)
def rotation_certificate(momentum):
    word, alpha, time = rotation_word(momentum)
    if not word:
        return {'pulses': 0, 'abort_slots': 0, 'error': F(),
                'target_rotation_error': F(), 'fixed_angle_error': F(), 'source_control_error': F()}
    target_c, target_r, target_i = target_entries(momentum)
    st, ct = iv.sine_cosine(time)
    sa, ca = iv.sine_cosine(alpha)
    # -i exp(i alpha) sin(t) = sin(alpha)sin(t)-i cos(alpha)sin(t).
    entry_error = ((ct-target_c).absolute_upper()+
                   (sa*st-target_r).absolute_upper()+
                   (-ca*st-target_i).absolute_upper())
    # The full Fock extension is diag(1,R,1). Its Frobenius norm is at
    # most twice the sum of the displayed entry errors.
    rotation_error = 2*entry_error
    template = templates.phase(0, -alpha)+templates.real_hopping(time)+templates.phase(0, alpha)
    end, fixed, variables = cert.consume(word, 0, template)
    if end != len(word) or variables:
        raise ArithmeticError('unconsumed source spin-rotation event')
    port = cert.controls_bound(word, {})
    return {'pulses': len(word),
            'abort_slots': sum(6-(1 << len(e['targets'])) for e in word),
            'error': rotation_error+fixed+port, 'target_rotation_error': rotation_error,
            'fixed_angle_error': fixed, 'source_control_error': port}


def certificate():
    rows = [rotation_certificate(k) for k in product(range(8), repeat=3)]
    error = sum((row['error'] for row in rows), F())
    fft = native_fft.certificate()
    total = (F(fft['forward_complete_Fock_source_error'])+error+
             F(fft['all_mode_occupation_pointer_source_error']))
    if total >= F(1, 10**8):
        raise ArithmeticError('complete detector source error exceeds declared ceiling')
    return {'spin_rotations': sum(row['pulses'] > 0 for row in rows),
            'spin_rotation_native_pulses': sum(row['pulses'] for row in rows),
            'spin_rotation_first_export_abort_slots': sum(row['abort_slots'] for row in rows),
            'spin_rotation_complete_Fock_operator_error': str(error),
            'detector_native_pulses': fft['native_pulses_per_Fourier_transform']+
                sum(row['pulses'] for row in rows)+fft['occupation_pointer_native_pulses'],
            'detector_first_export_abort_slots': fft['first_export_abort_slots_per_Fourier_transform']+
                sum(row['abort_slots'] for row in rows)+fft['occupation_pointer_first_export_abort_slots'],
            'complete_detector_operator_error': str(total),
            'complete_output_unhalved_L1_error_upper': str(2*total),
            'error_scope': 'all Fock inputs and arbitrary untouched references; full source exports retained',
            'mode_order': '2*((kx*8+ky)*8+kz)+component; rotation pairs are adjacent',
            'pointer_assumption': 'fresh blank pointer for every occupation read',
            'full_dictionary': 'all occupation strings and all first source/export abort records',
            'N1_energy_labels': 1024,
            'higher_N_outcomes_are_interacting_energy_labels': False,
            'full_2_to_1024_native_state_array_executed': False,
            'complete_source_error_projected_to_N1': False}


def local_execution():
    """Four-column source executions, including vacuum and double occupancy."""
    unique = {}
    for k in product(range(8), repeat=3):
        unique.setdefault(spectrum.direction(k), k)
    count = 0
    for k in unique.values():
        word, _, _ = rotation_word(k)
        actual = compiler.execute(list(word), 2) if word else np.eye(4)
        expected = np.eye(4, dtype=complex)
        expected[1:3, 1:3] = spectrum.rotation(k)
        if np.max(abs(actual-expected)) > 2e-13:
            raise ArithmeticError('complete four-column spin rotation execution failed')
        count += 1
    return {'distinct_complete_local_blocks_executed': count,
            'columns_per_local_block': 4,
            'vacuum_and_double_occupation_retained': True,
            'execution_residual_ceiling': '2e-13',
            'execution_residual_is_certificate': False}
