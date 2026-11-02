"""Literal source words for joint clocks on the protected spin encoding."""
from fractions import Fraction as F
from model.certified_reference import certificate as cert
from model.certified_reference.continuation import divide
from model.native_spin_observer import circuits as spin
from model.native_spectral_detector.circuits import hadamard_template, controlled_phase_template


def scale(angle, coefficient):
    return cert.PiMultiple(angle.value * coefficient) if type(angle) is cert.PiMultiple else angle * coefficient


def controlled_pauli(label, angle, control):
    if type(control) is not int or not 0 <= control < len(label) or label[control] != 'I':
        raise ValueError('a separate blank control factor is required')
    with_control = label[:control] + 'Z' + label[control + 1:]
    return cert.pauli(label, divide(angle, 2)) + cert.pauli(with_control, -divide(angle, 2))


def label(axis, selected, sites, bits):
    result = ['I'] * bits
    for site in selected:
        if not 0 <= site < sites:
            raise ValueError('actual spin site required')
        result[2 * site] = axis
        if axis != 'Z':
            result[2 * site + 1] = 'X'
    return ''.join(result)


def controlled_exchange(time, control, sites=2, seam=(0, 1), bits=7):
    word = []
    for axis in 'XYZ':
        word += controlled_pauli(label(axis, seam, sites, bits), scale(time, F(-1, 4)), control)
    return word + controlled_pauli('I' * bits, scale(time, F(1, 4)), control)


def combined_dimer_reflection(control=4, bits=7):
    # exp(-i*pi*h/J)=SWAP for h=J(S.S-1/4), including the scalar phase.
    word = controlled_exchange(cert._pi(1), control, bits=bits)
    word += controlled_pauli(label('Z', (0, 1), 2, bits), cert._pi(1, 2), control)
    return word + controlled_pauli('I' * bits, cert._pi(-1, 2), control)


def dimer_evolution(duration, steps, control, field=F(1, 4), bits=7):
    if type(duration) is not F or duration < 0 or type(steps) is not int or steps < 1:
        raise ValueError('exact nonnegative duration and positive step count required')
    dt = duration / steps
    one = controlled_exchange(dt, control, bits=bits)
    for site, sign in ((0, 1), (1, -1)):
        one += controlled_pauli(label('X', (site,), 2, bits), sign * field * dt / 2, control)
    return one * steps


def phase_origin(duration, energy, control, bits=7):
    return controlled_pauli('I' * bits, duration * energy, control)


def positive_fourier(clock_bits):
    if not clock_bits or len(set(clock_bits)) != len(clock_bits):
        raise ValueError('distinct clock bits required')
    word = []
    for high in reversed(range(len(clock_bits))):
        word += hadamard_template(clock_bits[high])
        for low in reversed(range(high)):
            word += spin.remap(controlled_phase_template(F(1, 1 << (high - low))),
                               (clock_bits[low], clock_bits[high]))
    for low in range(len(clock_bits) // 2):
        a, b = clock_bits[low], clock_bits[-1 - low]
        word += cert.cnot(a, b) + cert.cnot(b, a) + cert.cnot(a, b)
    return word


def dimer_clock(origin, steps=(32, 64)):
    if type(origin) is not F or len(steps) != 2:
        raise ValueError('certified rational origin and two step counts required')
    word = hadamard_template(4) + combined_dimer_reflection() + hadamard_template(4)
    for bit in (5, 6):
        word += hadamard_template(bit)
    for j, count in enumerate(steps):
        duration = F(1 << j)
        word += dimer_evolution(duration, count, 5 + j)
        word += phase_origin(duration, origin, 5 + j)
    return word + positive_fourier((5, 6))
