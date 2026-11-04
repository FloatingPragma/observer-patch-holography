"""Complete Gaussian core readout in the original native pulse language."""
from collections import Counter
from fractions import Fraction as F
from itertools import product

import numpy as np
from scipy.linalg import expm

from ..certified_reference import certificate as cert
from ..centered_native_preparation import tree
from ..native_spin_observer import circuits, execute, precision


def rotations(dimensions=3):
    if type(dimensions) is not int or dimensions not in (1, 3):
        raise ValueError('the exact decoder uses the three-site chain or cube')
    sites = tuple(product(range(3), repeat=dimensions))
    result = []
    for axis in range(dimensions):
        for transverse in product(range(3), repeat=dimensions-1):
            line = []
            for coordinate in range(3):
                x = list(transverse)
                x.insert(axis, coordinate)
                line.append(sites.index(tuple(x)))
            for spin in range(2):
                result.append((2*line[0]+spin, 2*line[1]+spin, F(-1, 2)))
                result.append((2*line[1]+spin, 2*line[2]+spin, F(1, 4)))
    return tuple(result)


def givens(modes, left, right, pi_angle):
    if (type(modes) is not int or type(left) is not int or type(right) is not int
            or not 0 <= left < right < modes or type(pi_angle) is not F):
        raise ValueError('ordered CAR modes and an exact multiple of pi required')
    # a_p^* a_q-a_q^* a_p = i(X_p Z...Y_q-Y_p Z...X_q)/2.
    word = []
    for a, b, sign in [('X', 'Y', 1), ('Y', 'X', -1)]:
        label = ['I']*modes
        label[left], label[right] = a, b
        label[left+1:right] = ['Z']*(right-left-1)
        word += cert.pauli(label, cert.PiMultiple(sign*pi_angle/2))
    return word


def decode(dimensions=3):
    return [event for left, right, angle in rotations(dimensions)
            for event in givens(2*3**dimensions, left, right, angle)]


def readout(dimensions=3):
    modes = 2*3**dimensions
    center = 3**dimensions//2
    forward = decode(dimensions)
    return (forward+cert.cnot(2*center, modes)
            +cert.cnot(2*center+1, modes+1)+cert.inverse(forward))


def one_body_decoder(dimensions=3):
    size = 2*3**dimensions
    result = np.eye(size)
    for left, right, angle in rotations(dimensions):
        theta = float(angle)*np.pi
        rotation = np.array([[np.cos(theta), np.sin(theta)],
                             [-np.sin(theta), np.cos(theta)]])
        result[[left, right]] = rotation@result[[left, right]]
    return result


def analytic_error(template):
    word = list(circuits.emit(template))
    _, angle, variables = cert.consume(word, 0, template)
    if variables:
        raise ArithmeticError('decoder contains an uncertified variable angle')
    native = tree.port_bound(Counter((event['kind'], event['angle']) for event in word))
    return {'pulses': len(word), 'word_sha256': cert.word_digest(word),
            'fixed_binary64_angle_error': str(angle),
            'fixed_binary64_port_error': str(native),
            'complete_isometry_error': str(angle+native),
            'complete_unhalved_trace_error': str(2*(angle+native))}


def _annihilator(modes, mode):
    value = np.zeros((2**modes, 2**modes))
    for mask in range(2**modes):
        if mask & (1 << mode):
            value[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return value


def chain_target():
    modes = [_annihilator(6, mode) for mode in range(6)]
    gaussian = np.eye(64, dtype=complex)
    for left, right, angle in rotations(1):
        generator = modes[left].T@modes[right]-modes[right].T@modes[left]
        gaussian = expm(float(angle)*np.pi*generator)@gaussian
    initial = np.zeros((256, 64), complex)
    initial[:64] = np.eye(64)/8
    target = np.zeros_like(initial)
    for record in range(4):
        projector = np.diag([int(((mask >> 2) & 3) == record) for mask in range(64)])
        target[record*64:(record+1)*64] = gaussian.conj().T@projector@gaussian/8
    return initial, target


def native_execution():
    template = readout(1)
    initial, target = chain_target()
    # The inner fused contraction accepts 32 columns. Reassemble two pieces
    # of the same linear isometry; no density matrix or reference dephasing
    # intervenes between these amplitude columns.
    runs = [execute.complete(circuits.emit(template), 8, initial[:, start:start+32])
            for start in (0, 32)]
    success = np.column_stack([run.success for run in runs])
    failures = sum(run.failure_weights for run in runs)
    failure_probability = float(failures.sum())
    residual = sum(run.normalization_residual for run in runs)
    distance = float(np.sqrt(np.linalg.norm(success-target)**2+failure_probability))
    bound = F(analytic_error(template)['complete_isometry_error'])
    if distance > float(bound)+2e-11 or residual > 2e-11:
        raise ArithmeticError('complete native decoder misses the full-CAR target')
    if success.shape != (256, 64) or failures.shape != (len(template), 4):
        raise ArithmeticError('native decoder lost carrier or failure outputs')
    return {'field_modes': 6, 'full_field_rows': 64, 'quantum_reference_columns': 64,
            'work_rows_including_two_pointer_bits': 256, 'linear_reference_column_batches': 2,
            'native_pulses': len(template), 'complete_logical_pointer_outcomes': 4,
            'physical_two_pointer_outcomes': 36,
            'ideal_zero_reserved_pointer_outcomes': 32,
            'native_first_failure_slots': runs[0].native_failure_slots,
            'complete_output_error_diagnostic': distance,
            'complete_first_failure_probability': failure_probability,
            'normalization_residual': residual,
            'reference_columns_dephased': False,
            'success_postselected_or_renormalized': False,
            'all_first_failure_weights_retained_during_execution': True}


def certificate():
    cube = readout()
    chain = readout(1)
    for dimensions in (1, 3):
        orbital = np.array([1., 0., 1.])/np.sqrt(2)
        phi = orbital
        for _ in range(dimensions-1):
            phi = np.kron(phi, orbital)
        zero = np.kron(phi[:, None], np.eye(2))
        target = np.zeros_like(zero)
        center = 3**dimensions//2
        target[2*center:2*center+2] = np.eye(2)
        if np.linalg.norm(one_body_decoder(dimensions)@zero-target) > 1e-12:
            raise ArithmeticError('the Gaussian decoder fails its source zero orbital')
    if len(rotations()) != 108 or len(cube) != 24630 or len(chain) != 278:
        raise ArithmeticError('complete native decoder resource count changed')
    bits = precision.bits_for(len(cube), F(1, 10**9))
    commands = {}
    for kind, _, angle in cube:
        key = kind+':'+str(angle.value)
        if key not in commands:
            commands[key] = precision.dyadic_command(kind, angle, bits)
    return {'cube_field_modes': 54, 'cube_full_field_dimension': 2**54,
            'cube_Givens_rotations_per_decode': 108,
            'cube_decode_native_pulses': len(decode()),
            'cube_complete_four_outcome_readout_native_pulses': len(cube),
            'complete_Jordan_Wigner_parity_strings_retained': True,
            'cube_readout_compiled_and_local_identities_verified': True,
            'cube_full_Fock_native_readout_executed': False,
            'chain_complete_native_error_certificate': analytic_error(chain),
            'cube_complete_native_error_certificate': analytic_error(cube),
            'refined_dyadic_fraction_bits': bits,
            'refined_complete_unhalved_trace_error_upper': str(F(58*len(cube), 2**bits)),
            'distinct_refined_commands': commands,
            'refined_rational_commands_run_through_binary64_API': False,
            'two_blank_memory_bits_and_addressed_native_pulses_required': True,
            'interleaved_native_repair_protection_proved': False,
            'simultaneous_free_Hamiltonian_drift_compensated': False,
            'interacting_spectral_projectors_measured_exactly': False}
