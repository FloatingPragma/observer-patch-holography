"""Source-executed static quasiparticles on the colored native two-particle lift.

The internal Hamiltonian and charge chart are supplied by the separate
Bogoliubov certificate. Native record j selects a declared duration multiplier
1+j/4 and an electromagnetic phase. This is a designed joint instrument,
not an exterior power of a nonunitary Kraus operator.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np

from ..common_source import record_code
from ..defect_protection import service as cover
from ..exchange_sector import statistics as native, statistics_service as source
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler


BITS = 7
SIZE = 1 << BITS
SOURCE_PAIR = (0, 1)
SEAM = (0, 5)
ADMISSION = 2e-7


def _parameters(energies, charges, duration, electromagnetic):
    e, q = np.asarray(energies), np.asarray(charges)
    if (e.ndim != 1 or q.shape != e.shape or not 2 <= len(e) <= SIZE
            or e.dtype.kind not in 'if' or q.dtype.kind not in 'if'
            or not np.isfinite(e).all() or not np.isfinite(q).all()
            or np.min(e) < 0):
        raise ValueError('two through128 finite nonnegative energies and matching charges required')
    if any(isinstance(v, (bool, str)) or not math.isfinite(float(v))
           for v in (duration, electromagnetic)) or float(duration) <= 0:
        raise ValueError('positive finite duration and finite electromagnetic phase required')
    return e.astype(float), q.astype(float), float(duration), float(electromagnetic)


def phase_word(energies, charges, duration, electromagnetic, current, *, vacuum_energy=0.):
    e, q, t, angle = _parameters(energies, charges, duration, electromagnetic)
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('an exact signed native current required')
    phases = np.zeros(SIZE)
    if isinstance(vacuum_energy, (bool, str)) or not math.isfinite(float(vacuum_energy)):
        raise ValueError('a finite specified vacuum energy required')
    # On the two-quasiparticle sector the vacuum scalar occurs once. Assign
    # half to each one-particle factor, retaining its full complex phase.
    phases[:len(e)] = -t*(1+current/4)*(e+float(vacuum_energy)/2)+angle*current*q
    return compiler.compile_diagonal(phases), phases


def phase_gate(energies, charges, duration, electromagnetic, current, *, vacuum_energy=0.):
    word, phases = phase_word(energies, charges, duration, electromagnetic, current,
                             vacuum_energy=vacuum_energy)
    actual = compiler.execute(word, BITS)
    target = np.diag(np.exp(1j*phases))
    if np.linalg.norm(actual-target, 2) > ADMISSION:
        raise ArithmeticError('complete128-state quasiparticle phase source execution failed')
    return actual, word


def charge_probe(charges, angle=1., quadrature=math.pi/2):
    """Execute a neutral-ancilla QND charge Ramsey gate on all256 columns."""
    q = np.asarray(charges, dtype=float)
    if (q.ndim != 1 or not 2 <= len(q) <= SIZE or not np.isfinite(q).all()
            or not math.isfinite(angle) or not math.isfinite(quadrature)):
        raise ValueError('finite charge chart and charge probe phases required')
    padded = np.zeros(SIZE)
    padded[:len(q)] = q
    phases = np.zeros(2*SIZE)
    phases[SIZE:] = angle*padded+quadrature
    word = ([compiler.pulse('Y', -math.pi/4, [BITS])]
            + compiler.compile_diagonal(phases)
            + [compiler.pulse('Y', math.pi/4, [BITS])])
    actual = compiler.execute(word, BITS+1)
    z = np.exp(1j*phases[SIZE:])
    a, b = np.diag((1+z)/2), np.diag((z-1)/2)
    target = np.block([[a, b], [b, a]])
    if np.linalg.norm(actual-target, 2) > ADMISSION:
        raise ArithmeticError('complete neutral-ancilla charge probe changed')
    return actual, word, {
        'source_full_input_columns': 2*SIZE,
        'operator_error_below': str(ADMISSION),
        'angle': angle, 'quadrature': quadrature,
        'charge_values': sorted(set(map(float, q))),
        'pointer_zero_probability_by_charge': [
            {'charge': float(v), 'probability': (1+math.cos(angle*v+quadrature))/2}
            for v in sorted(set(q))],
        'neutral_pointer_preserves_charge': True,
        'measured_EM_field_response': False}


def fredkin_word():
    """Control on local bit0; swap bits1,2, with all eight inputs retained."""
    return circuit.basis_swap_word(3, 5, 3)


def encoder_word():
    """|0> psi -> (|0>psi-|1>SWAP psi)/sqrt2 on two seven-bit registers."""
    word = [compiler.pulse('Y', -math.pi/4, [14])]
    for bit in range(BITS):
        word += compiler.remap(fredkin_word(), [14, bit, BITS+bit])
    word += [compiler.pulse('Z', math.pi/2, [14]),
             compiler.pulse('I', -math.pi/2, [14])]
    return word


def encoder_certificate():
    actual = compiler.execute(fredkin_word(), 3)
    target = np.eye(8, dtype=complex)
    target[[3, 5]] = target[[5, 3]]
    error = float(np.linalg.norm(actual-target, 2))
    if error > ADMISSION:
        raise ArithmeticError('controlled internal swap lost its source truth table')
    rotation = compiler.execute([compiler.pulse('Y', -math.pi/4, [0])], 1)
    rotation_target = np.array([[1., -1.], [1., 1.]])/math.sqrt(2)
    corrected_z = compiler.execute([compiler.pulse('Z', math.pi/2, [0]),
                                    compiler.pulse('I', -math.pi/2, [0])], 1)
    if (np.linalg.norm(rotation-rotation_target, 2) > ADMISSION
            or np.linalg.norm(corrected_z-np.diag([1., -1.]), 2) > ADMISSION):
        raise ArithmeticError('encoder one-bit factors changed')
    return {'controlled_swap_full_input_columns': 8,
            'controlled_swap_source_error_below': str(ADMISSION),
            'one_bit_rotation_and_corrected_Z_full_columns_checked': 2,
            'local_factor_operator_error_below': str(ADMISSION),
            'controlled_swaps_per_encode_or_decode': BITS,
            'encoder_resources': compiler.resources(encoder_word()),
            'decoder_resources': compiler.resources(compiler.inverse(encoder_word()))}


def ramsey_word(first, second):
    if (type(first) is not int or type(second) is not int
            or not 0 <= first < second < SIZE):
        raise ValueError('two increasing distinct padded internal addresses required')
    return compiler._gray_rotation(first, second, math.pi/4, BITS)


def _x_word(address, offset):
    result = []
    for bit in range(BITS):
        if address >> bit & 1:
            result += [compiler.pulse('X', math.pi/2, [offset+bit]),
                       compiler.pulse('I', -math.pi/2, [offset+bit])]
    return result


def prepare_internal(first, second, partner):
    """Actually execute all encoder pulses on the full32768-address state."""
    if type(partner) is not int or not 0 <= partner < SIZE:
        raise ValueError('a padded partner address required')
    read = ramsey_word(first, second)
    word = (_x_word(first, 0)+_x_word(partner, BITS)
            + compiler.inverse(read)+encoder_word())
    blank = np.zeros(1 << 15, complex)
    blank[0] = 1
    actual = compiler.execute(word, 15, blank)
    # Amplitude matrices use (first particle, second particle); the physical
    # little-endian register has first particle in its seven low wires.
    matrices = actual.reshape(2, SIZE, SIZE).transpose(0, 2, 1)
    rotation = compiler.execute(read, BITS)
    a = rotation.conj().T[:, first]
    b = np.eye(SIZE, dtype=complex)[:, partner]
    target = np.outer(a, b)
    if (np.linalg.norm(matrices[0]-target/math.sqrt(2)) > ADMISSION
            or np.linalg.norm(matrices[1]+target.T/math.sqrt(2)) > ADMISSION):
        raise ArithmeticError('colored negative-exchange encoder failed complete state comparison')
    return matrices.copy(), word, rotation


def prepare_source():
    permutation = source._witness_permutation((SOURCE_PAIR,))
    blank = np.zeros(512, complex)
    blank[0] = 1
    value = source._permutation(permutation, blank)
    target = np.eye(512, dtype=complex)[:, source._address(SOURCE_PAIR)]
    if np.linalg.norm(value-target) > ADMISSION:
        raise ArithmeticError('native source preparation failed')
    columns = np.zeros((1024, 2), complex)
    columns[:512, 0] = value
    columns[512:, 1] = value
    return columns, record_code._permutation_resources(permutation)


def tensor_norm_squared(vectors, matrices):
    """Norm of sum_t vectors[:,t] tensor matrices[t], with no rank truncation."""
    flat = np.asarray(matrices).reshape(len(matrices), -1)
    # Evaluate in bounded blocks as sums of squared amplitudes. This avoids
    # negative cancellation artifacts of a Gram contraction, without clipping.
    total = 0.
    for start in range(0, flat.shape[1], 256):
        amplitudes = vectors @ flat[:, start:start+256]
        total += float(np.vdot(amplitudes, amplitudes).real)
    return total


def decoder_factors(native_columns, matrices):
    """Exact contraction of the ideal inverse encoder, including its failure bit.

    Both decoded orientation outcomes are retained. On the exact negative
    colored code the second is zero; floating source residual is not erased.
    The report bounds replacing this ideal map by its locally checked source
    word; this function does not execute that large-register word directly.
    """
    c0, c1 = matrices
    b0, b1 = native_columns[:512], native_columns[512:]
    internal = np.stack((c0, c1, c0.T, c1.T))
    first = np.column_stack((b0, -b1))/math.sqrt(2)
    second = np.column_stack((-b0, -b1))/math.sqrt(2)
    return (first, second), internal


def pointer_probabilities(vectors, internal):
    flat = internal.reshape(len(internal), -1)
    probabilities = np.zeros(SIZE*SIZE)
    for start in range(0, len(probabilities), 256):
        amplitudes = vectors @ flat[:, start:start+256]
        probabilities[start:start+256] = np.sum(abs(amplitudes)**2, axis=0)
    return probabilities.reshape(SIZE, SIZE)


def _sum_resources(*rows):
    keys = source.COUNT_KEYS
    return {key: sum(row[key] for row in rows) for key in keys}


def execute(energies, charges, *, ramsey_indices, partner, duration=F(1),
            electromagnetic=F(1, 8), vacuum_energy=0.):
    """All C/E raw outputs of a colored-particle native-feedback Ramsey test.

    The supplied chart is an invariant static positive-energy quasiparticle
    sector, not the bare particle-number sector of a pairing Hamiltonian.
    Returned factors determine every final source amplitude for every pointer
    outcome; complete computational probabilities are included separately.
    """
    e, q, t, angle = _parameters(energies, charges, duration, electromagnetic)
    if (type(ramsey_indices) is not tuple or len(ramsey_indices) != 2
            or not all(type(i) is int and 0 <= i < len(e) for i in (*ramsey_indices, partner))):
        raise ValueError('valid quasiparticle Ramsey and partner indices required')
    first, second = ramsey_indices
    if q[first] != q[second] or e[first] == e[second]:
        raise ValueError('Ramsey modes must have equal charge and distinct energy')
    if abs(q[first]+q[partner]) > 1e-12:
        raise ValueError('the supplied pair preparation must have zero total charge')
    enc = encoder_certificate()
    matrices, internal_preparation, read = prepare_internal(first, second, partner)
    native_input, source_preparation = prepare_source()
    read_word = ramsey_word(first, second)
    pointer_word = compiler.cnot(0, 1)
    pointer = compiler.execute(pointer_word, 2)
    if np.linalg.norm(pointer-np.asarray([[1, 0, 0, 0], [0, 0, 0, 1],
                                        [0, 0, 1, 0], [0, 1, 0, 0]]), 2) > ADMISSION:
        raise ArithmeticError('complete neutral occupation pointer changed')
    pointer_cost = compiler.resources(pointer_word*14)
    # Nine local factors implement the decoder (seven Fredkins, one rotation,
    # one scalar-corrected Z). Fourteen checked two-bit gates implement the
    # pointer copies. Tensoring identities does not increase operator error;
    # telescoping products gives (1+epsilon)^23-1 without assuming the computed
    # matrices are exactly unitary. This bounds only the two ideal contractions
    # replacing local source executions, not a hardware error or interval proof.
    contraction_delta = math.expm1(23*math.log1p(ADMISSION))
    gates = {j: phase_gate(e, q, t, angle, j, vacuum_energy=vacuum_energy) for j in (-1, 0, 1)}
    charge_unitary, charge_word, charge_certificate = charge_probe(q)
    charge_maps = (charge_unitary[:SIZE, :SIZE], charge_unitary[SIZE:, :SIZE])
    arrays, laws = {}, {}
    for law in ('C', 'E'):
        rows, total = [], 0.
        native_costs = {(row['coin'], row['raw']): row for row in cover._resource_rows(law, SEAM)}
        for row in cover.repair_branches(law, SEAM, native_input):
            raw, coin, current = row['record_code'], row['coin'], row['current']
            active = not row['abort']
            unitary, word = gates[current] if active else (np.eye(SIZE), [])
            evolved = np.stack([unitary@matrix@unitary.T for matrix in matrices])
            vectors, internal = decoder_factors(row['amplitude'], evolved)
            # Current+1 moves 0 across spectator1. The record therefore also
            # identifies which sorted endpoint carries the Ramsey particle.
            transported_slot = int(current == 1)
            internal = np.stack([read@matrix if transported_slot == 0 else matrix@read.T
                                 for matrix in internal])
            # The fresh neutral charge pointer starts in zero. Keep both
            # output amplitudes of its actually executed full256-column gate.
            charged_internal = np.stack([
                [operator@matrix if transported_slot == 0 else matrix@operator.T
                 for matrix in internal] for operator in charge_maps])
            probabilities = np.stack([[pointer_probabilities(v, m) for m in charged_internal]
                                      for v in vectors])
            mass = float(probabilities.sum())
            total += mass
            tag = f'{law}_{coin}_{raw}'
            arrays[tag+'_source_factors'] = np.stack(vectors)
            arrays[tag+'_internal_factors'] = charged_internal
            arrays[tag+'_all_pointer_probabilities'] = probabilities
            # A zero ideal branch remains in every array, including any
            # numerical padding/complement amplitude produced by the source.
            selected = float(probabilities[0, :, first, :].sum() if transported_slot == 0
                             else probabilities[0, :, :, first].sum())
            exact_row = native.branch(law, SEAM, coin, raw, -1)
            expected_weight = .5 if native.INDEX[SOURCE_PAIR] in exact_row.inputs else 0.
            effective = t*(1+current/4) if active else 0.
            expected = expected_weight*math.cos(effective*(e[second]-e[first])/2)**2
            if abs(selected-expected) > ADMISSION or abs(mass-expected_weight) > ADMISSION:
                raise ArithmeticError('complete colored Ramsey observation differs from its independent phase law')
            charge_expected = expected_weight*(1+math.cos(float(q[first])+math.pi/2))/2
            if abs(float(probabilities[:, 0].sum())-charge_expected) > ADMISSION:
                raise ArithmeticError('same-process neutral-pointer charge readout changed')
            phase_cost = compiler.resources(word)
            count = _sum_resources(source_preparation, compiler.resources(internal_preparation),
                                   native_costs[coin, raw], phase_cost, phase_cost,
                                   enc['decoder_resources'], compiler.resources(read_word),
                                   compiler.resources(charge_word), pointer_cost)
            rows.append({'coin': coin, 'raw': raw, 'abort': row['abort'],
                         'current': current, 'effective_duration': effective,
                         'internal_Ramsey_slot': transported_slot,
                         'ideal_native_weight': expected_weight,
                         'ideal_selected_joint_probability': expected,
                         'ideal_charge_pointer_zero_joint_probability': charge_expected,
                         'joint_probability_error_below': str(ADMISSION),
                         'array_prefix': tag, 'complete_internal_pointer_outcomes': 4*SIZE*SIZE,
                         'source_addresses_per_pointer_outcome': 512,
                         'resources': source._failure_account(count)})
        if abs(total-1) > ADMISSION:
            raise ArithmeticError('complete native-plus-colored-pointer instrument lost probability')
        laws[law] = rows
    report = {'schema': 'oph.colored_quasiparticle_dynamics.v1',
              'internal_dimension': len(e), 'padded_internal_dimension': SIZE,
              'source_pair': list(SOURCE_PAIR), 'native_seam': list(SEAM),
              'ramsey_indices': list(ramsey_indices), 'opposite_charge_partner': partner,
              'specified_energies': list(map(float, e)), 'specified_charges': list(map(float, q)),
              'specified_vacuum_energy': float(vacuum_energy),
              'spatial_carrier': {
                  'native_ports': 12, 'quasiparticle_modes_per_port': len(e),
                  'underlying_CAR_modes': 12*len(e),
                  'fixed_particle_number': 2, 'same_port_collision_states_in_this_packet': False,
                  'replicated_local_vacuum_scalar': 12*float(vacuum_energy),
                  'chosen_common_reference_scalar': float(vacuum_energy),
                  'added_reference_counterterm': -11*float(vacuum_energy)},
              'duration': t, 'electromagnetic_phase': angle,
              'scheduler_rule': 'effective_duration=duration*(1+native_current/4)',
              'generator_rule': 'U_j=exp(-i*t*(1+j/4)*(E+E_vac/2)+i*electromagnetic_phase*j*Q); U_j tensor U_j retains the single full vacuum scalar',
              'colored_repair': 'native crossing acts as -SWAP on the two internal fibers',
              'same_internal_unitary_on_each_particle': True,
              'encoder': enc, 'charge_probe': charge_certificate, 'laws': laws,
              'field_state_representation': 'exact finite sum of source/internal tensor factors; no Schmidt truncation',
              'pointer_fields': ['decoded_orientation', 'neutral_charge_pointer',
                                 'sorted_first_internal', 'sorted_second_internal'],
              'proper_code_failure_probability_exact': '0',
              'full_internal_preparation_source_execution_qubits': 15,
              'local_pointer_source_execution_input_columns': 4,
              'ideal_contracted_source_factors': {
                  'decoder_factors': 9, 'internal_pointer_copies': 14,
                  'local_operator_error_below': str(ADMISSION),
                  'normalized_input_operator_error_below': contraction_delta,
                  'normalized_input_probability_L1_error_below':
                      2*contraction_delta+contraction_delta**2,
                  'normalized_input_probability_TV_error_below':
                      contraction_delta+contraction_delta**2/2,
                  'bound_kind': 'numerically admitted local operators plus analytic telescoping; not interval arithmetic'},
              'source_instrument_arbitrary_reference_extension_used': True,
              'static_background_quasiparticle_number_two_required': True,
              'record_infers_energy_difference_not_absolute_rest_energy': True,
              'computed_absolute_energies_are_supplied_action_data': True,
              'native_positions_identified_with_physical_spacetime': False,
              'dynamical_Higgs_or_pair_creating_interactions_on_this_sector': False,
              'physical_species_identification_complete': False,
              'floating_execution_is_not_interval_or_hardware_certificate': True}
    return report, arrays
