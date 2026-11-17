"""One complete finite repair/propagation/Higgs/clock source experiment.

These are selected four-mode source factors of the spatial action, with a
fixed opposite-charge spectator and external boundary links omitted. This
is not an invariant subspace or an execution of the full Fock1080 dynamics.
All source and clock records, collision holds, padding and reference columns
are retained. Stored factors reconstruct every complex output amplitude.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from tools import strict_json
from ..common_source.chiral_candidate import gaussian_source
from ..defect_field_clock import clock
from ..defect_field_interactions import interaction
from ..defect_field_benchmark import quasiparticles as qp
from ..defect_protection import packet as protection_packet
from ..observer.source_scalar import compiler
from . import full_domain
from .dynamics_packet import compare_arrays, compare_json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DURATION = F(1)
PHASE = math.pi/2
ADMISSION = 2e-7
ADDITIVE_RESOURCE_KEYS = (
    'native_pulses', 'single_code_pulses', 'pair_code_pulses',
    'identity_phase_pulses', 'proper_code_load_events', 'proper_code_export_events',
    'proper_code_service_events')


def ladders(modes=4):
    result = []
    for mode in range(modes):
        matrix = np.zeros((1 << modes, 1 << modes), complex)
        for mask in range(1 << modes):
            if mask >> mode & 1:
                matrix[mask ^ (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
        result.append(matrix)
    return tuple(result)


def charge():
    """Actual component(27,78) charges: +1,-1 at each spatial site."""
    return np.diag([sum((1 if i % 2 == 0 else -1)*((mask >> i) & 1)
                        for i in range(4)) for mask in range(16)])


def current(raw):
    if type(raw) is not int or not 0 <= raw < 16:
        raise ValueError('one complete raw native record required')
    return 1 if raw in (1, 5) else -1 if raw in (2, 6) else 0


def ideal_chart():
    """Independent Pauli construction of C's corrected QP-to-bare chart."""
    Z = np.diag([1., -1.])
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]], complex)
    data = interaction.bare_component_chart((27, 78))
    # This actual component has real positive pairing. Its Higgs creation
    # phase convention is checked in C's full calibrated channel certificate.
    if data['creation_phase_sign'] != 1:
        raise ArithmeticError('actual charged-lepton Higgs convention changed')
    original = qp.quasiparticles()
    if (original['pairing'][27, 78].imag != 0 or original['pairing'][27, 78].real <= 0
            or not np.allclose(np.diag(original['original_nambu_charge'])[[27, 78]], [1., -1.],
                               atol=1e-12, rtol=0)):
        raise ArithmeticError('the actual selected component no longer has the declared mass phase and charges')
    rotation = expm(1j*math.pi/8*(np.kron(X, Y)+np.kron(Y, X)))
    phase = expm(-1j*math.pi/8*(np.kron(Z, np.eye(2))+np.kron(np.eye(2), Z)))
    return rotation@phase


def routed_hopping_word(coefficient, duration, left, right):
    route = []
    swap = gaussian_source.local_gate('fswap').word
    for neighbor in range(right-1, left, -1):
        route += compiler.remap(swap, [neighbor, neighbor+1])
    gate = gaussian_source.local_gate('hopping', coefficient, duration)
    return route+compiler.remap(gate.word, [left, left+1])+compiler.inverse(route)


def field_gate(signed_current, duration=DURATION, phase=PHASE):
    effective = F(duration)*clock.drive(signed_current)
    chart = interaction.bare_component_chart((27, 78))
    R = ideal_chart()
    if np.linalg.norm(chart['source_unitary']-R, 2) > ADMISSION:
        raise ArithmeticError('source and independent corrected QP chart disagree')
    Rword = compiler.remap(chart['source_word'], [0, 1])+compiler.remap(chart['source_word'], [2, 3])
    hop = (routed_hopping_word(.5j, float(effective), 0, 2)
           +routed_hopping_word(-.5j, float(effective), 1, 3))
    pair = interaction.pair_word(chart['energy'], phase, float(effective))
    word = Rword+hop+compiler.inverse(Rword)+compiler.remap(pair, [0, 1])+compiler.remap(pair, [2, 3])
    actual = compiler.execute(word, 4)
    a = ladders()
    Hhop = sum((z*a[i].conj().T@a[j]+z.conjugate()*a[j].conj().T@a[i]
                for i, j, z in ((0, 2, .5j), (1, 3, -.5j))), np.zeros((16, 16), complex))
    pairH = interaction.pair_hamiltonian(chart['energy'], phase)
    Hpair = np.kron(np.eye(4), pairH)+np.kron(pairH, np.eye(4))
    jointR = np.kron(R, R)
    target = expm(-1j*float(effective)*Hpair)@jointR.conj().T@expm(-1j*float(effective)*Hhop)@jointR
    if np.linalg.norm(actual-target, 2) > ADMISSION:
        raise ArithmeticError('full16-column joint field source differs from independent CAR product')
    Q = charge()
    if np.linalg.norm(actual@Q-Q@actual, 2) > ADMISSION:
        raise ArithmeticError('joint field evolution changed actual electric charge')
    return actual, target, {'current': signed_current, 'effective_duration': str(effective),
        'source_input_columns': 16, 'operator_error_below': str(ADMISSION),
        'resources': compiler.resources(word), 'source_word': word}


def preparation():
    """Actual128-address state, including untouched spectator and reference."""
    target = np.zeros(128)
    # Active bits0..3; spectator(+,-) bits4,5; which-path reference bit6.
    target[1+(1 << 5)] = target[4+(1 << 5)+(1 << 6)] = 1/math.sqrt(2)
    word = compiler.compile_real_state(target)
    blank = np.zeros(128, complex)
    blank[0] = 1
    actual = compiler.execute(word, 7, blank)
    if np.linalg.norm(actual-target) > ADMISSION:
        raise ArithmeticError('source joint field/reference preparation changed')
    matrix = actual.reshape(2, 4, 16).transpose(2, 1, 0).reshape(16, 8)
    reference_charge = np.diag(np.repeat([0., 1., -1., 0.], 2))
    if np.linalg.norm(charge()@matrix+matrix@reference_charge) > ADMISSION:
        raise ArithmeticError('prepared charged active field and spectator are not neutral')
    return matrix, {'source_qubits': 7, 'all_preparation_addresses': 128,
                    'resources': compiler.resources(word)}


def repair(law, policy='coherent'):
    if policy == 'coherent':
        from . import coherent_extension
        implementation = coherent_extension
    elif policy == 'flagged':
        implementation = full_domain
    else:
        raise ValueError('coherent or flagged collision policy required')
    actual, target, resources = implementation.source_isometry(2, 2, law)
    # Physical source output index: field +16*coin +32*raw.
    K = actual.reshape(16, 2, 16, 16).transpose(1, 0, 2, 3).reshape(32, 16, 16)
    ideal = target.reshape(16, 2, 16, 16).transpose(1, 0, 2, 3).reshape(32, 16, 16)
    if np.linalg.norm(sum(k.conj().T@k for k in K)-np.eye(16), 2) > ADMISSION:
        raise ArithmeticError('complete source repair effects changed')
    if max(np.linalg.norm(k@charge()-charge()@k, 2) for k in K) > ADMISSION:
        raise ArithmeticError('source repair changed actual electric charge')
    return K, ideal, resources


def field_output(arrays, law, first, second):
    K = arrays[law+'_repair_operators']
    U = arrays['field_unitaries']
    return (U[current(second % 16)+1]@K[second]
            @U[current(first % 16)+1]@K[first]@arrays['field_preparation'])


def build():
    initial, prep_record = preparation()
    gates = [field_gate(j) for j in (-1, 0, 1)]
    arrays = {'field_preparation': initial,
              'field_unitaries': np.stack([row[0] for row in gates]),
              'ideal_field_unitaries': np.stack([row[1] for row in gates])}
    sequences = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    clock_rows, clocks = [], []
    for first, second in sequences:
        run = clock.history((first, second), DURATION)
        # The exact rational inverse consumes a declared portable readout grid.
        # Two rounded cells enter each binary marginal, adding at most one
        # grid unit to its error. Quantum outputs themselves remain untouched.
        grid = F(1, 10**9)
        rounded = {key: F(round(value*10**9), 10**9)
                   for key, value in run['probabilities'].items()}
        inferred = clock.infer_time(rounded, probability_error=clock.SERVICE_ERROR+grid)
        interval = list(map(F, inferred['effective_time_GeV_inverse_interval']))
        expected = DURATION*(clock.drive(first)+clock.drive(second))
        if not interval[0] <= expected <= interval[1]:
            raise ArithmeticError('record-only clock inverse misses common field action time')
        clocks.append(run)
        clock_rows.append({'currents': [first, second], 'inferred_time': inferred,
            'inverse_probability_grid': str(grid),
            'rounded_inverse_inputs': [{'outcome': list(key), 'probability': str(value)}
                                       for key, value in sorted(rounded.items())],
            'specified_time_for_comparison_only': str(expected), 'source_cost': run['source_cost'],
            'full_input_columns': 4, 'purifier_columns': 4})
    outcomes = [(x, y) for x in (0, 1) for y in (0, 1)]
    arrays['clock_outputs'] = np.array([[run['outputs'][key] for key in outcomes] for run in clocks])
    arrays['clock_operators'] = np.array([[run['work_kraus'][key] for key in outcomes] for run in clocks])
    arrays['clock_preparation'] = clocks[0]['prepared_clock_and_purifiers']
    clock_effects = [sum(k.conj().T@k for k in row) for row in arrays['clock_operators']]
    clock_references = [sum(c.conj().T@c for c in row) for row in arrays['clock_outputs']]
    clock_probabilities = np.sum(abs(arrays['clock_outputs'])**2, axis=(2, 3))
    initial_reference = np.kron(initial.conj().T@initial,
        arrays['clock_preparation'].conj().T@arrays['clock_preparation'])
    report_laws = {}
    collision = np.diag([int(mask & 3 == 3 or mask >> 2 == 3) for mask in range(16)])
    for policy, law in ((policy, law) for policy in ('coherent', 'flagged') for law in ('C', 'E')):
        K, idealK, repair_resources = repair(law, policy)
        key = policy+'_'+law
        arrays[key+'_repair_operators'], arrays[key+'_ideal_repair_operators'] = K, idealK
        probabilities = np.zeros((32, 32, 4))
        field_effects = np.zeros((9, 16, 16), complex)
        final_reference = np.zeros((32, 32), complex)
        number_distribution = np.zeros(5)
        final_field = np.zeros((16, 16), complex)
        first_collision, second_collision = 0., 0.
        created_collision = 0.
        source_operator_error_squared = 0.
        for first in range(32):
            j1 = current(first % 16)
            A1 = arrays['field_unitaries'][j1+1]@K[first]
            first_state = A1@initial
            created_collision += float(np.linalg.norm(collision@first_state)**2)
            if first % 16 == full_domain.COLLISION_RAW:
                first_collision += float(np.linalg.norm(first_state)**2)
            for second in range(32):
                j2 = current(second % 16)
                sequence = sequences.index((j1, j2))
                A = arrays['field_unitaries'][j2+1]@K[second]@A1
                idealA = (arrays['ideal_field_unitaries'][j2+1]@idealK[second]
                          @arrays['ideal_field_unitaries'][j1+1]@idealK[first])
                source_operator_error_squared += float(np.linalg.norm(A-idealA, 'fro')**2)
                state = A@initial
                mass = float(np.linalg.norm(state)**2)
                probabilities[first, second] = mass*clock_probabilities[sequence]
                field_effects[sequence] += A.conj().T@A
                final_reference += np.kron(state.conj().T@state, clock_references[sequence])
                final_field += state@state.conj().T
                for mask in range(16):
                    number_distribution[mask.bit_count()] += float(np.vdot(state[mask], state[mask]).real)
                if second % 16 == full_domain.COLLISION_RAW:
                    second_collision += mass
        joint_effect = sum(np.kron(field_effects[i], clock_effects[i]) for i in range(9))
        if np.linalg.norm(joint_effect-np.eye(64), 2) > ADMISSION:
            raise ArithmeticError('two-event field-and-clock operator effect is not complete')
        if np.linalg.norm(final_reference-initial_reference, 2) > ADMISSION:
            raise ArithmeticError('unread complete process changed untouched reference marginal')
        if abs(float(probabilities.sum())-1) > ADMISSION:
            raise ArithmeticError('complete two-event joint outcome tree lost mass')
        if math.sqrt(source_operator_error_squared) > ADMISSION:
            raise ArithmeticError('complete field history isometry differs from its ideal CAR source law')
        if first_collision > 1e-14 or created_collision <= 1e-8:
            raise ArithmeticError('created collision witness absent or present before evolution')
        if policy == 'coherent' and second_collision > 1e-14:
            raise ArithmeticError('coherent policy invented a collision flag')
        arrays[key+'_all_joint_probabilities'] = probabilities
        arrays[key+'_joint_effect'] = joint_effect
        arrays[key+'_untouched_reference_marginal'] = final_reference
        arrays[key+'_final_active_number_distribution'] = number_distribution
        arrays[key+'_final_field_density'] = final_field
        report_laws[key] = {'policy': policy, 'law': law, 'events': 2, 'native_raw_slots_per_event': 32,
            'native_histories': 1024, 'clock_outcomes_per_history': 4,
            'total_complete_outcomes': 4096,
            'first_collision_probability': first_collision,
            'created_second_collision_probability': created_collision,
            'second_collision_flag_probability': second_collision,
            'final_active_number_distribution': number_distribution.tolist(),
            'joint_probability_sum': float(probabilities.sum()),
            'joint_operator_effect_error_below': str(ADMISSION),
            'complete_field_history_isometry_Frobenius_error_below': str(ADMISSION),
            'untouched_reference_error_below': str(ADMISSION),
            'repair_resources_per_event': repair_resources,
            'complete_resources_by_current_sequence': [
                {'currents': list(sequence), 'additive_resources': {
                    resource: prep_record['resources'][resource]+2*repair_resources[resource]
                    +gates[sequence[0]+1][2]['resources'][resource]
                    +gates[sequence[1]+1][2]['resources'][resource]
                    +clocks[i]['source_cost'][resource]
                    for resource in ADDITIVE_RESOURCE_KEYS},
                 'field_preparation_calls': 1, 'clock_preparation_calls': 1,
                 'fresh_native_coin_bits': 2, 'fresh_native_raw_record_bits': 8,
                 'fresh_clock_readout_bits': 2,
                 'all_native_and_proper_code_failure_slots_retained': True}
                for i, sequence in enumerate(sequences)]}
    comparison = {law: float(np.sum(abs(np.linalg.eigvalsh(
        arrays['coherent_'+law+'_final_field_density']-arrays['flagged_'+law+'_final_field_density'])))/2)
        for law in ('C', 'E')}
    # Direct-sum history outputs make the admitted field isometry norm bound
    # independent of history count. Controlled clock isometries use B's full
    # input-column bound, including all four clock outputs.
    clock_error = float(clock.SERVICE_ERROR)
    instrument_error = ADMISSION*(1+clock_error)+clock_error
    preparation_error = ADMISSION*(1+clock_error)+clock_error
    output_error = instrument_error+(1+instrument_error)*preparation_error
    return {'schema': 'oph.joint_repair_field_clock.v1',
        'internal_component': [27, 78], 'actual_internal_charges': [1, -1],
        'local_energy_GeV': interaction.bare_component_chart((27, 78))['energy'],
        'base_interval_GeV_inverse': str(DURATION), 'Higgs_phase': PHASE,
        'field_law': 'repair K_h; QP-to-bare R; routed hopping; inverse R; local Higgs gates; repeat; read clock',
        'clock_drive': 'both hands and every field factor use tau*(1+j/4), once per event',
        'field_preparation': prep_record, 'field_gates': [row[2] for row in gates],
        'clock_sequences': clock_rows, 'laws': report_laws,
        'primary_collision_policy': 'coherent existing-raw0 hold',
        'flagged_collision_policy_is_control': True,
        'unread_final_field_trace_distance_between_policies': comparison,
        'policy_comparison_scope': 'C can preserve N1/N3 coherence through raw0; E already separates binary odd occupancy into raw4/5/6 on this odd-parity witness',
        'field_output_reconstruction': 'U_j2 K_h2 U_j1 K_h1 field_preparation',
        'joint_output_reconstruction': 'field_output tensor clock_outputs[sequence,outcome]; spectator/reference columns retained',
        'active_Fock_dimension': 16, 'all_active_particle_number_sectors_retained': 5,
        'untouched_spectator_Fock_dimension': 4, 'which_path_reference_dimension': 2,
        'clock_work_dimension': 4, 'clock_purifier_dimension': 4,
        'complete_operator_check_dimension': 64,
        'initial_total_particle_number_ideal': 2,
        'initial_and_evolved_total_electric_charge_ideal': 0,
        'proper_code_failure_slots_counted_and_ideal_zero': True,
        'complete_joint_instrument_operator_error_below': instrument_error,
        'prepared_complete_output_vector_error_below': output_error,
        'all_outcome_total_variation_error_below': output_error+output_error**2/2,
        'bound_scope': 'analytic composition conditional on full-column numerical source admissions, not outward-rounded certification',
        'source_factor_experiment_not_full1080_mode_execution': True,
        'full_lattice_boundary_links_omitted_not_an_invariant_subspace_claim': True,
        'exact_ordered_factors_not_exponential_of_summed_action': True,
        'all_records_and_small_floating_padding_outputs_retained': True,
        'full_physical_particle_identification_claimed': False,
        'numerical_admission_not_outward_or_hardware_certificate': True}, arrays


def bindings():
    values = protection_packet.source_bindings()
    paths = [HERE/'joint.py', HERE/'test_joint.py', HERE/'JOINT.md',
             Path(full_domain.__file__), HERE/'coherent_extension.py',
             Path(interaction.__file__), Path(clock.__file__), Path(qp.__file__),
             Path(gaussian_source.__file__), ROOT/'model/common_source/calibrated_matter.py',
             ROOT/'model/common_source/calibrated_matter_parameters.json']
    for path in paths:
        values[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(values.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    pins = bindings()
    report, arrays = build()
    report['source_sha256'] = pins
    if bindings() != pins:
        raise RuntimeError('joint source changed during execution')
    output, archive = HERE/'joint_results.json', HERE/'joint_results.npz'
    if args.write:
        output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
        np.savez_compressed(archive, **arrays)
    elif args.check:
        compare_json(strict_json.load(output), report)
        with np.load(archive, allow_pickle=False) as saved:
            compare_arrays({key: saved[key] for key in saved.files}, arrays)
    print('Verified complete two-event joint repair, field and clock source experiment.')


if __name__ == '__main__':
    main()
