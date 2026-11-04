"""Execute a fixed source interaction, readout, scheduler and native repair."""
from collections import Counter
from fractions import Fraction as F
import hashlib
import math

import numpy as np

from model.common_source import scalar_coupling
from model.common_source.slow_geometry import source as carrier
from model.observer import native
from model.observer.source_scalar import compiler
from model.observer.source_scalar.experiment import sources as quantum_sources
from model.relative_observer_geometry import clock_reconstruction

SEAMS = ((0, 1), (0, 5))
ANGLE = math.atan2(3, 4)


def selected_seam(bit, word):
    """Eight fair scheduler bits; every word, including idle, is retained."""
    if type(bit) is not int or bit not in (0, 1):
        raise ValueError('a committed binary matter record is required')
    if type(word) is not int or not 0 <= word < 256:
        raise ValueError('an eight-bit scheduler word is required')
    if word >= 240:
        return None
    if bit == 0 and word == 8:
        return SEAMS[0]
    if bit == 1 and word == 0:
        return SEAMS[1]
    return carrier().seams()[word//8]


def record_probabilities(energy):
    if type(energy) is not F or not 0 <= energy <= 1:
        raise ValueError('an exact electric-energy expectation in [0,1] is required')
    p1 = F(1, 50)+F(24, 25)*energy
    return 1-p1, p1


def scheduler_weights(energy, *, conditioned_active=False):
    if type(conditioned_active) is not bool:
        raise ValueError('the active-clock convention must be explicit')
    probability = record_probabilities(energy)
    weights = {edge: F(0) for edge in carrier().seams()}
    weights[None] = F(0)
    for bit in (0, 1):
        for word in range(256):
            weights[selected_seam(bit, word)] += probability[bit]/256
    if conditioned_active:
        return {edge: p/(1-weights[None]) for edge, p in weights.items() if edge is not None}
    return weights


def _resources(word):
    result = compiler.resources(word)
    result.pop('zero_angle_requests', None)
    return result


def quantum_execution():
    """Production contractions use emitted native words on every input column."""
    quantum_sources()
    interaction = scalar_coupling.scalar_subwords(2, 4*ANGLE, coupling=F(1, 4))['interaction']
    preparation = [compiler.pulse('X', math.pi/4, [2])]
    u = compiler.execute(preparation+interaction, 4)
    pointer_word = compiler.cnot(2, 4)
    pointer_full = compiler.execute(pointer_word, 5)
    pointers = tuple(pointer_full[bit*16:(bit+1)*16, :16] for bit in (0, 1))
    inputs = np.eye(16, dtype=complex)[:, [0, 1]]
    effects = []
    max_error = 0.
    for bit in (0, 1):
        actual = pointers[bit]@u@inputs
        target = np.zeros((16, 2), complex)
        for charge in (0, 1):
            z = 1-2*charge
            target[charge+4*bit, charge] = ((4+3*z) if bit == 0 else 1j*(4-3*z))/math.sqrt(50)
        max_error = max(max_error, float(np.linalg.norm(actual-target, 2)))
        effects.append(actual.conj().T@actual)
    max_error = max(max_error, float(np.linalg.norm(sum(effects)-np.eye(2), 2)))
    if max_error > 1e-9:
        raise ArithmeticError('emitted source interaction/pointer failed its exact target control')
    coin_word = [compiler.pulse('X', math.pi/4, [0])]+compiler.cnot(0, 1)
    coin = compiler.execute(coin_word, 2)
    target_coin = np.array([1., 0., 0., 1j])/math.sqrt(2)
    if np.linalg.norm(coin[:, 0]-target_coin, 2) > 1e-9:
        raise ArithmeticError('native fair-coin word failed its retained two-bit target')
    charged_preparation = [compiler.pulse('X', math.pi/2, [0])]
    charged = compiler.execute(charged_preparation, 4)[:, 0]
    if abs(charged[1]-1j) > 1e-9 or np.linalg.norm(np.delete(charged, 1)) > 1e-9:
        raise ArithmeticError('charged endpoint preparation failed')
    return {
        'register': 'two Gauss-sector charged bits and two dimension-two scalar sites',
        'fixed_interaction': 'H_int=(1/4) X_scalar0 Z_charge0',
        'fixed_interaction_time': '4 atan(3/4)',
        'source_subword': 'model.common_source.scalar_coupling.scalar_subwords(2,4 atan(3/4),1/4)[interaction]',
        'source_control_receipt_pins_verified': True,
        'scalar_preparation': '(|0>+i|1>)/sqrt(2)',
        'pointer': 'executed native occupation pointer on scalar bit2',
        'effective_charge_Kraus': ['(4 I+3 Z)/sqrt(50)', 'i(4 I-3 Z)/sqrt(50)'],
        'terminal_work_columns_times_sqrt50': [
            {'pointer_and_scalar_bit': 0, 'charge_diagonal': [[7, 0], [1, 0]]},
            {'pointer_and_scalar_bit': 1, 'charge_diagonal': [[0, 1], [0, 7]]}],
        'record1_effect_diagonal': ['1/50', '49/50'],
        'electric_energy': 'E^2=n_charge0 on the supplied charge1=0 sector',
        'electric_energy_unconditional_expectation_preserved': True,
        'complete_pointer_outputs_and_terminal_work_retained': True,
        'all_input_columns_executed': True,
        'floating_control_error_threshold': '1/1000000000',
        'floating_control_is_the_exact_algebraic_proof': False,
        'polar_stabilization_used_by_this_numerical_control': False,
        'physical_or_rounded_control_complement_failure_instrument_certified': False,
        'exact_ideal_proper_code_has_zero_complement_failure': True,
        'interaction_and_scalar_preparation_resources': _resources(preparation+interaction),
        'pointer_resources': _resources(pointer_word),
        'one_fair_coin_word_resources': _resources(coin_word),
        'fresh_fair_coin_words_per_attempt': 9,
        'charged_endpoint_preparation_resources': _resources(charged_preparation),
        'physical_fresh_state_and_storage_costs_derived': False,
    }


def native_execution():
    """Enumerate every native work input and complete controller history."""
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    digest = hashlib.sha256()
    counts = [Counter(selected_seam(bit, word) for word in range(256)) for bit in (0, 1)]
    conditional_means = {}
    for state in states:
        state_means = []
        for bit in (0, 1):
            totals = [0]*12
            for word in range(256):
                seam = selected_seam(bit, word)
                for coin in (0, 1):
                    if seam is None:
                        output, current = state, 0
                    else:
                        output, record = native.native(state, seam, coin)
                        if record[0] != coin or record[1] != state[seam[0]]-output[seam[0]]:
                            raise ArithmeticError('native committed current or placement coin changed')
                        current = record[1]
                        expected = list(state)
                        expected[seam[0]] -= current
                        expected[seam[1]] += current
                        if tuple(expected) != output:
                            raise ArithmeticError('retained current does not reconstruct the load increment')
                    if output not in index:
                        raise ArithmeticError('native total-two work state left its carrier')
                    for i, value in enumerate(output):
                        totals[i] += value
                    row = (index[state], bit, word, coin, seam, current, index[output])
                    digest.update((repr(row)+'\n').encode())
            state_means.append(tuple(F(value, 512) for value in totals))
        conditional_means[state] = state_means
    from .geometry import laplacian
    endpoints = []
    for energy in (F(0), F(1, 2), F(1)):
        probs = record_probabilities(energy)
        theta = F(3, 25)*(1-2*energy)
        matrix = laplacian(theta)
        for state, means in conditional_means.items():
            actual = [sum(probs[b]*means[b][i] for b in (0, 1)) for i in range(12)]
            expected = [state[i]-sum(matrix[i][j]*state[j] for j in range(12))/64
                        for i in range(12)]
            if actual != expected:
                raise ArithmeticError('actual source histories disagree with the spectral response matrix')
        weights = scheduler_weights(energy, conditioned_active=True)
        clock_rows = clock_reconstruction.clock_records(weights)
        recovered, _ = clock_reconstruction.recover(clock_rows)
        if recovered != weights:
            raise ArithmeticError('local pair-clock inference lost the actual seam law')
        endpoints.append({'electric_energy': str(energy), 'theta': str(theta),
                          'conditional_active_weights': {str(e): str(p) for e, p in weights.items()},
                          'all_132_pair_clock_laws_recover_same_weights': True})
    x = (2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    y = (0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    correlated = [sum(record_probabilities(energy)[b]*conditional_means[state][b][i]/2
                      for energy, state in ((F(0), x), (F(1), y)) for b in (0, 1))
                  for i in range(12)]
    average = [(a+b)/2 for a, b in zip(map(F, x), map(F, y))]
    base = laplacian()
    product = [average[i]-sum(base[i][j]*average[j] for j in range(12))/64
               for i in range(12)]
    difference = [a-b for a, b in zip(correlated, product)]
    if max(map(abs, difference)) != F(3, 800):
        raise ArithmeticError('correlated probe control failed to distinguish the marginal mean law')
    return {
        'native_work_states': len(states),
        'complete_histories_per_native_work_state': 1024,
        'complete_native_histories_checked': len(states)*1024,
        'complete_history_sha256': digest.hexdigest(),
        'conditional_scheduler_slot_counts': [{str(e): count for e, count in counter.items()}
                                              for counter in counts],
        'native_idle_probability': '1/16',
        'all_thirty_seams_have_positive_probability': True,
        'zero_current_and_idle_records_retained': True,
        'all_initial_load_mean_identities_checked_exactly': True,
        'observer_clock_comparisons': endpoints,
        'fresh_probe_mean_response': 'M=I-[L+theta(D_(0,1)-D_(0,5))]/64',
        'energy_to_theta': 'theta=(3/25)(1-2<E^2>)',
        'kinetic_form_energy_derivative': '-3/800 * [(f0-f1)^2-(f0-f5)^2]',
        'correlated_probe_countercontrol': {
            'joint_input': 'equal mixture of (charge0,2 at port0) and (charge1,2 at port1)',
            'same_native_and_charge_marginals_as_product_control': True,
            'actual_minus_product_mean': [str(x) for x in difference],
            'maximum_load_mean_difference': '3/800'},
    }
