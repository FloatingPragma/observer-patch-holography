"""Independent contractions of the frozen particle measurement construction.

Producer run/verify, CommonClock and operational channel helpers are not used.
The declared source and literal M6/native primitives are shared dependencies.
"""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import copy
import hashlib
import json
import math
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
NEXT = ROOT / 'openwave-tuning/next'
PART = HERE.parent / 'particle'
CLOCK = HERE.parent / 'clock'
OP = NEXT / 'ROW1_NATIVE/routes/quantum_operational'
RER = ROOT / 'reverse-engineering-reality'
sys.path[:0] = [str(RER), str(NEXT / 'C10')]
from joint_source import JointNativeSource
from model.observer.source_scalar import compiler
from model.observer import native
from model.observer.source_discrimination import circuit
from model.native_interface import pointer
from model.native_clock import configuration

FROZEN = {
    PART / 'PARTICLE_CONTRACT.json': '67c406d20f796193d2a88f3c8e30f463fbda1cf7ce8a93cf75f353bf8a8f6cf4',
    PART / 'RESULT.json': '9bcce002927e744b8cf2e37d78f2f26ceb0e442a36a2058a6d78541291c171b1',
    OP / 'OPERATIONAL_QUANTUM_REVIEW.json': '8e7b32459b12e1dc099e5233b9ec1f230087073d3ea17993e5a283df6c050d4e',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, tolerance=2e-12, label='array'):
    error = float(np.max(abs(np.asarray(actual) - np.asarray(expected))))
    require(error < tolerance, f'{label}: {error}')
    return error


def tr(rho):
    return float(np.trace(rho).real)


def pins(document, base=ROOT):
    for path, expected in document['source_sha256'].items():
        require(digest(base / path) == expected, 'changed dependency: ' + path)


def square_apply(gate, targets, bits, columns):
    """Explicit basis-index lift, independent of compiler tensor reshapes."""
    targets = list(targets)
    lifted = np.zeros((1 << bits, 1 << bits), complex)
    for source in range(1 << bits):
        local = sum(((source >> target) & 1) << (len(targets)-1-j)
                    for j, target in enumerate(targets))
        for dest in range(len(gate)):
            full = source
            for j, target in enumerate(targets):
                full = (full & ~(1 << target)) | (((dest >> (len(targets)-1-j)) & 1) << target)
            lifted[full, source] = gate[dest, local]
    return lifted @ columns


def rectangular_apply(gate, targets, bits, columns):
    """Retain a rectangular output and every untouched input coordinate."""
    targets = list(targets)
    spectators = [j for j in range(bits-1, -1, -1) if j not in targets]
    output = np.zeros((len(gate) * (1 << len(spectators)), columns.shape[1]), complex)
    for source in range(1 << bits):
        local = sum(((source >> target) & 1) << (len(targets)-1-j)
                    for j, target in enumerate(targets))
        rest = sum(((source >> target) & 1) << (len(spectators)-1-j)
                   for j, target in enumerate(spectators))
        for dest in range(len(gate)):
            output[(dest << len(spectators)) + rest] += gate[dest, local] * columns[source]
    return output


def replay(word, bits, initial):
    state = initial.copy()
    aborts = []
    unitary_error = 0.
    for event, row in enumerate(word):
        full = compiler.native_pulse(row['kind'], float(Fraction(str(row['angle']))))
        unitary_error = max(unitary_error, close(full.conj().T @ full, np.eye(6), label='literal pulse unitarity'))
        d = 1 << len(row['targets'])
        for complement, outrow in enumerate(full[d:, :d]):
            aborts.append((event, complement, rectangular_apply(outrow[None, :], row['targets'], bits, state)))
        state = square_apply(full[:d, :d], row['targets'], bits, state)
    effect = sum(k.conj().T @ k for _, _, k in aborts)
    close(state.conj().T @ state + effect, initial.conj().T @ initial, label='complete source word')
    return state, aborts, unitary_error


def native_replay(operational, archive, result):
    word = operational['decision_oracles']['three_site']['admitted_word']['word']
    initial = np.zeros((32, 8), complex)
    initial[:8] = np.eye(8)
    oracle, aborts, unitary_error = replay(word, 5, initial)
    expected = np.zeros_like(oracle)
    for bits in range(8):
        expected[bits + 8 * bits.bit_count(), bits] = 1
    oracle_error = close(oracle, expected, 2e-11, 'population count oracle')
    pword = operational['native_pointer_source_word']
    pinit = np.eye(8, dtype=complex)[:, [pointer.READY_CODE, pointer.READY_CODE + 4]]
    coupled, paborts, _ = replay(pword, 3, pinit)
    require(len(word) == 114 and len(pword) == 6, 'source pulse count')
    pointer_records = {}
    for law in ('C', 'E'):
        grouped, metadata = defaultdict(dict), {}
        for z in range(4):
            carrier = circuit.decode(z)
            for coin in (0, 1):
                output, (_, current) = native.native(carrier, pointer.SEAM, coin)
                label = current if law == 'C' else carrier[0] - carrier[1]
                key = coin, label
                grouped[key].setdefault(output, np.zeros((2, 2), complex))
                grouped[key][output] += coupled[[z, z + 4]] / math.sqrt(2)
                metadata[key] = dict(law=law, coin=coin, label=label, current=current,
                                     decoded_pointer_bit=int(2*current + 2*coin - 1 < 0))
        pointer_records[law] = [(metadata[key], list(grouped[key].values())) for key in sorted(grouped)]
    channels = {}
    matrix_unit_error = 0.
    for law, records in pointer_records.items():
        branches = {q: [] for q in range(4)}
        fine = []
        for meta1, ks1 in records:
            for meta2, ks2 in records:
                q = meta1['decoded_pointer_bit'] + 2*meta2['decoded_pointer_bit']
                indices = []
                for k1 in ks1:
                    for k2 in ks2:
                        final = square_apply(k2, [4], 5, square_apply(k1, [3], 5, oracle))
                        for output in range(4):
                            indices.append([q, len(branches[q])])
                            branches[q].append(final[8*output:8*(output+1)])
                fine.append(dict(first=meta1, second=meta2, coarse_sum=q, kraus_indices=indices))
        require(fine == operational['native_pointer_fine_records'][law], 'fine pointer records changed')
        if law == 'C':
            require(fine == result['fine_native_readout_records'], 'particle fine records changed')
        channels[law] = {}
        for q, ks in branches.items():
            channels[law][q] = sum(np.einsum('ai,bj->abij', k, k.conj()) for k in ks)
            p = np.diag([int(i.bit_count() == q) for i in range(8)])
            matrix_unit_error = max(matrix_unit_error, close(channels[law][q], np.einsum('ai,bj->abij', p, p), 2e-11, 'coarse matrix units'))
            if law == 'C':
                close(np.asarray(ks), archive[f'C_Q{q}_local_kraus'], label='actual pointer Kraus archive')
    families = {'oracle': [(dict(event=e, complement=c), k) for e, c, k in aborts],
                'pointer_first': [], 'pointer_second': []}
    for event, complement, k in paborts:
        families['pointer_first'].append((dict(event=event, complement=complement), rectangular_apply(k, [3], 5, oracle)))
        for ir, (meta, ks) in enumerate(pointer_records['C']):
            for hidden, first in enumerate(ks):
                prefix = square_apply(first, [3], 5, oracle)
                families['pointer_second'].append((dict(event=event, complement=complement,
                    first_native_record=ir, first_record=meta, first_hidden_output=hidden), rectangular_apply(k, [4], 5, prefix)))
    effects, abort_error = {}, 0.
    for name, family in families.items():
        rows = result['export_failure_kraus_families'][name]
        require(len(family) == len(rows), 'abort slot lost')
        for (meta, k), row in zip(family, rows):
            require(all(row[key] == value for key, value in meta.items()), 'abort metadata changed')
            require(len(k) == row['retained_local_output_dimension'], 'abort output dimension')
            # Export amplitudes are of floating-point size. Absolute 1e-12
            # tolerances would accept deleting every matrix in this family.
            # Products of roundoff leakage can be below 1e-30; changing
            # contraction order changes their relative error by order one.
            # The 1e-16 scale floor gives a 2e-28 absolute bound, sixteen
            # orders tighter than the ordinary density-state tolerance.
            relative = float(np.max(abs(k - archive[row['local_kraus_key']])) / max(1e-16, np.max(abs(k))))
            require(relative < 2e-12, f'abort Kraus changed: {name} {meta}, scaled error={relative}')
            abort_error = max(abort_error, relative)
        effects[name] = sum(k.conj().T @ k for _, k in family)
    complete = sum(np.einsum('aaij->ji', ch) for ch in channels['C'].values()) + sum(effects.values())
    complete_error = close(complete, np.eye(8), label='complete instrument effect')
    reset = Fraction(1)
    for entry in operational['finite_native_reset_contract']:
        cert = entry['certificate']
        values = configuration.success_probabilities(cert['seed_scans'], cert['retries'], full=True,
                                                      mode='reset', target=tuple(cert['target']))
        require(list(map(str, values)) == entry['complete_input_success_probabilities'], 'native reset table changed')
        require(str(values[0]) == entry['success_probability'], 'initial reset input changed')
        reset *= values[0] ** 2
        require('timeout' in entry['all_fine_histories']['failure_retains'], 'timeout histories omitted')
    require(str(reset) == result['reset_success_exact_per_detector'], 'four-reset success law')
    require(str(1-reset) == result['reset_timeout_exact_per_detector'], 'reset timeout law')
    require(operational['finite_native_reset_contract'] == result['reset_contract'], 'reset history contract changed')
    return channels['C'], effects, reset, dict(source_pulses=114, pointer_source_pulses=12,
        complete_abort_maps=sum(map(len, families.values())), oracle_max_error=oracle_error,
        literal_pulse_unitary_error=unitary_error, coarse_matrix_unit_max_error=matrix_unit_error,
        complete_effect_max_error=complete_error, abort_map_max_scaled_error=abort_error,
        abort_amplitude_scale_floor=1e-16, abort_scaled_error_tolerance=2e-12,
        reset_timeout_probability=float(1-reset))


class IndependentClock:
    """Row-major Liouvillians assembled directly from the declared H and K."""
    def __init__(self, model):
        self.model = model
        self.sectors = []
        self.ops = [sparse.csr_matrix(k) / math.sqrt(3)
                    for seam in model.field.operators_by_law['MIX'] for _, k in seam]
        self.complete_effect_error = close(sum((k.conj().T @ k).toarray() for k in self.ops), np.eye(model.n), label='native clock effect')
        for ids, energies, vectors in model.spectral.values():
            h = sparse.csr_matrix(model.hamiltonian[np.ix_(ids, ids)])
            eye = sparse.eye(len(ids), format='csr')
            coherent = -1j*(sparse.kron(h, eye) - sparse.kron(eye, h.T))
            jump = sum(sparse.kron(k[ids][:, ids], k[ids][:, ids].conjugate(), format='csr') for k in self.ops)
            self.sectors.append((ids, energies, vectors, coherent.tocsr(), jump))

    def counts(self, batch, duration=.37, rate=.8):
        batch = np.asarray(batch)
        out = np.zeros((len(batch), 3, self.model.n, self.model.n), complex)
        for ids, energies, vectors, coherent, jump in self.sectors:
            block = batch[:, ids[:, None], ids]
            if not np.any(block):
                continue
            size = len(ids)**2
            x = block.reshape(len(batch), size).T
            killed = coherent - rate*sparse.eye(size, format='csr')
            # Independent row-major construction; N1 uses a triangular
            # generating-coefficient ODE, and tail is the full complement.
            gen = sparse.bmat([[killed, None], [rate*jump, killed]], format='csr')
            y = expm_multiply(duration*gen, np.vstack([x, np.zeros_like(x)]))
            total = expm_multiply(duration*(killed + rate*jump), x)
            states = (y[:size], y[size:], total-y[:size]-y[size:])
            for count, state in enumerate(states):
                out[:, count, ids[:, None], ids] = state.T.reshape(block.shape)
        return out


def exact_join(contract, clock, result):
    for key in ('source_id', 'source_dimension', 'clock_id', 'source_H_array_sha256_le_float64',
                'preparation_array_sha256_le_complex128'):
        require(contract[key] == clock[key], 'typed source/clock mismatch: '+key)
    require(contract['preparation_id'] == 'C10_global_reference_complete_fair_reservoir_v1', 'preparation ID')
    require(contract['reference_id'] == 'actual_C10_allnumber_global_interacting_ground', 'reference ID')
    expected = [key for row in clock['preparation_labels'] for key in row['complete_final_count_state_archive_keys']]
    require([row['clock_archive_key'] for row in result['clock_inputs']] == expected, 'complete input keys')
    rows = result['terminal_records']
    require(len(rows) == 576, 'terminal count')
    require({(r['clock_input_index'], r['first_regional_charge'], r['second_count_class'], r['second_regional_charge']) for r in rows}
            == {(i, a, n, b) for i in range(12) for a in range(4) for n in range(3) for b in range(4)}, 'complete terminal labels')
    require(len(result['first_detector_records']) == 12 and len(result['continued_clock_and_second_detector_records']) == 144, 'detector rows')
    for row in result['first_detector_records'] + result['continued_clock_and_second_detector_records']:
        require(set(row['stops']) == {'reset_timeout', 'oracle', 'pointer_first', 'pointer_second'}, 'failure flag deleted')
    require(contract['measurement_order']['global_Q_H_commute'] and not contract['measurement_order']['regional_Q_H_commute'], 'measurement order changed')
    forbidden = ('stable_particle_from_probability_conservation', 'asymptotic_lifetime_bound',
                 'physical_particle_species_identified', 'predicted_mass', 'physical_clock_calibrated', 'continuum_limit_transferred')
    require(not any(contract['claim_boundary'][k] for k in forbidden), 'unsupported physical award')


def run():
    for path, pin in FROZEN.items():
        require(digest(path) == pin, 'frozen input changed: '+str(path))
    contract, result = read(PART/'PARTICLE_CONTRACT.json'), read(PART/'RESULT.json')
    clock_contract, operational = read(CLOCK/'CLOCK_CONTRACT.json'), read(OP/'OPERATIONAL_QUANTUM_REVIEW.json')
    for doc in (contract, result, clock_contract, operational):
        pins(doc)
    require(digest(PART/contract['state_archive']) == contract['state_archive_sha256'], 'particle archive custody')
    require(digest(CLOCK/clock_contract['state_archive']) == clock_contract['state_archive_sha256'], 'clock archive custody')
    require(digest(OP/operational['archive']) == operational['archive_sha256'], 'operational archive custody')
    exact_join(contract, clock_contract, result)
    archive = np.load(PART/'particle_states.npz', allow_pickle=False)
    ca = np.load(CLOCK/'common_clock_states.npz', allow_pickle=False)
    source = np.load(NEXT/'C10/joint_source_states.npz', allow_pickle=False)
    model = JointNativeSource()
    h, numbers = model.hamiltonian, model.numbers
    close(h, source['hamiltonian'], label='regenerated source H')
    require(hashlib.sha256(np.asarray(h, dtype='<f8').tobytes()).hexdigest() == contract['source_H_array_sha256_le_float64'], 'actual H array ID')
    reference, ground, _ = model.reference()
    prep = [k @ reference @ k.T for _, k in model.source.reservoir_instrument()]
    close(sum(prep), ca['input_state'], label='actual preparation')
    require(hashlib.sha256(np.asarray(sum(prep), dtype='<c16').tobytes()).hexdigest() == contract['preparation_array_sha256_le_complex128'], 'actual preparation array ID')
    masks = np.array([mask for mask, _ in model.states])
    regional = sum((masks >> site) & 1 for site in (0, 1, 3))
    close(regional, archive['regional_occupation'], label='regional labels')
    # Reconstruct occupations from F and B before comparing to carrier masks.
    div = np.zeros((384, 6), int)
    edges = ((0,3), (1,3), (1,4), (2,4), (2,5), (0,5))
    records = np.array([f for _, f in model.states])
    for i, (a, b) in enumerate(edges):
        div[:, a] += records[:, i]
        div[:, b] -= records[:, i]
    decoded = (np.array([1,1,1,0,0,0]) - div - (3-numbers)[:, None]*np.eye(6, dtype=int)[0]) % 6
    close(decoded, ((masks[:, None] >> np.arange(6)) & 1), label='record-only Gauss decoder')
    indices = np.zeros((8, 48), int)
    seen = defaultdict(int)
    for i, (mask, _) in enumerate(model.states):
        local = sum(((mask >> site) & 1) << j for j, site in enumerate((0,1,3)))
        rest = sum(((mask >> site) & 1) << j for j, site in enumerate((2,4,5)))
        indices[local, 6*rest+seen[mask]] = i
        seen[mask] += 1
    require(set(seen.values()) == {6} and len(set(indices.flat)) == 384, 'coherent transfer isometry')
    close(indices, archive['coherent_code_indices'], label='coherent code indices')
    channels, effects, reset_exact, native_checks = native_replay(operational, archive, result)
    print('Native pulse, C/E matrix-unit, 594 abort-map, and exact reset checks passed.', flush=True)
    reset, timeout = float(reset_exact), float(1-reset_exact)
    flat = indices.ravel()

    def measure(rho):
        blocks = rho[np.ix_(flat, flat)].reshape(8, 48, 8, 48)
        reduced = np.trace(blocks, axis1=1, axis2=3)
        stops = {'reset_timeout': timeout*tr(rho)}
        stops.update({name: reset*tr(effect @ reduced) for name, effect in effects.items()})
        outputs = []
        for q in range(4):
            local = reset*np.einsum('abij,irjs->arbs', channels[q], blocks, optimize=True).reshape(384,384)
            out = np.zeros_like(rho, dtype=complex)
            out[np.ix_(flat, flat)] = local
            outputs.append(out)
        close(sum(map(tr, outputs)) + sum(stops.values()), tr(rho), label='detector full law')
        return outputs, stops

    def check_record(row, rho, outputs, stops):
        close(row['conclusive_probabilities'], list(map(tr, outputs)), label='conclusive record')
        check_stops(row['stops'], stops)

    def check_stops(expected, stops):
        for key, value in stops.items():
            # Zero physical inputs can acquire sub-roundoff numerical traces.
            # Scale by the relevant effect, with a 1e-15 input-weight floor.
            effect_scale = timeout if key == 'reset_timeout' else reset*float(np.linalg.norm(effects[key],2))
            tolerance = 2e-11*max(abs(value), effect_scale*1e-15, 1e-300)
            require(abs(expected[key]-value) <= tolerance, f'stop probability {key}: replay={value}, archived={expected[key]}, tolerance={tolerance}')

    independent_clock = IndependentClock(model)
    first_clock = independent_clock.counts(prep).reshape(12,384,384)
    clock_error = 0.
    first_states, first_stops = [], []
    work_first = 0.
    for i, (row, state) in enumerate(zip(result['clock_inputs'], first_clock)):
        clock_error = max(clock_error, close(state, ca[row['clock_archive_key']], label='first complete clock state'))
        close(tr(state), row['probability'], label='clock input probability')
        expected_meta = clock_contract['preparation_labels'][i//3]
        require(row['first_count_class'] == i%3 and all(row[k] == expected_meta[k] for k in ('reservoir_initial','reservoir_final')), 'clock input record labels')
        outputs, stops = measure(state)
        require(result['first_detector_records'][i]['clock_input_index'] == i, 'first detector input label')
        require(result['first_detector_records'][i]['input_state']['clock_archive_key'] == row['clock_archive_key'], 'first detector actual input')
        check_record(result['first_detector_records'][i], state, outputs, stops)
        first_states.extend(outputs)
        first_stops.append(stops)
        work_first += tr(h @ (sum(outputs)-reset*state))
    continued = independent_clock.counts(first_states).reshape(144,384,384)
    print('Independent first and second count-split states computed.', flush=True)
    offsets = archive['number_block_offsets']
    expected_offsets = np.cumsum([0]+[len(ids)**2 for ids,*_ in model.spectral.values()])
    close(offsets, expected_offsets, label='number compression offsets')
    packed = archive['continued_number_blocks']
    require(packed.shape == (144,33264), 'complete continuation archive shape')
    coefficients = archive['terminal_energy_coefficients']
    require(coefficients.shape == (576,348,2,2), 'complete energy coefficient shape')
    groups = result['complete_energy_groups']
    spectral_error = 0.
    bynumber = defaultdict(list)
    for gi, group in enumerate(groups):
        q, d = group['source_number'], group['degeneracy']
        ids = archive[group['position_indices_archive_key']]
        basis = archive[group['basis_archive_key']]
        require(group['global_charge'] == q-3 and basis.shape == (len(ids), d), 'energy group labels')
        close(ids, np.flatnonzero(numbers == q), label='spectral carrier indices')
        spectral_error = max(spectral_error, close(h[np.ix_(ids,ids)] @ basis, group['energy']*basis, 2e-11, 'actual energy eigenspace'))
        close(group['excitation_energy'], group['energy']-ground, label='actual reference subtraction')
        close(group['calibrated_excitation_GeV'], group['excitation_energy']*contract['calibration']['GeV_per_source_energy'], label='single display calibration')
        bynumber[q].append((gi, group, ids, basis))
    for q, rows in bynumber.items():
        full_basis = np.hstack([x[3] for x in rows])
        close(full_basis.conj().T @ full_basis, np.eye(len(full_basis)), 2e-11, 'complete grouped energy basis')
        require(len(rows) == sum(1 for i,e in enumerate(model.spectral[q][1]) if i == 0 or abs(e-model.spectral[q][1][i-1]) >= 1e-9), 'whole energy degeneracies')
    maximum_coefficient_error = continuation_error = 0.
    minimum_continuation_eigenvalue = 0.
    maximum_offnumber_entry = 0.
    terminal_probabilities = []
    second_stops = []
    work_second = work_clock = 0.
    actual_degenerate_offdiagonal = 0.
    min_coefficient_eigenvalue = float(np.linalg.eigvalsh(coefficients).min())
    require(min_coefficient_eigenvalue > -2e-12, 'negative compact branch state')
    for ci, state in enumerate(continued):
        row = result['continued_clock_and_second_detector_records'][ci]
        require((row['clock_input_index'], row['first_regional_charge'], row['second_count_class'], row['number_block_state_index']) == (ci//12,(ci//3)%4,ci%3,ci), 'continued record labels')
        for q, (ids, *_rest) in model.spectral.items():
            stored = packed[ci, offsets[q]:offsets[q+1]].reshape(len(ids),len(ids))
            continuation_error = max(continuation_error, close(state[np.ix_(ids,ids)], stored, label='all continued number blocks'))
            minimum_continuation_eigenvalue = min(minimum_continuation_eigenvalue, float(np.linalg.eigvalsh(stored).min()))
        maximum_offnumber_entry = max(maximum_offnumber_entry, float(np.max(abs(state*(numbers[:,None] != numbers[None,:])))))
        outputs, stops = measure(state)
        check_record(row, state, outputs, stops)
        second_stops.append(stops)
        work_second += tr(h @ (sum(outputs)-reset*state))
        for q2, terminal in enumerate(outputs):
            ti = 4*ci+q2
            terminal_probabilities.append(tr(terminal))
            termrow = result['terminal_records'][ti]
            require((termrow['terminal_index'], termrow['clock_input_index'], termrow['first_regional_charge'], termrow['second_count_class'], termrow['second_regional_charge'])
                    == (ti,ci//12,(ci//3)%4,ci%3,q2), 'ordered terminal labels')
            close(tr(terminal), result['terminal_records'][ti]['probability'], label='terminal probability')
            require(result['terminal_records'][ti]['number_block_state_index'] == ci, 'terminal source state index')
            for q, rows in bynumber.items():
                full_basis = np.hstack([x[3] for x in rows])
                ids = rows[0][2]
                transformed = full_basis.conj().T @ terminal[np.ix_(ids,ids)] @ full_basis
                start = 0
                for gi, group, _, basis in rows:
                    d = group['degeneracy']
                    block = transformed[start:start+d,start:start+d]
                    maximum_coefficient_error = max(maximum_coefficient_error, close(block, coefficients[ti,gi,:d,:d], label='all energy branch coefficients'))
                    if d == 1:
                        close(coefficients[ti,gi,1:,:], 0, label='unused coefficient padding')
                        close(coefficients[ti,gi,:,1:], 0, label='unused coefficient padding')
                    else:
                        actual_degenerate_offdiagonal = max(actual_degenerate_offdiagonal, float(abs(block[0,1])))
                    start += d
    for i, state in enumerate(first_states):
        work_clock += tr(h @ (sum(continued[3*i:3*i+3])-state))
        maximum_offnumber_entry = max(maximum_offnumber_entry, float(np.max(abs(state*(numbers[:,None] != numbers[None,:])))))
    require(minimum_continuation_eigenvalue > -2e-12, 'negative continued branch')
    require(maximum_offnumber_entry < 1e-13, 'number-block compression discards coherence')
    close(sum(continued), archive['continued_ensemble_state'], label='continued ensemble')
    all_stops = sum(sum(row.values()) for row in first_stops+second_stops)
    probabilities = np.trace(coefficients, axis1=2, axis2=3).real
    terminals = np.array(terminal_probabilities)
    close(probabilities.sum(axis=1), terminals, label='all terminal probabilities')
    close(terminals.reshape(12,4,3,4), archive['regional_transition_probabilities'], label='regional joint probabilities')
    close(terminals.sum()+all_stops, 1., label='complete unconditional measurement law')
    close(all_stops, result['unconditional']['all_apparatus_stop_probability'], label='all apparatus stops')
    transitions = terminals.reshape(12,4,3,4)
    close(terminals.sum(), result['unconditional']['primary_conclusive_probability'], label='primary conclusive marginal')
    close(transitions.sum(axis=(0,1,3)), result['unconditional']['second_clock_count_probabilities_after_two_detector_successes'], label='complete count marginal')
    changed = sum(transitions[:,a,:,b].sum() for a in range(4) for b in range(4) if a != b)
    close(changed, result['unconditional']['regional_count_changed_probability'], label='regional transport marginal')
    close(terminals.sum()-changed, result['unconditional']['regional_count_unchanged_probability'], label='unchanged regional marginal')
    charge_groups = [gi for gi,g in enumerate(groups) if g['global_charge'] == 1]
    charge_transition = probabilities[:,charge_groups].sum(axis=1).reshape(12,4,3,4).sum(axis=(0,2))
    charge_report = result['positive_charge_transport']
    close(charge_transition, charge_report['regional_transition_matrix'], label='charged regional transport matrix')
    close(charge_transition.sum(), charge_report['unconditional_probability'], label='charged unconditional probability')
    close(np.trace(charge_transition), charge_report['regional_count_unchanged_probability'], label='charged unchanged regional probability')
    close(charge_transition.sum()-np.trace(charge_transition), charge_report['regional_count_changed_probability'], label='charged regional transport probability')
    for name, value in [('detector1_success_weighted_source_energy_change', work_first),
                        ('detector2_success_weighted_source_energy_change', work_second),
                        ('continued_native_clock_source_energy_change', work_clock)]:
        close(value, result['unconditional'][name], label='conditional source energy ledger')
    selected = result['selected_positive_excitation']
    ti, gi = selected['terminal_index'], selected['energy_group_index']
    allowed = [i for i,g in enumerate(groups) if g['global_charge'] == 1 and g['excitation_energy'] > 1e-10]
    require(probabilities[ti,gi] == probabilities[:,allowed].max(), 'selected branch maximal unconditional probability')
    group = groups[gi]
    ids, basis = archive[group['position_indices_archive_key']], archive[group['basis_archive_key']]
    excitation = np.zeros((384,384), complex)
    d = group['degeneracy']
    excitation[np.ix_(ids,ids)] = basis @ coefficients[ti,gi,:d,:d] @ basis.conj().T
    selected_error = close(excitation, archive['selected_excitation_state'], label='selected actual subnormalized poststate')
    close(tr(excitation), selected['probability'], label='selected probability')
    close(h @ excitation, group['energy']*excitation, label='selected H eigenstate')
    close((numbers-4)[:,None]*excitation, 0, label='selected actual charge support')
    after = independent_clock.counts([excitation])[0]
    witness_probability = 0.
    for n, state in enumerate(after):
        row = selected['continued_count_and_regional_records'][n]
        require(row['count_class'] == n, 'selected count label')
        close(state, archive[row['pre_detector_state_key']], label='selected complete same-clock continuation')
        outputs, stops = measure(state)
        close(list(map(tr,outputs)), row['conclusive_regional_charge_probabilities'], label='selected regional poststates')
        for q, output in enumerate(outputs):
            close(output, archive[f"{row['pre_detector_state_key']}_Q{q}"], label='selected detector archived state')
        check_stops(row['stops'],stops)
        witness_probability += sum(map(tr,outputs))+sum(stops.values())
    close(witness_probability, tr(excitation), label='selected continuation full law')

    def energy_channel(state):
        out = np.zeros_like(state)
        for group in groups:
            ids, basis = archive[group['position_indices_archive_key']], archive[group['basis_archive_key']]
            c = basis.conj().T @ state[np.ix_(ids,ids)] @ basis
            out[np.ix_(ids,ids)] += basis @ c @ basis.conj().T
        return out

    def regional_channel(state):
        return state*(regional[:,None] == regional[None,:])

    def distance(a,b):
        # Inputs in this audit are number-block diagonal. Hermitian sector
        # eigenspectra replace a full384 singular-value decomposition.
        return sum(float(abs(np.linalg.eigvalsh((a-b)[np.ix_(ids,ids)])).sum()) for ids,*_ in model.spectral.values())/2

    actual = ca['final_state']
    er, re = regional_channel(energy_channel(actual)), energy_channel(regional_channel(actual))
    close(er, archive['energy_then_regional_poststate'], label='energy then regional')
    close(re, archive['regional_then_energy_poststate'], label='regional then energy')
    order_distance = distance(er,re)
    close(order_distance, result['mutations']['regional_energy_order_trace_distance'], label='measurement order disturbance')
    global_commutator = float(np.max(abs(h*(numbers[None,:]-numbers[:,None]))))
    regional_commutator = max(float(np.linalg.norm((h*(regional[None,:]-regional[:,None]))[np.ix_(ids,ids)],2)) for ids,*_ in model.spectral.values())
    close(global_commutator, 0, label='global Q/H commute')
    close(regional_commutator, result['mutations']['regional_Q_H_commutator_operator_norm'], label='regional Q/H commutator')
    change = distance(sum(after), excitation)/tr(excitation)
    rate_zero = independent_clock.counts([excitation], rate=0)[0].sum(axis=0)
    intervention = distance(sum(after), rate_zero)/tr(excitation)
    close(change, selected['conditional_trace_distance_after_continuation'], label='conditional poststate change')
    close(intervention, selected['native_rate_on_off_conditional_trace_distance'], label='rate intervention')
    charged_gap = min(g['energy']-ground for g in groups if g['global_charge'] == 1)
    global_gap = min(g['energy']-ground for g in groups if g['energy']-ground > 1e-9)
    close(charged_gap, .3480195548854885, label='charged spectral gap')
    close(global_gap, .2492618078492539, label='whole-source spectral gap')
    require(charged_gap > global_gap, 'charged gap cannot replace global gap')
    print('All 576 regional records and 200448 compact energy states independently checked.', flush=True)

    mutations = {}
    def rejected(name, test):
        try:
            test()
        except ValueError:
            mutations[name] = True
        else:
            raise ValueError('mutation accepted: '+name)
    # Real altered instruments on coherent source inputs, not boolean flags.
    control = np.zeros((8,8), complex)
    control[np.ix_([1,2],[1,2])] = .5
    good = np.einsum('abij,ij->ab', channels[1],control)
    fine = np.diag(np.diag(good))
    coarse_mutation_distance = float(abs(np.linalg.eigvalsh(good-fine)).sum()/2)
    rejected('fine_local_measurement_then_forget', lambda: close(fine,good,label='wrong coarse instrument'))
    doublet = next((i,g) for i,g in enumerate(groups) if g['degeneracy'] == 2)
    dgi, dg = doublet
    dvec = archive[dg['basis_archive_key']]
    c = np.full((2,2),.5,complex)
    whole, fine_energy = dvec @ c @ dvec.conj().T, dvec @ np.diag(np.diag(c)) @ dvec.conj().T
    energy_mutation_distance = float(abs(np.linalg.eigvalsh(whole-fine_energy)).sum()/2)
    rejected('rank_one_energy_then_forget', lambda: close(whole,fine_energy,label='wrong energy instrument'))
    changed_result = copy.deepcopy(result)
    changed_result['first_detector_records'][0]['stops'].pop('oracle')
    rejected('drop_export_flag_even_when_tiny', lambda: exact_join(contract,clock_contract,changed_result))
    rejected('drop_reset_timeout_probability', lambda: close(terminals.sum(),1.,label='drop all failures'))
    sample_abort = next(archive[r['local_kraus_key']] for family in result['export_failure_kraus_families'].values() for r in family if np.any(archive[r['local_kraus_key']]))
    rejected('zero_rectangular_abort_map', lambda: require(np.max(abs(sample_abort-0))/np.max(abs(sample_abort)) < 2e-12, 'abort map erased'))
    tail_probability = float(terminals.reshape(12,4,3,4)[:,:,2,:].sum())
    rejected('drop_complete_count_tail', lambda: close(terminals.sum()+all_stops-tail_probability,1.,label='tail deleted'))
    rejected('success_only_renormalization', lambda: close(terminals.sum()/terminals.sum()+all_stops,1.,label='success-only law'))
    wrong_clock = copy.deepcopy(contract)
    wrong_clock['clock_id'] = 'passive192CAR'
    rejected('wrong_clock_identifier', lambda: exact_join(wrong_clock,clock_contract,result))
    changed_rate = independent_clock.counts([excitation], rate=.4)[0]
    wrong_clock_distance = distance(changed_rate.sum(axis=0),after.sum(axis=0))/tr(excitation)
    rejected('wrong_clock_actual_half_rate', lambda: close(changed_rate,after,label='actual wrong clock'))
    calibration = contract['calibration']
    def calibration_check(scale):
        close(scale*charged_gap, calibration['supplied_anchor_GeV'], 1e-15, 'energy unit mismatch')
    calibration_check(calibration['GeV_per_source_energy'])
    rejected('MeV_GeV_factor1000_mismatch', lambda: calibration_check(1000*calibration['GeV_per_source_energy']))
    require(all(mutations.values()), 'mutation test failure')
    source_paths = set(FROZEN) | {Path(__file__), HERE/'PROOF.md', HERE/'README.md', PART/'particle_states.npz', PART/'particle_chain.py',
        CLOCK/'CLOCK_CONTRACT.json', CLOCK/'common_clock_states.npz', CLOCK/'common_clock.py',
        OP/'operational_states.npz', NEXT/'C10/joint_source_states.npz', NEXT/'C11/CHARGED_ENERGY_UNITS.json'}
    for doc in (contract,clock_contract,operational):
        source_paths |= {ROOT / p for p in doc['source_sha256']}
    return dict(schema='oph.independent_particle_audit.v1', reviewer='worker09_independent_particle_measurement_auditor',
        producer_verification_counted_as_independent=False,
        audit_method='Independent basis-index M6 pulse and literal C/E pointer replay; direct row-major block Liouvillians; all compact grouped-energy matrices; no producer run/verify/channel/clock helpers.',
        source_sha256={str(p.relative_to(ROOT)):digest(p) for p in sorted(source_paths)},
        checks=dict(frozen_source_preparation_reference_clock_join=True, literal_native_pulse_replay=True,
            complete_first_failure_maps_and_reset_law=True, arbitrary_reference_coarse_instrument=True,
            all12_clock_inputs_and144_continuations=True, all576_regional_and200448_energy_states=True,
            whole_degenerate_Luders_coherences=True, complete_failure_zero_and_tail_outcomes=True,
            actual_selected_poststate_same_clock_continuation=True, regional_energy_noncommutation=True,
            charged_and_whole_spectrum_gaps_distinguished=True, no_mass_stability_or_clock_award=True),
        native_replay=native_checks,
        complete_records=dict(clock_inputs=12, second_detector_inputs=144, regional_outcomes=576,
            grouped_energy_outcomes=200448, energy_groups=348, energy_doublets=36,
            exactly_zero_energy_probabilities=int(np.count_nonzero(probabilities == 0)),
            energy_outcomes_below_1e_minus14_probability=int(np.count_nonzero(probabilities < 1e-14)),
            total_probability=float(terminals.sum()+all_stops), all_stop_probability=all_stops,
            complete_second_count_tail_probability=tail_probability),
        errors=dict(first_clock_max_entry=clock_error, continued_clock_max_entry=continuation_error,
            all_energy_coefficients_max_entry=maximum_coefficient_error,
            actual_eigenspace_max_entry=spectral_error, selected_state_max_entry=selected_error,
            minimum_compact_state_eigenvalue=min_coefficient_eigenvalue,
            minimum_continued_state_eigenvalue=minimum_continuation_eigenvalue,
            maximum_discarded_offnumber_entry=maximum_offnumber_entry,
            maximum_actual_retained_doublet_offdiagonal=actual_degenerate_offdiagonal),
        witnesses=dict(charge_one_gap=charged_gap, whole_source_gap=global_gap,
            regional_energy_order_trace_distance=order_distance, regional_H_commutator_operator_norm=regional_commutator,
            selected_conditional_continuation_distance=change, native_rate_on_off_conditional_distance=intervention),
        mutation_rejections=mutations,
        mutation_witnesses=dict(wrong_coarse_instrument_trace_distance=coarse_mutation_distance,
            fine_energy_forget_trace_distance=energy_mutation_distance,
            half_rate_clock_conditional_trace_distance=wrong_clock_distance,
            dropped_reset_stop_probability=all_stops, dropped_second_count_tail_probability=tail_probability),
        interpretation='Finite simulation-generated ordered regional charge, energy and transport records on the declared source. Proper-code transfer, routing, pulse angles, spectral apparatus and clock parameters are supplied services. The one GeV anchor is calibration. No physical species, mass, stability, lifetime, time conversion or continuum transfer follows.',
        official_rows_changed=[])


if __name__ == '__main__':
    report = run()
    (HERE/'AUDIT_RESULT.json').write_text(json.dumps(report,indent=2)+'\n')
    particle = read(PART/'PARTICLE_CONTRACT.json')
    audit_contract = dict(schema='oph.independent_particle_audit_contract.v1',
        audit_id='worker09_C10_particle_measurement_independent_v1',
        audit_author='worker09_independent_particle_measurement_auditor',
        reviewed_producer='worker04_particle_chain',
        admission='PASS_FINITE_DECLARED_SOURCE_STATE_CLOCK_PARTICLE_INSTRUMENT',
        independence=dict(separate_author_from_particle_producer=True,
            producer_verify_used_as_independent_evidence=False,
            producer_run_or_measurement_helpers_imported=False,
            CommonClock_imported=False,
            shared_declared_primitives=['JointNativeSource', 'M6 source pulse response', 'literal native C/E event', 'native reset policy']),
        **{key:particle[key] for key in ('source_id','source_dimension','source_H_array_sha256_le_float64',
            'preparation_id','preparation_array_sha256_le_complex128','reference_id','clock_id','particle_instrument_id')},
        particle_contract=str((PART/'PARTICLE_CONTRACT.json').relative_to(ROOT)),
        particle_contract_sha256=digest(PART/'PARTICLE_CONTRACT.json'),
        particle_result_sha256=digest(PART/'RESULT.json'),
        particle_state_archive_sha256=digest(PART/'particle_states.npz'),
        operational_readout_receipt_sha256=digest(OP/'OPERATIONAL_QUANTUM_REVIEW.json'),
        result='AUDIT_RESULT.json', result_sha256=digest(HERE/'AUDIT_RESULT.json'),
        source_sha256=report['source_sha256'],
        checks=report['checks'], mutation_rejections=report['mutation_rejections'],
        admissible_domain=dict(source='actual full384 C10 source and complete fair-reservoir preparation',
            clock='same supplied T=.37 and lambda=.8 complete marked Poisson law, all zero-current attempts counted',
            detector='replayed 114-pulse regional population count and two literal C-pointer readouts, exact reset and all abort flags',
            energy='declared grouped joint global-Q/H Luders apparatus retaining whole degeneracies',
            records='all 12 initial count states, 144 continuations, 576 regional records, 200448 grouped energy states',
            transport='finite ordered regional count change; no tagged worldline or individual particle trajectory'),
        claim_boundary=particle['claim_boundary'], official_rows_changed=[])
    (HERE/'AUDIT_CONTRACT.json').write_text(json.dumps(audit_contract,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('checks','complete_records','errors','witnesses','mutation_rejections')},indent=2))
