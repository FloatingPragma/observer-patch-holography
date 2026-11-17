"""Learn charge weights through a finite native current-record census.

On each fixed tree edge, stop at the first coin different from the first
coin, or retain a deadline abort.  The first current and the current at that
flip determine the pre-edge contrast.  Recorded transfers then convert all
tree contrasts into contrasts of the original state.  Total occupation two
fixes its additive constant.  This is a destructive native measurement,
not access to hidden occupations or a supplied projective work-bit read.

Each completed edge has probability alpha=1-2**(1-h), independent of state
and reference.  Its two compressed coin orders have weight alpha/2.  A
timeout has two constant-coin words of weight 2**(-h), with their actual
disturbed outputs.  Complete successful record effects therefore identify
one original basis state; the coarse charge effects are alpha**11 Pi_q.
The source-derived preparation controls coherences and reference memory.

Bernoulli martingale Bernstein bounds use the unknown target probability
only inside their proof.  Algebraic inversion yields an empirical radius;
no known charge weights, independent copies or supplied spectral floor
enter estimation.  Seeded kernel sampling is reproducible software evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import reference_tomography as old, local_preparation as local
from ..native_interface import construction as interface
from ..observer import native, quantum_mixing

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('native_census_results.json')
SCHEMA = 'oph.native_modular_reference.native_census.v1'
TRAINING_TRIALS = 1 << 20
DEADLINE = 32
GRID = 1 << 53
SAMPLER_GATE = F(1, 10**10)
DIMENSIONS = (45, 20, 1)
LABELS = ('success_q0', 'success_q1', 'success_q2', 'herald_failure',
          'warmup_timeout', 'completion_timeout', 'census_timeout', 'numerical_abort')
# A fixed source-edge tree, selected before any record is observed.
TREE = ((0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (2, 10),
        (2, 3), (2, 4), (2, 6), (1, 8), (1, 9))
SOURCE_PATHS = tuple(sorted(set((*local.SOURCE_PATHS, *interface.SOURCE_PATHS,
    'model/native_modular_reference/native_census.py',
    'model/native_modular_reference/reference_tomography.py',
    'model/native_equilibration/common_variation.py'))))


def _deadline(value):
    old._integer(value, 2, 64, 'census deadline')
    return value


def success_probability(deadline=DEADLINE):
    return (1-F(1, 2**(_deadline(deadline)-1)))**len(TREE)


@lru_cache(maxsize=1)
def _basis():
    return native.states(binary=True)


def decode(records):
    """Decode only eleven compressed completed-edge coin/current records.

    Each record contains (first_coin, first_current, flip_current).  Repeated
    copies of the first coin have zero current and need not be stored here.
    Neither native states nor source reference weights are decoder inputs.
    """
    if type(records) is not tuple or len(records) != len(TREE):
        raise ValueError('eleven completed edge records required')
    displacement = [0]*12
    differences = []
    for (a, b), record in zip(TREE, records, strict=True):
        if (type(record) is not tuple or len(record) != 3
                or any(type(v) is not int for v in record)
                or record[0] not in (0, 1) or any(j not in (-1, 0, 1) for j in record[1:])):
            raise ValueError('original native coin/current triples required')
        coin, first, flip = record
        direction = 2*coin-1
        if (first not in (0, -direction) or flip not in (0, direction)
                or first != 0 and flip != direction):
            raise ValueError('impossible binary native edge record')
        contrast = first if first else (direction if flip else 0)
        differences.append((a, b, contrast-displacement[a]+displacement[b]))
        displacement[a] -= first+flip
        displacement[b] += first+flip
    potential = {0: 0}
    while len(potential) < 12:
        before = len(potential)
        for a, b, contrast in differences:
            if a in potential:
                wanted = potential[a]-contrast
                if b in potential and potential[b] != wanted:
                    raise ValueError('inconsistent tree contrasts')
                potential[b] = wanted
            elif b in potential:
                potential[a] = potential[b]+contrast
        if len(potential) == before:
            raise ArithmeticError('census tree is disconnected')
    base = F(2-sum(potential.values()), 12)
    initial = tuple(base+potential[i] for i in range(12))
    final = tuple(initial[i]+displacement[i] for i in range(12))
    if any(x not in (0, 1) for x in (*initial, *final)) or sum(final) != 2:
        raise ValueError('records do not reconstruct a binary total-two state')
    return tuple(map(int, initial)), tuple(map(int, final))


def execute_record(state, first_coins):
    """Source-side compressed successful execution; decode owns observation.

    This simulator receives a state; the decoder above receives only records.
    Actual deadline waiting is integrated exactly by the finite instrument.
    """
    if (type(state) is not tuple or any(type(x) is not int for x in state)
            or state not in _basis()):
        raise ValueError('original binary native state required')
    if (type(first_coins) is not tuple or len(first_coins) != len(TREE)
            or any(type(c) is not int or c not in (0, 1) for c in first_coins)):
        raise ValueError('eleven original binary coins required')
    records = []
    for seam, coin in zip(TREE, first_coins, strict=True):
        state, (_, first) = native.native(state, seam, coin)
        state, (_, flip) = native.native(state, seam, 1-coin)
        records.append((coin, first, flip))
    return tuple(records), state


@lru_cache(maxsize=1)
def _transitions():
    states = native.states(binary=True)
    lookup = {x: i for i, x in enumerate(states)}
    result = []
    for seam in TREE:
        first, completed = [], []
        for coin in (0, 1):
            single = [native.native(x, seam, coin)[0] for x in states]
            first.append(tuple(lookup[x] for x in single))
            completed.append(tuple(lookup[native.native(x, seam, 1-coin)[0]] for x in single))
        result.append((tuple(first), tuple(completed)))
    return tuple(result)


def _push(counts, maps):
    output = np.zeros_like(counts)
    for mapping in maps:
        np.add.at(output, np.asarray(mapping), counts)
    return output


@lru_cache(maxsize=1)
def _count_data():
    """Integer population numerators; no fine-history enumeration."""
    states = native.states(binary=True)
    charge = np.array([sum(x[p] for p in local.CORE) for x in states])
    successes = [np.diag((charge == q).astype(np.int64)) for q in range(3)]
    prefix = np.eye(66, dtype=np.int64)
    aborts = []
    for first, completed in _transitions():
        aborts.append(_push(prefix, first).tobytes())
        prefix = _push(prefix, completed)
        successes = [_push(m, completed) for m in successes]
    if not np.all(sum(successes).sum(axis=0) == 2**len(TREE)):
        raise ArithmeticError('census population instrument lost a column')
    return tuple(m.tobytes() for m in successes), tuple(aborts)


def population_maps(deadline=DEADLINE):
    """Actual finite census maps, including every timeout's disturbed state."""
    alpha = 1-F(1, 2**(_deadline(deadline)-1))
    success, aborts = _count_data()
    unpack = lambda b: np.frombuffer(b, dtype=np.int64).reshape(66, 66)
    maps = {f'success_q{q}': float(alpha**11/2**11)*unpack(data)
            for q, data in enumerate(success)}
    maps['census_timeout'] = sum((float(alpha**i*(1-alpha)/2**(i+1))*unpack(data)
                                 for i, data in enumerate(aborts)), np.zeros((66, 66)))
    return maps


@lru_cache(maxsize=88)
def _edge_branches(law, seam, first_coin, *, complete):
    """Exact compressed partial injections of the actual C/E instruments."""
    branches = interface.fixed_maps(law, seam)
    initial = [b for b in branches if b.coin == first_coin]
    if not complete:
        return tuple((b.current, None, b.mapping) for b in initial)
    later = [b for b in branches if b.coin == 1-first_coin]
    result = []
    for before in initial:
        for after in later:
            lookup = dict(after.mapping)
            mapping = tuple((x, lookup[y]) for x, y in before.mapping if y in lookup)
            if mapping:
                result.append((before.current, after.current, mapping))
    return tuple(result)


def _channel_signature(mappings):
    """Exact integer superoperator of conditional-coin partial injections."""
    signature = Counter()
    for mapping in mappings:
        for left, out_left in mapping:
            for right, out_right in mapping:
                signature[left, right, out_left, out_right] += 1
    return signature


def _source_edge_audit(law, seam):
    states = _basis()
    lookup = {x: i for i, x in enumerate(states)}
    expected = defaultdict(list)
    currents = {}
    for coin in (0, 1):
        for i, state in enumerate(states):
            output, (_, current) = native.native(state, seam, coin)
            label = current if law == 'C' else state[seam[0]]-state[seam[1]]
            expected[coin, label].append((i, lookup[output]))
            currents[coin, label] = current
    actual = interface.fixed_maps(law, seam)
    if len(actual) != len(expected) or len({(b.coin, b.label) for b in actual}) != len(actual):
        raise ArithmeticError('native census fine-record partition changed')
    for branch in actual:
        key = branch.coin, branch.label
        if (key not in expected or branch.current != currents[key]
                or branch.mapping != tuple(expected[key])):
            raise ArithmeticError('native census output or fine-record partition changed')
        matrix = np.zeros((66, 66), complex)
        for i, j in expected[key]:
            matrix[j, i] = 1/math.sqrt(2)
        if type(branch.matrix) is not np.ndarray or not np.array_equal(branch.matrix, matrix):
            raise ArithmeticError('native census Kraus phases or amplitudes changed')
    for coin in (0, 1):
        first = tuple(b.mapping for b in actual if b.coin == coin)
        repeated = []
        for before in first:
            for after in first:
                later = dict(after)
                mapping = tuple((i, later[j]) for i, j in before if j in later)
                if mapping:
                    repeated.append(mapping)
        if _channel_signature(first) != _channel_signature(repeated):
            raise ArithmeticError('same-coin quantum compression changes a matrix unit')
    return True


def execute(density, *, law='C', deadline=DEADLINE):
    """Complete finite census on quantum inputs and untouched references.

    Successful fine effects identify the original basis state, so their
    reference blocks propagate through the certified population matrices.
    Timeout branches retain quantum coherences left by their actual prefix.
    """
    from . import local_history as history
    density = history._density(density)
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('one original native C/E law required')
    alpha = 1-F(1, 2**(_deadline(deadline)-1))
    ref = len(density)//66
    def push(rho, seam, complete):
        result = np.zeros_like(rho, dtype=complex)
        for coin in (0, 1):
            for _, _, mapping in _edge_branches(law, seam, coin, complete=complete):
                before, after = zip(*mapping, strict=True)
                columns = np.array([r+ref*i for i in before for r in range(ref)])
                rows = np.array([r+ref*i for i in after for r in range(ref)])
                result[np.ix_(rows, rows)] += rho[np.ix_(columns, columns)]/2
        return result
    continuing = density
    timeout = np.zeros_like(density, dtype=complex)
    for seam in TREE:
        timeout += float(1-alpha)*push(continuing, seam, False)
        continuing = float(alpha)*push(continuing, seam, True)
    maps = population_maps(deadline)
    tensor = density.reshape(66, ref, 66, ref)
    result = {'census_timeout': timeout}
    for q in range(3):
        output = np.zeros((66, ref, 66, ref), complex)
        for y in range(66):
            output[y, :, y, :] = sum((maps[f'success_q{q}'][y, x]*tensor[x, :, x, :]
                                      for x in range(66)), np.zeros((ref, ref), complex))
        result[f'success_q{q}'] = output.reshape(density.shape)
    if any(not np.isfinite(m).all() for m in result.values()):
        raise ArithmeticError('nonfinite complete census output')
    # The independently propagated unread successful channel must agree with
    # the fine-record dephasing theorem, including reference coherences.
    if np.linalg.norm(continuing-sum(result[f'success_q{q}'] for q in range(3)), 'fro') > 1e-12:
        raise ArithmeticError('census singleton-effect channel disagrees with native composition')
    return result


def _compute_identification():
    """Source support proof on every edge and all 4,356 matrix units.

    Every edge branch is a partial injection. Its record fixes the contrast,
    and repetitions of the first coin neither change state nor split a
    branch. The tree-decoder identity then proves successful fine records
    have singleton input support: all off-diagonal effects vanish, including
    arbitrary untouched reference blocks. No floating channel norm is used.
    """
    states = native.states(binary=True)
    if len(TREE) != 11 or len(set(TREE)) != 11 or any(e not in quantum_mixing.seams() for e in TREE):
        raise ArithmeticError('eleven distinct native source seams required')
    connected = {0}
    for _ in range(12):
        for a, b in TREE:
            if a in connected or b in connected:
                connected.update((a, b))
    if len(connected) != 12:
        raise ArithmeticError('source tree does not span the carrier')
    laws = {}
    for law in ('C', 'E'):
        branch_checks = matrix_checks = 0
        for seam in TREE:
            _source_edge_audit(law, seam)
            for coin in (0, 1):
                for complete in (False, True):
                    branches = _edge_branches(law, seam, coin, complete=complete)
                    seen = [0]*66
                    for first, flip, mapping in branches:
                        inputs, outputs = zip(*mapping, strict=True)
                        if len(set(inputs)) != len(inputs) or len(set(outputs)) != len(outputs):
                            raise ArithmeticError('native census branch is not a partial injection')
                        for x, y in mapping:
                            seen[x] += 1
                            once, (_, current) = native.native(states[x], seam, coin)
                            repeated, repeated_record = native.native(once, seam, coin)
                            if repeated != once or repeated_record[1] != 0 or current != first:
                                raise ArithmeticError('same-coin compression changes a native path')
                            if complete:
                                terminal, (_, second) = native.native(once, seam, 1-coin)
                                if terminal != states[y] or second != flip:
                                    raise ArithmeticError('compressed opposite-coin output changed')
                                contrast = first if first else ((2*coin-1) if flip else 0)
                                if contrast != states[x][seam[0]]-states[x][seam[1]]:
                                    raise ArithmeticError('current pair does not determine contrast')
                            elif once != states[y]:
                                raise ArithmeticError('compressed timeout output changed')
                        branch_checks += 1
                    if seen != [1]*66:
                        raise ArithmeticError('compressed branch partition is incomplete')
                    # Every matrix-unit effect is diagonal: a partial injection
                    # has distinct outputs, while its diagonal partition is I.
                    lookup_maps = tuple(dict(b[2]) for b in branches)
                    for i in range(66):
                        for j in range(66):
                            effect = sum(int(i in m and j in m and m[i] == m[j]) for m in lookup_maps)
                            if effect != int(i == j):
                                raise ArithmeticError('complete edge effect differs from identity')
                            matrix_checks += 1
        laws[law] = {'partial_injection_branches': branch_checks,
                     'matrix_unit_effect_checks': matrix_checks,
                     'same_coin_quantum_channel_idempotence': True,
                     'original_fine_kraus_matrices_and_terminal_outputs_checked': True,
                     'arbitrary_untouched_references': True}
    # Exhaust all original states and all coin orders in the decoder audit.
    # This checks the observation algorithm; instrument propagation uses the
    # integer matrix recurrence above, not this enumeration.
    for x in states:
        for word in range(1 << len(TREE)):
            records, final = execute_record(x, tuple((word >> i) & 1 for i in range(11)))
            if decode(records) != (x, final):
                raise ArithmeticError('record-only decoder differs from source execution')
    return {'tree': TREE, 'deadline': DEADLINE, 'laws': laws,
            'successful_original_states_checked': 66*(1 << 11),
            'charge_effect_factor': success_probability(),
            'charge_multiplicities': DIMENSIONS,
            'coarse_effect': 'alpha^11 times original-core-charge projector',
            'timeout_retains_disturbed_work': True,
            'controller_reads_hidden_occupations': False,
            'full_coherent_state_tomography': False}


@lru_cache(maxsize=1)
def _identification_bytes():
    """Immutable cached exhaustive decoder evidence; never return its object."""
    return json.dumps(old._json(_compute_identification()), sort_keys=True).encode('utf-8')


def exact_identification():
    # Repeat the live source branch audits even when exhaustive decoder proof
    # bytes are cached. Runtime changes to phases, amplitudes, fine records or
    # quantum compression cannot hide behind a previously checked receipt.
    for law in ('C', 'E'):
        for seam in TREE:
            _source_edge_audit(law, seam)
    row = json.loads(_identification_bytes())
    row['tree'] = tuple(tuple(edge) for edge in row['tree'])
    row['charge_multiplicities'] = tuple(row['charge_multiplicities'])
    row['charge_effect_factor'] = F(row['charge_effect_factor'])
    return row


def resources(trials=TRAINING_TRIALS):
    old._integer(trials, 1, 1 << 40, 'trial count')
    source_error, q, maximum = old._record_contract()
    # This is the same source-derived expected-work proof as local.budget;
    # obtaining expected counters does not evaluate a reference density.
    warm_fraction = F(len(local.incident.incident_edges(local.CORE)), 30)
    _, complete_fraction, gamma = local._schedule()
    warm_blocks = complete_blocks = 0
    while F(33)/warm_fraction/2**warm_blocks > F(1, 2**24):
        warm_blocks += 1
    while F(33)/complete_fraction/2**complete_blocks > F(1, 2**24):
        complete_blocks += 1
    warm = F(16*warm_blocks+1)/warm_fraction
    complete = (q+source_error)*math.ceil(1/gamma)*complete_blocks/complete_fraction
    expected = F(19, 3)*warm+F(34, 5)*complete
    # One finite logical layout, not measured Python/interpreter memory.
    # On valid source paths, displacement=current-initial and reconstructed
    # potential=initial_i-initial_root each lie in {-1,0,1}.  Their sum lies
    # in [-12,12]. Duplicate working registers make this a conservative
    # reusable layout; the fixed tree/program and training counters are extra.
    registers = {'eleven_compressed_coin_current_triples': 11*5,
                 'twelve_signed_displacements': 12*2,
                 'twelve_reconstructed_potentials': 12*2,
                 'discovered_port_mask': 12,
                 'within_edge_attempt_counter': DEADLINE.bit_length(),
                 'tree_edge_counter': len(TREE).bit_length(),
                 'working_coin': 1, 'working_current': 2,
                 'terminal_flag': 3, 'decoded_core_charge': 2,
                 'reconstructed_root_occupation': 1, 'signed_potential_sum': 5}
    return {'preparation_native_attempts_upper': maximum,
            'preparation_expected_native_attempts_upper': expected,
            'unread_core_attempts_upper': 1,
            'census_native_attempts_upper': len(TREE)*DEADLINE,
            'census_expected_native_attempts_upper': F(3*len(TREE)),
            'source_native_attempts_per_trial_upper': maximum+1+len(TREE)*DEADLINE,
            'source_expected_native_attempts_per_trial_upper': expected+1+3*len(TREE),
            'census_controller_register_bits': registers,
            'census_controller_register_bits_subtotal': sum(registers.values()),
            'training_terminal_counter_bits': len(LABELS)*trials.bit_length(),
            'training_trial_counter_bits': trials.bit_length(),
            'memory_counts_are_logical_layout_not_measured_python_memory': True,
            'census_subtotal_excludes': ('tree_and_program_tables', 'interpreter_workspace',
                'private_fine_tapes', 'training_counts', 'source_preparation_controller',
                'quantum_registers', 'host_sampler_arrays_and_receipt_state_counts'),
            'host_receipt_terminal_state_counts_are_observer_inputs': False}


def infer_counts(counts, *, sampling_bias=F(), confidence_bits=10):
    """Target-blind rational confidence box from all eight terminal counts.

    Let p be an ideal complete-trial charge probability, eps the uniform
    conditional bias, f its frequency, c=2k/N, A=eps+2k/(3N). Bernstein
    implies x=|f-p| <= A+sqrt(c(p+eps)). Since p<=f+x, the positive root is
    A+c/2+sqrt(c(f+eps+A)+c^2/4). Each radius is rounded upward rationally.
    Three two-sided bounds fail with probability at most 6*2^(-k).
    """
    if type(counts) is not dict or set(counts) != set(LABELS):
        raise ValueError('all eight complete terminal count categories required')
    if any(type(v) is not int or v < 0 for v in counts.values()):
        raise ValueError('terminal counts must be nonnegative original integers')
    n = sum(counts.values())
    old._integer(n, 1, 1 << 40, 'trial count')
    old._integer(confidence_bits, 4, 64, 'confidence bits')
    sampling_bias = old._rational(sampling_bias)
    if sampling_bias > 1:
        raise ValueError('sampling bias must be at most one')
    epsilon = old._record_contract()[0]+sampling_bias
    frequencies = tuple(F(counts[f'success_q{q}'], n) for q in range(3))
    c, a = F(2*confidence_bits, n), epsilon+F(2*confidence_bits, 3*n)
    radii = tuple(old._up(a+c/2+old._sqrt_upper(c*(f+epsilon+a)+c*c/4)) for f in frequencies)
    centers = tuple(f/g for f, g in zip(frequencies, DIMENSIONS, strict=True))
    intervals = tuple(((f-r)/g, (f+r)/g) for f, r, g in zip(frequencies, radii, DIMENSIONS, strict=True))
    budget = resources(n)
    result = {'status': 'unresolved', 'trials': n,
              'central_subnormalized_weights': centers,
              'subnormalized_weight_intervals': intervals,
              'charge_probability_radii': radii, 'observed_success_mass': sum(frequencies),
              'systematic_halftrace_upper': epsilon, 'confidence_bits': confidence_bits,
              'confidence_failure_upper': F(6, 2**confidence_bits),
              'census_success_factor': success_probability(),
              **budget,
              'source_native_attempts_all_trials_upper': n*budget['source_native_attempts_per_trial_upper'],
              'source_expected_native_attempts_all_trials_upper': n*budget['source_expected_native_attempts_per_trial_upper'],
              'independent_input_copies_required': False,
              'source_or_hardware_execution_authenticated_by_counts': False}
    if any(lo <= 0 for lo, _ in intervals):
        result['reason'] = 'complete records do not certify three positive weights'
        return result
    total = sum(frequencies)
    normalized = tuple(w/total for w in centers)
    result.update(status='faithful', central_weights=normalized,
                  moments=(20*normalized[1]+2*normalized[2], normalized[2]),
                  normalization_scalar_ignored_only_in_modular_channel=True)
    return result


def _raw_kernel(law):
    prep = local.population_maps(law=law)
    core = old.population_counts()[0]/2
    census = population_maps()
    entered = core@prep['success']
    arrays = [census[label]@entered for label in LABELS[:3]]
    arrays += [prep[label] for label in LABELS[3:6]]
    arrays += [census['census_timeout']@entered]
    raw = np.vstack(arrays)
    if raw.shape != (462, 66) or not np.isfinite(raw).all() or np.any(raw < 0):
        raise ArithmeticError('invalid complete finite census population kernel')
    residual = float(np.max(np.abs(raw.sum(axis=0)-1)))
    if residual > 5e-10:
        raise ArithmeticError('finite census kernel lost mass')
    return raw, residual


def quantized_kernel(law='C'):
    raw, residual = _raw_kernel(law)
    cdf = np.cumsum(raw, axis=0)
    result = np.zeros((528, 66), dtype=np.uint64)
    for x in range(66):
        previous = 0
        for y in range(462):
            scaled = F(float(cdf[y, x]))*GRID
            threshold = min(GRID, -(-scaled.numerator//scaled.denominator))
            if threshold < previous:
                raise ArithmeticError('original binary census CDF is not monotone')
            result[y, x] = threshold-previous
            previous = threshold
        result[462+x, x] = GRID-previous
    _validate_kernel(result)
    return result, {'binary_population_mass_residual_diagnostic': residual,
                    'grid': GRID, 'numerical_abort_is_separate_terminal': True}


def _ideal_distribution():
    """Independent rational source target, used only by sampler validation."""
    q = old._record_contract()[1]
    sigma = local.reference()
    raw = local.incident.reference(local.CORE)[0]
    successes, aborts = _count_data()
    unpack = lambda b: np.frombuffer(b, dtype=np.int64).reshape(66, 66)
    alpha = 1-F(1, 2**(DEADLINE-1))
    apply = lambda m: [sum((F(int(c))*w for c, w in zip(row, sigma, strict=True)), F()) for row in m]
    result = []
    for data in successes:
        result += [q*alpha**11*x/2**11 for x in apply(unpack(data))]
    result += [F(1, 66)-q*w for w in raw]
    result += [F()]*(2*66)
    aborted = [F()]*66
    for i, data in enumerate(aborts):
        for y, x in enumerate(apply(unpack(data))):
            aborted[y] += q*alpha**i*(1-alpha)*x/2**(i+1)
    result += aborted+[F()]*66
    if len(result) != 528 or sum(result) != 1 or min(result) < 0:
        raise ArithmeticError('independent census target is incomplete')
    return tuple(result)


def _validate_kernel(kernel):
    if (type(kernel) is not np.ndarray or kernel.shape != (528, 66) or kernel.dtype != np.uint64
            or np.any(kernel > GRID) or any(sum(map(int, kernel[:, x])) != GRID for x in range(66))):
        raise ValueError('complete original uint64 census terminal/state kernel required')


def audit_kernel(kernel):
    _validate_kernel(kernel)
    ideal = _ideal_distribution()
    target = tuple(sum(ideal[66*i:66*(i+1)], F()) for i in range(3))
    worst, biases = F(), [F()]*3
    for x in range(66):
        column = tuple(F(int(v), GRID) for v in kernel[:, x])
        worst = max(worst, sum((abs(a-b) for a, b in zip(column, ideal, strict=True)), F())/2)
        for q in range(3):
            biases[q] = max(biases[q], abs(sum(column[66*q:66*(q+1)], F())-target[q]))
    if any(b > worst for b in biases):
        raise ArithmeticError('charge bias exceeds complete sampler variation')
    return {'exact_dyadic_columns_checked': 66, 'complete_total_variation_upper': worst,
            'charge_mean_bias_upper': tuple(biases),
            'bernoulli_second_moment_equals_mean': True,
            'entire_multitrial_path_distance_claimed': False}


def _portable_sampler_certificate():
    return {'exact_dyadic_columns_checked': 66,
            'complete_total_variation_upper': F(SAMPLER_GATE.numerator, SAMPLER_GATE.denominator),
            'charge_mean_bias_upper': tuple(F(SAMPLER_GATE.numerator, SAMPLER_GATE.denominator)
                                            for _ in range(3)),
            'bernoulli_second_moment_equals_mean': True,
            'entire_multitrial_path_distance_claimed': False}


def sample_kernel(kernel, trials, *, seed=20261025, initial=0, chunk_size=1 << 16):
    _validate_kernel(kernel)
    old._integer(trials, 1, 1 << 30, 'trial count')
    old._integer(seed, 0, (1 << 63)-1, 'seed')
    old._integer(initial, 0, 65, 'initial state')
    old._integer(chunk_size, 1, 1 << 20, 'chunk size')
    common = kernel.min(axis=1)
    common_mass = int(common.sum())
    common_cdf = np.cumsum(common)
    residual_cdf = np.cumsum(kernel-common[:, None], axis=0)
    rng = np.random.default_rng(seed)
    counts = np.zeros(528, dtype=np.int64)
    digest, prefix, state, exceptional = hashlib.sha256(), [], initial, 0
    for start in range(0, trials, chunk_size):
        draws = rng.integers(0, GRID, size=min(chunk_size, trials-start), dtype=np.uint64)
        selected = np.searchsorted(common_cdf, draws, side='right')
        for position in np.flatnonzero(draws >= common_mass):
            prior = state if position == 0 else int(selected[position-1]) % 66
            selected[position] = np.searchsorted(residual_cdf[:, prior], int(draws[position])-common_mass, side='right')
            exceptional += 1
        if np.any(selected >= 528):
            raise ArithmeticError('census sampler omitted terminal mass')
        state = int(selected[-1]) % 66
        counts += np.bincount(selected, minlength=528)
        digest.update(selected.astype('<u2').tobytes())
        if len(prefix) < 32:
            prefix.extend((LABELS[int(i)//66], int(i)%66) for i in selected[:32-len(prefix)])
    return {'counts': {label: int(counts[66*i:66*(i+1)].sum()) for i, label in enumerate(LABELS)},
            'terminal_state_counts': counts.reshape(8, 66).tolist(),
            'record_and_terminal_state_sha256': digest.hexdigest(), 'prefix': prefix,
            'trials': trials, 'seed': seed, 'initial_state': initial, 'final_state': state,
            'common_subkernel_probability': F(common_mass, GRID), 'residual_draws': exceptional,
            'same_carrier_reused': True, 'microscopic_attempt_tapes_executed': False,
            'seeded_PRNG_demonstration_proves_independent_uniform_randomness': False}


def execute_trials(trials=TRAINING_TRIALS, *, seed=20261025, law='C'):
    kernel, numerical = quantized_kernel(law)
    audit = audit_kernel(kernel)
    sample = sample_kernel(kernel, trials, seed=seed)
    certificate = _portable_sampler_certificate()
    if (audit['complete_total_variation_upper'] > SAMPLER_GATE
            or any(b > SAMPLER_GATE for b in audit['charge_mean_bias_upper'])
            or sample['common_subkernel_probability'] < 1-SAMPLER_GATE):
        raise ArithmeticError('actual census sampler fails portable source-comparison gates')
    del sample['common_subkernel_probability']
    sample['common_subkernel_probability_lower'] = 1-SAMPLER_GATE
    return {'law': law, 'numerical_kernel': numerical, 'sampler_certificate': certificate,
            'execution': sample, 'inference': infer_counts(sample['counts'], sampling_bias=SAMPLER_GATE)}


def run(trials=TRAINING_TRIALS, *, seed=20261025):
    row = execute_trials(trials, seed=seed)
    inferred = row['inference']
    if inferred['status'] == 'faithful':
        target = old._record_contract()[1]*success_probability()
        charges = [sum(x[p] for p in local.CORE) for x in native.states(binary=True)]
        weights = tuple(next(w for w, c in zip(local.reference(), charges, strict=True) if c == q) for q in range(3))
        row['validation_only_known_target_in_confidence_box'] = all(
            lo <= target*w <= hi for w, (lo, hi) in zip(weights, inferred['subnormalized_weight_intervals'], strict=True))
    return {'schema': SCHEMA,
            'identification': exact_identification(), 'continuing_execution': row,
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
            'scope': scope()}


def scope():
    return {'fixed_known_reference_weights_used_by_estimator': False,
            'full_coherent_state_tomography': False, 'independent_copies_required': False,
            'native_current_census_is_destructive': True,
            'finite_physical_service_accuracy_established': False,
            'confidence_theorem_domain': 'actual native records or ideal uniform dyadic-kernel draws',
            'fixed_seed_PRNG_run': 'reproducible illustration, not an independence certificate'}


def compare_receipt(supplied, actual=None, *, trials=TRAINING_TRIALS, seed=20261025):
    """Replay the declared experiment; a caller's second receipt is no oracle.

    Custom runs require explicit trial and seed arguments, independently of
    the retained receipt. Never infer the experiment from untrusted metadata.
    """
    old._integer(trials, 1, 1 << 30, 'sample trial count')
    old._integer(seed, 0, (1 << 63)-1, 'sample seed')
    owned = old._json(run(trials=trials, seed=seed))
    return (_compare_receipt(supplied, owned)
            and (actual is None or _compare_receipt(actual, owned)))


def _compare_receipt(supplied, actual):
    """Compare with an already source-owned replay (used by the CLI)."""
    def canonical(value):
        return json.dumps(value, sort_keys=True, allow_nan=False, separators=(',', ':'))
    identification = old._json(exact_identification())
    custody = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}
    cleaned = []
    for original in (supplied, actual):
        value = deepcopy(original)
        try:
            if (type(value) is not dict or set(value) != {'schema', 'identification', 'continuing_execution', 'source_sha256', 'scope'}
                    or value['schema'] != SCHEMA
                    or canonical(value['identification']) != canonical(identification)
                    or canonical(value['scope']) != canonical(scope())
                    or canonical(value['source_sha256']) != canonical(custody)):
                return False
            row = value['continuing_execution']
            if type(row['law']) is not str or row['law'] not in ('C', 'E'):
                return False
            numerical = row['numerical_kernel']
            if (set(numerical) != {'binary_population_mass_residual_diagnostic', 'grid', 'numerical_abort_is_separate_terminal'}
                    or type(numerical['grid']) is not int or numerical['grid'] != GRID
                    or numerical['numerical_abort_is_separate_terminal'] is not True):
                return False
            diagnostic = row['numerical_kernel']['binary_population_mass_residual_diagnostic']
            if type(diagnostic) is not float or not math.isfinite(diagnostic) or not 0 <= diagnostic <= 5e-10:
                return False
            row['numerical_kernel']['binary_population_mass_residual_diagnostic'] = 'checked'
            if row['sampler_certificate'] != old._json(_portable_sampler_certificate()):
                return False
            if row['execution']['common_subkernel_probability_lower'] != str(1-SAMPLER_GATE):
                return False
            sample = row['execution']
            expected_sample_keys = {'counts', 'terminal_state_counts', 'record_and_terminal_state_sha256',
                'prefix', 'trials', 'seed', 'initial_state', 'final_state', 'residual_draws',
                'same_carrier_reused', 'microscopic_attempt_tapes_executed',
                'seeded_PRNG_demonstration_proves_independent_uniform_randomness',
                'common_subkernel_probability_lower'}
            if (type(sample) is not dict or set(sample) != expected_sample_keys
                    or sample['same_carrier_reused'] is not True
                    or sample['microscopic_attempt_tapes_executed'] is not False
                    or sample['seeded_PRNG_demonstration_proves_independent_uniform_randomness'] is not False):
                return False
            counts, states = sample['counts'], sample['terminal_state_counts']
            if (type(states) is not list or len(states) != len(LABELS)
                    or any(type(block) is not list or len(block) != 66
                           or any(type(n) is not int or n < 0 for n in block) for block in states)
                    or any(sum(block) != counts[label] for label, block in zip(LABELS, states, strict=True))
                    or type(sample['trials']) is not int or sum(counts.values()) != sample['trials']
                    or any(type(sample[key]) is not int or not 0 <= sample[key] < 66
                           for key in ('initial_state', 'final_state'))):
                return False
            old._integer(sample['trials'], 1, 1 << 30, 'sample trial count')
            old._integer(sample['seed'], 0, (1 << 63)-1, 'sample seed')
            old._integer(sample['residual_draws'], 0, sample['trials'], 'residual draw count')
            digest, prefix = sample['record_and_terminal_state_sha256'], sample['prefix']
            if (type(digest) is not str or len(digest) != 64 or any(x not in '0123456789abcdef' for x in digest)
                    or type(prefix) is not list or len(prefix) != min(32, sample['trials'])
                    or any(type(record) is not list or len(record) != 2
                           or type(record[0]) is not str or record[0] not in LABELS
                           or type(record[1]) is not int or not 0 <= record[1] < 66 for record in prefix)):
                return False
            expected = old._json(infer_counts(counts, sampling_bias=SAMPLER_GATE))
            if canonical(row['inference']) != canonical(expected):
                return False
            expected_row_keys = {'law', 'numerical_kernel', 'sampler_certificate', 'execution', 'inference'}
            if expected['status'] == 'faithful':
                expected_row_keys.add('validation_only_known_target_in_confidence_box')
                scale = old._record_contract()[1]*success_probability()
                charges = tuple(sum(x[p] for p in local.CORE) for x in _basis())
                weights = tuple(next(w for w, c in zip(local.reference(), charges, strict=True) if c == q) for q in range(3))
                valid = all(F(lo) <= scale*w <= F(hi) for w, (lo, hi) in zip(
                    weights, expected['subnormalized_weight_intervals'], strict=True))
                if row['validation_only_known_target_in_confidence_box'] is not valid:
                    return False
            if set(row) != expected_row_keys:
                return False
        except (KeyError, TypeError, ValueError):
            return False
        cleaned.append(value)
    try:
        return json.dumps(cleaned[0], sort_keys=True, allow_nan=False) == json.dumps(cleaned[1], sort_keys=True, allow_nan=False)
    except (ValueError, TypeError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--check', type=Path)
    group.add_argument('--write', type=Path)
    args = parser.parse_args()
    payload = old._json(run())
    if args.check:
        if not _compare_receipt(strict_json.load(args.check), payload):
            raise SystemExit('native census receipt drift')
    else:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(payload, indent=2, sort_keys=True)+'\n')
    print('Finite native-record charge census and continuing inference verified.')


if __name__ == '__main__':
    main()
