"""A complete finite native training tape on one continuing source carrier.

Only committed coin/record counts cross the inference interface. The simulator
privately uses endpoint weights to draw each coin by integer rejection sampling;
it never replaces a realized count by its expectation. Native basis transitions
are read afresh and checked against both endpoint arithmetic and the local
source-dilation permutation. Repeated fixed-seam events are then evaluated in
vectorized chunks using those checked transitions, without resetting the source.

This is an exact basis-trajectory replay, not a coherent execution of its many
compiled source pulses. The separate complete-effect certificate establishes
the input-independent coin law for arbitrary quantum inputs and untouched R.
The PRNG seed is reproducibility data, never a probability guarantee.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from . import local_bias as local, mixing
from ..observer import native, record_recovery as records
from ..observer.source_discrimination import circuit
from ..common_source import native_routing as routing, record_code


ROOT = Path(__file__).resolve().parents[2]
SAMPLES_PER_SEAM = 1 << 20
DEFAULT_SEED = 20261029
CHUNK_SIZE = 1 << 16
LAWS = ('C', 'E')
STREAM_PREFIX = b'OPH-native-endpoint-training-tape-v1\0'


def _integer(value, lower, upper, name):
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError(f'{name} must be an original integer in [{lower},{upper}]')
    return value


def _initial(value):
    if value is None:
        value = (1, 1)+(0,)*10
    if (type(value) is not tuple or len(value) != 12
            or any(type(x) is not int or x not in (0, 1) for x in value)
            or sum(value) != 2):
        raise ValueError('an original binary total-two source configuration required')
    return value


def _immutable(value):
    return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)


def transition_tables():
    """Fresh complete H78 tables, independently checked before any fast replay.

    For each law, seam, coin and original state there is exactly one output
    and one raw address. Within a fixed raw address the map is injective.
    Therefore sum_raw K*(coin,raw) K(coin,raw)=p_coin I, including arbitrary R.
    This is the reason inference may use the fresh coin record without reading
    the source. It is not an assertion of source/archive independence.
    """
    states = tuple(native.states())
    seams = tuple(records.seams())
    if (len(states) != 78 or len(set(states)) != 78
            or len(seams) != 30 or len(set(seams)) != 30
            or any(type(e) is not tuple or len(e) != 2 or e[0] >= e[1] for e in seams)
            or circuit.RECORD_BITS != {'C': 2, 'E': 3}
            or circuit.CURRENT_CODES != {0: 0, 1: 1, -1: 2}):
        raise ArithmeticError('native state, seam or raw-record alphabet changed')
    index = {x: i for i, x in enumerate(states)}
    next_state = np.empty((2, 30, 2, 78), dtype=np.uint8)
    raw_record = np.empty_like(next_state)
    raw_domains = {}
    for li, law in enumerate(LAWS):
        domains = np.zeros((30, 2, 1 << circuit.RECORD_BITS[law]), dtype=np.uint16)
        for ei, (a, b) in enumerate(seams):
            for coin in (0, 1):
                inverse = circuit.record_values(law, coin)
                images = [set() for _ in range(1 << circuit.RECORD_BITS[law])]
                for i, state in enumerate(states):
                    output, record = records.recorded_step(law, state, (a, b), coin)
                    total = state[a]+state[b]
                    expected = list(state)
                    expected[a], expected[b] = (total+coin)//2, (total+1-coin)//2
                    current = state[a]-expected[a]
                    semantic = current if law == 'C' else state[a]-state[b]
                    raw = {0: 0, 1: 1, -1: 2}[current]
                    if law == 'E':
                        raw += 4*(total % 2)
                    # Independently apply the actual local source permutation:
                    # C's disjoint swaps, followed by E's parity-copy CNOT.
                    source_address = circuit.SEAM_STATES.index((state[a], state[b]))
                    for left, right in circuit.REPAIR_SWAPS['C', coin]:
                        if source_address in (left, right):
                            source_address = right if source_address == left else left
                    if law == 'E' and (source_address >> coin) & 1:
                        source_address ^= 1 << 5
                    local_output = source_address & 7
                    if (tuple(expected) != output or record != (a, b, coin, semantic)
                            or inverse.get(raw) != semantic
                            or source_address >> 3 != raw
                            or local_output >= len(circuit.SEAM_STATES)
                            or circuit.SEAM_STATES[local_output] != (expected[a], expected[b])):
                        raise ArithmeticError('literal native/source record differs from endpoint oracle')
                    j = index[output]
                    if j in images[raw]:
                        raise ArithmeticError('complete native raw block is not a partial injection')
                    images[raw].add(j)
                    next_state[li, ei, coin, i] = j
                    raw_record[li, ei, coin, i] = raw
                    domains[ei, coin, raw] += 1
                if sum(map(len, images)) != 78:
                    raise ArithmeticError('complete raw effects do not sum to the coin identity')
        raw_domains[law] = domains.tolist()
    if not np.array_equal(next_state[0], next_state[1]):
        raise ArithmeticError('C and E must have the same endpoint population map')
    # The fast repeated-seam reduction is checked on every original H78 input,
    # including first-event collisions. No target-state substitution is used.
    for ei in range(30):
        for initial in range(78):
            settled = next_state[0, ei, :, initial]
            for li in range(2):
                for previous in (0, 1):
                    for coin in (0, 1):
                        if next_state[li, ei, coin, settled[previous]] != settled[coin]:
                            raise ArithmeticError('fixed-seam two-state reduction is invalid')
    digest = hashlib.sha256(next_state.tobytes()+raw_record.tobytes()).hexdigest()
    certificate = {
        'native_carrier_dimension': 78, 'binary_replay_dimension': 66,
        'seams': [list(e) for e in seams], 'checked_basis_columns': 2*30*2*78,
        'complete_raw_slots': 30*(2*4+2*8),
        'raw_domain_sizes_H78': raw_domains,
        'transition_table_sha256': digest,
        'coin_effect_identity': 'sum_raw K(coin,raw)^* K(coin,raw)=p_coin I_78',
        'quantum_scope': 'every initial density and untouched reference; fresh coin independent of past records',
        'not_claimed': 'source, past records, and untouched reference need not be independent',
        'fixed_seam_reduction_checks': 30*78*2*2*2,
    }
    return {'states': states, 'seams': seams, 'next': _immutable(next_state),
            'raw': _immutable(raw_record), 'certificate': certificate}


def _events(tables, edge_index, initial, coins, ordinal):
    """Evaluate every event from fresh native tables; no synthetic histogram."""
    if (type(coins) is not np.ndarray or coins.dtype != np.uint8 or coins.ndim != 1
            or not len(coins) or np.any(coins > 1)):
        raise ValueError('a nonempty original uint8 vector of binary coins required')
    _integer(edge_index, 0, 29, 'edge index')
    _integer(initial, 0, 77, 'H78 state index')
    _integer(ordinal, 0, 30*SAMPLES_PER_SEAM, 'event ordinal')
    settled = tables['next'][0, edge_index, :, initial]
    before = np.empty(len(coins), dtype=np.uint8)
    before[0] = initial
    before[1:] = settled[coins[:-1]]
    laws = ((np.arange(len(coins), dtype=np.uint64)+ordinal) & 1).astype(np.uint8)
    raw = tables['raw'][laws, edge_index, coins, before]
    # Ten-bit encoding, stored as little-endian uint16 only for streaming.
    slots = (laws.astype(np.uint16)*480+16*edge_index+8*coins+raw).astype('<u2')
    return slots, laws, raw, int(settled[coins[-1]])


def native_word_costs(edges=None):
    """Fresh exact pulse costs for the declared routed endpoint service.

    The biased coin still uses one Y pulse, with a different angle. Routing
    words and local native repair lengths are counted literally. These costs
    describe that supplied implementation; the acquisition trajectory executes
    exact basis transitions and does not numerically replay these pulse words.
    Explicit requests may include any ascending distinct port pair. This only
    counts its routed service; it does not add that pair to the native scheduler.
    """
    native_seams = tuple(records.seams())
    seams = native_seams if edges is None else edges
    if (type(seams) is not tuple or not 1 <= len(seams) <= 66
            or any(type(edge) is not tuple or len(edge) != 2
                   or any(type(port) is not int for port in edge)
                   or not 0 <= edge[0] < edge[1] < 12 for edge in seams)
            or len(set(seams)) != len(seams)):
        raise ValueError('distinct original ascending port-pair tuples required')
    macro = routing._canonical_word()
    macro_count = len(macro)
    local_lengths = {law: tuple(len(circuit.repair_word(law, c)) for c in (0, 1)) for law in LAWS}
    if any(row[0] != row[1] for row in local_lengths.values()):
        raise ArithmeticError('native coin branches no longer have equal pulse counts')
    result = {}
    for edge in seams:
        expected_images = list(range(512))
        for state in native.states():
            expected = list(state)
            a, b = edge
            expected[0], expected[a] = expected[a], expected[0]
            expected[1], expected[b] = expected[b], expected[1]
            expected_images[circuit.encode(state)] = circuit.encode(expected)
        # The native menu uses the existing emitted route; an explicitly
        # requested nonnative pair uses the same declared full permutation
        # grammar, without claiming that it belongs to the native scheduler.
        images = (record_code._route_permutation(edge) if edge in native_seams
                  else routing._permutation(tuple(expected_images)))
        if tuple(images) != tuple(expected_images):
            raise ArithmeticError('compiled source route changed physical ports or padding')
        gray_edges = routing.permutation_edges(images)
        route_count = len(gray_edges)*macro_count+2*sum(pattern.bit_count() for _, pattern in gray_edges)
        result[edge] = {law: 1+local_lengths[law][0]+2*route_count for law in LAWS}
    return result


def resource_counts(samples_per_seam, *, chunk_size=CHUNK_SIZE):
    """Exact ideal services plus separately labelled literal source-word costs."""
    n = _integer(samples_per_seam, 1, SAMPLES_PER_SEAM, 'samples per seam')
    chunk = _integer(chunk_size, 1, SAMPLES_PER_SEAM, 'chunk size')
    rows = []; total_pulses = 0
    for ei, (edge, costs) in enumerate(native_word_costs().items()):
        c_events = (n+1)//2 if (ei*n) % 2 == 0 else n//2
        counts = (c_events, n-c_events)
        per_event = tuple(costs[law] for law in LAWS)
        total_pulses += sum(c*p for c, p in zip(counts, per_event, strict=True))
        rows.append({'seam': list(edge), 'C_events': counts[0], 'E_events': counts[1],
                     'compiled_C_pulses_per_event': per_event[0], 'compiled_E_pulses_per_event': per_event[1]})
    total = 30*n
    return {
        'native_endpoint_invocations': total, 'ideal_biased_coin_invocations': total,
        'deterministic_seam_selector_calls': total, 'deterministic_law_selector_calls': total,
        'bounded_uniform_integer_draws': total,
        'host_vectorized_rng_calls': 30*((n+chunk-1)//chunk),
        'host_random_bits': 'variable rejection-sampling cost, not a fixed bit budget',
        'native_input_preparations': 1, 'source_resets_during_training': 0,
        'blank_coin_register_initializations': total, 'blank_raw_register_initializations': total,
        'committed_coin_register_reads': total, 'committed_raw_register_reads': total,
        'record_encoding_bits': 10, 'record_writes': total,
        'FIFO_capacity': 1, 'FIFO_bits': 10, 'exported_records': total-1,
        'exported_record_bits': 10*(total-1), 'total_record_bits': 10*total,
        'source_address_bits': 9, 'peak_native_source_and_raw_qubits': 12,
        'six_level_processor_additional_for_compiled_realization': True,
        'event_counter_bits': total.bit_length(),
        'per_seam_coin_counter_bits': n.bit_length(),
        'all30_coin_counter_bits': 30*n.bit_length(),
        'diagnostic_histograms_and_stream_hash': 'host replay workspace; complete physical record storage counted above',
        'compiled_word_costs': {'native_pulses': total_pulses,
            'proper_code_load_calls': total_pulses, 'proper_code_export_calls': total_pulses,
            'per_seam': rows,
            'scope': 'literal supplied source-word counts only; not a numerical pulse replay or physical cost/error certificate'},
        'supplied_services': 'biased coin, input preparation, record read/commit/export, selectors and storage need physical implementations',
    }


def source_bindings():
    paths = ('model/local_repair_family/acquisition.py', 'model/local_repair_family/test_acquisition.py',
             'model/local_repair_family/local_bias.py', 'model/local_repair_family/mixing.py',
             'model/observer/native.py', 'model/observer/record_recovery.py',
             'model/observer/source_discrimination/circuit.py', 'model/common_source/native_routing.py',
             'model/common_source/record_code.py', 'model/observer/source_scalar/compiler.py',
             'model/common_source/slow_geometry.py', 'evidence/native_geometric_source_20260925/spec.json',
             'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
             'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def acquire(private_parameters, *, samples_per_seam=SAMPLES_PER_SEAM,
            seed=DEFAULT_SEED, initial_state=None, chunk_size=CHUNK_SIZE):
    """Return (JSON-safe replay report, inference-only exact count mapping).

    The second return is exactly {ascending_seam: (trials, coin1_count)}.
    It contains neither weights, source loads, initial/final state nor seed.
    The first return is private simulator evidence and complete-tape custody.
    Replay it with the same supplied parameters; hashes are never accepted as
    evidence without the owned fresh reconstruction.
    """
    key = mixing._parameters(private_parameters)
    ports = key[0]
    n = _integer(samples_per_seam, 1, SAMPLES_PER_SEAM, 'samples per seam')
    seed = _integer(seed, 0, (1 << 64)-1, 'seed')
    chunk = _integer(chunk_size, 1, SAMPLES_PER_SEAM, 'chunk size')
    initial = _initial(initial_state)
    tables = transition_tables()
    state_index = tables['states'].index(initial)
    binary = tuple(native.states(binary=True))
    generator = np.random.Generator(np.random.PCG64(seed))
    stream = hashlib.sha256(STREAM_PREFIX)
    counts = {}; blocks = []; ordinal = 0; fifo_slot = None
    for ei, edge in enumerate(tables['seams']):
        a, b = edge
        denominator, numerator = ports[a]+ports[b], ports[a]
        block_start = state_index
        raw_counts = {'C': np.zeros((2, 4), dtype=np.int64), 'E': np.zeros((2, 8), dtype=np.int64)}
        ones = 0; block_hash = hashlib.sha256(STREAM_PREFIX)
        for offset in range(0, n, chunk):
            size = min(chunk, n-offset)
            # Generator.integers uses bounded integer sampling, not a floating
            # Bernoulli grid. Explicit dtype fixes the PRNG consumption format.
            draws = generator.integers(0, denominator, size=size, dtype=np.uint64)
            coins = (draws < numerator).astype(np.uint8)
            slots, laws, raw, state_index = _events(tables, ei, state_index, coins, ordinal)
            tape_bytes = slots.tobytes()
            stream.update(tape_bytes); block_hash.update(tape_bytes)
            ones += int(np.count_nonzero(coins))
            for li, law in enumerate(LAWS):
                choose = laws == li
                width = 1 << circuit.RECORD_BITS[law]
                addresses = coins[choose].astype(np.int64)*width+raw[choose]
                raw_counts[law] += np.bincount(addresses, minlength=2*width).reshape(2, width)
            fifo_slot = int(slots[-1]); ordinal += size
        if sum(int(row.sum()) for row in raw_counts.values()) != n:
            raise ArithmeticError('complete emitted raw counts lost an event')
        if sum(int(row[1].sum()) for row in raw_counts.values()) != ones:
            raise ArithmeticError('coin and joint raw-record marginals disagree')
        counts[edge] = (n, ones)
        blocks.append({'seam': list(edge), 'trials': n, 'coin1_count': ones,
                       'first_event_ordinal': ei*n, 'start_state_index_H78': block_start,
                       'terminal_state_index_H78': state_index,
                       'raw_counts': {law: value.tolist() for law, value in raw_counts.items()},
                       'stream_sha256': block_hash.hexdigest()})
    terminal = tables['states'][state_index]
    if terminal not in binary or ordinal != 30*n:
        raise ArithmeticError('native training left its carrier or fixed horizon')
    report = {
        'schema': 'oph.local_repair_family.acquisition.v1',
        'protocol': {'samples_per_seam': n, 'total_events': ordinal,
            'seam_order': [list(e) for e in tables['seams']],
            'law_schedule': 'global event ordinal even:C, odd:E; committed before each fresh coin',
            'fixed_training_schedule_overrides_continuation_selector': True,
            'record_encoding': 'uint10 slot=480*law_index+16*seam_index+8*coin+raw; stream uses little-endian uint16',
            'stream_prefix_hex': STREAM_PREFIX.hex(),
            'RNG': 'NumPy PCG64 / Generator.integers uint64 exact bounded-integer draws',
            'seed': seed, 'seed_supplies_statistical_guarantee': False},
        'table_certificate': tables['certificate'], 'blocks': blocks,
        'stream_sha256': stream.hexdigest(), 'retained_FIFO_slot': fifo_slot,
        'private_replay': {'port_weights': list(ports), 'initial_state': list(initial),
            'terminal_state': list(terminal), 'terminal_state_index': binary.index(terminal),
            'carrier_is_not_reset_between_seam_blocks': True,
            'basis_trajectory_is_not_arbitrary_quantum_execution': True},
        'resources': resource_counts(n, chunk_size=chunk),
        'scope': {'inference_input': 'only the separate seam -> (trials,coin1_count) mapping',
            'all_raw_slots_retained_in_histograms_including_zero': True,
            'complete_tape_reconstructible_by_seeded_native_replay': True,
            'unknown_quantum_input_coin_law': 'certified by complete raw effects, not by the basis trajectory',
            'physical_time_energy_or_precision_identified': False},
        'source_sha256': source_bindings(),
    }
    return report, counts
