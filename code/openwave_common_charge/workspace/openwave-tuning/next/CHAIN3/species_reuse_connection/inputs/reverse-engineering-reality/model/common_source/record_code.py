"""Sharp six-dimensional memory for the complete recorded native seam menu.

Every logical word is the uniform superposition of one perfect matching of
the nonedge graph.  The six matchings partition all 36 nonnative pairs.
The recorded C/E branch Gram matrices are scalar, although the unread work
channel does not preserve this code.  Recovery uses the existing source
compiler and its full-register Gray permutation macro.  All 512 rows,
raw abort records and arbitrary reference columns are retained.

This is active encoded memory under declared coherent control services,
not a physical field identification or a passive native invariant sector.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import math

import numpy as np

from ..observer import native
from ..observer.source_discrimination import circuit, experiment as recorded
from ..observer.source_scalar import compiler
from ..repair_scheduler import coins
from . import native_routing as routing, slow_geometry


MATCHINGS = (
    ((0, 2), (1, 3), (4, 6), (5, 7), (8, 10), (9, 11)),
    ((0, 3), (1, 2), (4, 7), (5, 6), (8, 11), (9, 10)),
    ((0, 4), (1, 6), (2, 8), (3, 11), (5, 10), (7, 9)),
    ((0, 6), (1, 4), (2, 9), (3, 10), (5, 8), (7, 11)),
    ((0, 8), (1, 11), (2, 5), (3, 7), (4, 10), (6, 9)),
    ((0, 9), (1, 10), (2, 7), (3, 5), (4, 8), (6, 11)),
)
SEAMS = tuple(slow_geometry.source().seams())
SUPPORTS = tuple(tuple(tuple(int(p in pair) for p in range(12))
                       for pair in matching) for matching in MATCHINGS)
DIMENSION = 6
SIZE = 512


def _seam(seam):
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(p) is not int for p in seam) or seam not in SEAMS):
        raise ValueError('an ordered native seam from the thirty-seam menu required')
    return seam


def _label(law, seam, coin, raw):
    law, seam = routing._candidate(law), _seam(seam)
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError('an integer native coin zero or one required')
    if type(raw) is not int or not 0 <= raw < 1 << circuit.RECORD_BITS[law]:
        raise ValueError('a raw integer record in the complete native menu required')
    return law, seam, coin, raw


def branch_images(law, seam, coin, raw):
    """Independent endpoint enumeration; all six rows have the same length."""
    law, seam, coin, raw = _label(law, seam, coin, raw)
    images = []
    a, b = seam
    for support in SUPPORTS:
        row = []
        for state in support:
            low, parity = divmod(state[a]+state[b], 2)
            left = low+coin*parity
            current = state[a]-left
            code = {0: 0, 1: 1, -1: 2}[current]
            if law == 'E':
                code += 4*parity
            if code == raw:
                final = list(state)
                final[a], final[b] = left, low+(1-coin)*parity
                row.append(circuit.encode(tuple(final)))
        images.append(tuple(row))
    if (len({len(row) for row in images}) != 1
            or len({i for row in images for i in row}) != sum(map(len, images))):
        raise ArithmeticError('native matching code lost its exact scalar Gram identity')
    return tuple(images)


def branch_weight(law, seam, coin, raw):
    """Exact outcome probability including fair coin and six code supports."""
    return F(len(branch_images(law, seam, coin, raw)[0]), 12)


def _encoder_permutation():
    return routing.complete_permutation(tuple(
        (logical+8*term, circuit.encode(state))
        for logical, support in enumerate(SUPPORTS)
        for term, state in enumerate(support)))


def _branch_permutation(law, seam, coin, raw):
    return routing.complete_permutation(tuple(
        (logical+8*term, address)
        for logical, row in enumerate(branch_images(law, seam, coin, raw))
        for term, address in enumerate(row)))


def _mixing_word(count, inverse=False):
    if type(count) is not int or count not in (1, 4, 5, 6):
        raise ValueError('a declared matching-support size required')
    if type(inverse) is not bool:
        raise ValueError('inverse must be an explicit Boolean')
    vector = [1/math.sqrt(count)]*count+[0.]*(8-count)
    word = compiler.remap(compiler.compile_real_state(vector), (3, 4, 5))
    return compiler.inverse(word) if inverse else word


@lru_cache(maxsize=2)
def _macro_charts(inverse):
    # These are executed source matrices, not ideal controlled-flip targets.
    macro, _ = routing._executed_macro(inverse)
    return tuple(routing._immutable(routing._chart(macro, target)) for target in range(9))


def _apply_permutation(images, state, inverse=False):
    """Apply the complete emitted source schedule, retaining numerical leakage.

    Reuse only the source-executed macro and its literal wire charts.  Each
    invocation is charged by permutation_resources even when numerically
    cached.  No ideal permutation matrix replaces executed amplitudes.
    """
    result = state
    charts = _macro_charts(inverse)
    for target, pattern in routing.permutation_edges(images, inverse=inverse):
        before = routing._conjugation_word(pattern)
        result = compiler.execute(before, 9, result)
        result = charts[target]@result
        result = compiler.execute(compiler.inverse(before), 9, result)
    return result


def _encode(value):
    return _apply_permutation(_encoder_permutation(),
                              compiler.execute(_mixing_word(6), 9, value))


def _finite(value):
    if not np.isfinite(value).all():
        raise ArithmeticError('matching-code execution produced nonfinite amplitudes')
    return routing._immutable(value)


@lru_cache(maxsize=1)
def code_matrix():
    """Execute the encoder on six unknown logical input columns."""
    return _finite(_encode(np.eye(SIZE, DIMENSION, dtype=complex)))


def ideal_code_matrix():
    """Independent declaration used only as a comparison, never preparation."""
    matrix = np.zeros((SIZE, DIMENSION), complex)
    for logical, support in enumerate(SUPPORTS):
        for state in support:
            matrix[circuit.encode(state), logical] = 1/math.sqrt(6)
    return routing._immutable(matrix)


def declared_recovery_matrix66(law, seam, coin, raw):
    """Algebraic recovery on exactly66 states, separate from source execution.

    First invert native movement by a complete endpoint swap.  Each of the
    six disjoint matching blocks then uses a Householder reflection taking
    its normalized surviving subset to the full uniform matching.  All30
    native-edge pair rows remain in the full matrix.  This provides a
    same-carrier extension; source execution uses the independently compiled
    full512 permutation/preparation extension in apply_recovery instead.
    """
    law, seam, coin, raw = _label(law, seam, coin, raw)
    if raw not in circuit.record_values(law, coin):
        raise ValueError('invalid raw native records halt instead of selecting recovery')
    rows = branch_images(law, seam, coin, raw)
    count = len(rows[0])
    if not count:
        return routing._immutable(np.eye(66))
    states = native.states(binary=True)
    index = {state: i for i, state in enumerate(states)}
    source_index = {circuit.encode(state): i for i, state in enumerate(states)}
    value = circuit.record_values(law, coin)[raw]
    current = value if law == 'C' else (value-(2*coin-1)*(abs(value) % 2))//2
    permutation = list(range(66))
    if current:
        a, b = seam
        for i, state in enumerate(states):
            previous = list(state)
            previous[a], previous[b] = previous[b], previous[a]
            permutation[i] = index[tuple(previous)]
    inverse_movement = np.zeros((66, 66), complex)
    inverse_movement[permutation, np.arange(66)] = 1
    reflection = np.eye(66, dtype=complex)
    for logical, outputs in enumerate(rows):
        difference = np.zeros(66, complex)
        difference[[permutation[source_index[x]] for x in outputs]] = 1/math.sqrt(count)
        difference[[index[state] for state in SUPPORTS[logical]]] -= 1/math.sqrt(6)
        norm_squared = float(np.vdot(difference, difference).real)
        if norm_squared:
            reflection -= 2*np.outer(difference, difference.conj())/norm_squared
    return routing._immutable(reflection@inverse_movement)


def encoding(state):
    """Supply a six-level logical state and six blank code bits."""
    routing._original_scalars(state)
    raw = np.asarray(state)
    if (raw.ndim not in (1, 2) or raw.shape[0] != DIMENSION
            or (raw.ndim == 2 and raw.shape[1] < 1)
            or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all()):
        raise ValueError('finite six-component logical amplitudes or reference columns required')
    value = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(value).all():
        raise ValueError('logical amplitudes must remain finite in execution precision')
    return _finite(code_matrix()@value)


def decoding(state):
    """Execute the inverse on every source row; syndrome is not selected."""
    value = routing._state(state)
    value = _apply_permutation(_encoder_permutation(), value, True)
    return _finite(compiler.execute(_mixing_word(6, True), 9, value))


def apply_recovery(law, seam, coin, raw, state):
    """Complete record-selected source unitary; invalid raw records halt."""
    law, seam, coin, raw = _label(law, seam, coin, raw)
    value = routing._state(state)
    if raw not in circuit.record_values(law, coin):
        raise ValueError('invalid raw native records halt instead of selecting recovery')
    count = len(branch_images(law, seam, coin, raw)[0])
    if not count:
        return _finite(value)
    value = _apply_permutation(_branch_permutation(law, seam, coin, raw), value, True)
    value = compiler.execute(_mixing_word(count, True), 9, value)
    return _finite(_encode(value))


def recovery_unitary(law, seam, coin, raw):
    """Source-execute all 512 input columns when a complete matrix is needed."""
    _label(law, seam, coin, raw)
    return apply_recovery(law, seam, coin, raw, np.eye(SIZE, dtype=complex))


def _route_permutation(seam):
    """Route the requested physical ports to source seam (0,1)."""
    a, b = _seam(seam)
    images = list(range(SIZE))
    for state in native.states():
        final = list(state)
        final[0], final[a] = final[a], final[0]
        final[1], final[b] = final[b], final[1]
        images[circuit.encode(state)] = circuit.encode(tuple(final))
    return routing._permutation(tuple(images))


def repair_branches(law, seam, state, *, recover=True):
    """Execute all30 native seams with every raw source record retained."""
    law, seam, value = routing._candidate(law), _seam(seam), routing._state(state)
    if type(recover) is not bool:
        raise ValueError('active recovery must be an explicit Boolean')
    route = _route_permutation(seam)
    value = _apply_permutation(route, value)
    branches = recorded.repair_branches(law, value)
    for branch in branches:
        coin, raw = branch['coin'], branch['record_code']
        amplitude = _apply_permutation(route, branch['amplitude'], True)
        applied = recover and not branch['abort']
        if applied:
            amplitude = apply_recovery(law, seam, coin, raw, amplitude)
        branch.update(amplitude=_finite(amplitude), seam=seam,
                      recovery_applied=applied,
                      code_branch_weight=branch_weight(law, seam, coin, raw),
                      recovery_resources=recovery_resources(law, seam, coin, raw)
                      if applied else dict.fromkeys(routing.COUNT_KEYS, 0))
    return branches


def _counts(*rows):
    return {key: sum(row[key] for row in rows) for key in routing.COUNT_KEYS}


@lru_cache(maxsize=1)
def _macro_resource_counts():
    # Cache only immutable per-invocation counts, never a caller-owned report.
    macro = compiler.resources(routing._canonical_word())
    one_x = compiler.resources([compiler.pulse('X', math.pi/2, [0])])
    return (tuple(macro[key] for key in routing.COUNT_KEYS),
            tuple(one_x[key] for key in routing.COUNT_KEYS))


def _permutation_resources(images):
    edges = routing.permutation_edges(images)
    conjugations = 2*sum(pattern.bit_count() for _, pattern in edges)
    macro, one_x = _macro_resource_counts()
    counts = {key: len(edges)*m+conjugations*x
              for key, m, x in zip(routing.COUNT_KEYS, macro, one_x, strict=True)}
    return {**counts, 'controlled_flip_invocations': len(edges),
            'conjugating_source_x_pulses': conjugations, 'whole_register_dimension': SIZE}


@lru_cache(maxsize=8)
def _mixing_resource_counts(count, inverse=False):
    resources = compiler.resources(_mixing_word(count, inverse))
    return tuple(resources[key] for key in routing.COUNT_KEYS)


def _mixing_resources(count, inverse=False):
    return dict(zip(routing.COUNT_KEYS, _mixing_resource_counts(count, inverse), strict=True))


def encoding_resources():
    counts = _counts(_permutation_resources(_encoder_permutation()),
                     _mixing_resources(6))
    return {**counts, 'supplied_logical_levels': 6, 'supplied_register_qubits': 3,
            'blank_code_qubits': 6, 'new_classical_records': 0,
            'full_register_dimension': SIZE}


def decoding_resources():
    return {**encoding_resources(), 'blank_code_qubits': 0,
            'output_rows_retained': SIZE, 'syndrome_postselection_applied': False}


def recovery_resources(law, seam, coin, raw):
    law, seam, coin, raw = _label(law, seam, coin, raw)
    count = len(branch_images(law, seam, coin, raw)[0])
    admitted = raw in circuit.record_values(law, coin)
    counts = dict.fromkeys(routing.COUNT_KEYS, 0)
    if admitted and count:
        counts = _counts(encoding_resources(),
                         _permutation_resources(_branch_permutation(law, seam, coin, raw)),
                         _mixing_resources(count, True))
    return {**counts, 'code_branch_weight_exact': str(F(count, 12)),
            'raw_abort': not admitted,
            'identity_on_valid_zero_weight_branch': admitted and not count,
            'blank_code_qubits': 0, 'new_classical_records': 0,
            'archive_record_rewrites': 0, 'input_dependent_normalization': False}


@lru_cache(maxsize=2)
def _native_repair_counts(law):
    resources = routing.routed_repair_resources(law, (0, 1))
    return tuple(resources[key] for key in routing.COUNT_KEYS)


def repair_resources(law, seam):
    law, seam = routing._candidate(law), _seam(seam)
    route = _permutation_resources(_route_permutation(seam))
    repair = dict(zip(routing.COUNT_KEYS, _native_repair_counts(law), strict=True))
    return {**_counts(route, route, repair), 'retained_classical_bits': 1+circuit.RECORD_BITS[law],
            'full_register_dimension': SIZE,
            'recovery_cost_is_record_dependent': True}


def _full_carrier_support_certificate():
    """Derive binary support from rank-one effects on the complete78 carrier.

    A rank-one positive effect compressed to k>1 logical dimensions can be
    scalar only if it is zero.  Each double occupation is isolated by a
    retained C coin/current outcome, or an E contrast outcome.  Consequently
    every exact k>1 code is orthogonal to all twelve double occupations.
    The binary carrier used in the next proof is therefore derived.
    """
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    doubles = {p: tuple(2*int(q == p) for q in range(12)) for p in range(12)}
    binary = set(native.states(binary=True))
    if (len(states) != 78 or len(binary) != 66
            or set(states) != binary | set(doubles.values())):
        raise ArithmeticError('complete total-two source carrier changed')
    laws = {}
    for law in ('C', 'E'):
        witnesses = Counter()
        effect_ranks = Counter()
        raw_outcomes = 0
        for seam in SEAMS:
            a, b = seam
            expected = {(coin, raw): [] for coin in (0, 1)
                        for raw in range(1 << circuit.RECORD_BITS[law])}
            # This endpoint oracle does not call native.native or coin_maps.
            for coin in (0, 1):
                for i, state in enumerate(states):
                    low, parity = divmod(state[a]+state[b], 2)
                    final = list(state)
                    final[a], final[b] = low+coin*parity, low+(1-coin)*parity
                    current = state[a]-final[a]
                    raw = {0: 0, 1: 1, -1: 2}[current]
                    if law == 'E':
                        raw += 4*parity
                    expected[coin, raw].append((i, index[tuple(final)]))
            actual = {key: [] for key in expected}
            seen = set()
            for coin, label, mapping in coins.coin_maps(law, seam, binary=False):
                current = label if law == 'C' else (label-(2*coin-1)*(abs(label) % 2))//2
                raw = {0: 0, 1: 1, -1: 2}[current]
                if law == 'E':
                    raw += 4*(abs(label) % 2)
                key = coin, raw
                if key in seen or circuit.record_values(law, coin).get(raw) != label:
                    raise ArithmeticError('full78 source raw record label is inconsistent')
                seen.add(key)
                actual[key] = list(mapping)
                effect_ranks[len(mapping)] += 1
            if any(len(actual[key]) != len(expected[key])
                   or set(actual[key]) != set(expected[key]) for key in expected):
                raise ArithmeticError('full78 source effects differ from the endpoint oracle')
            raw_outcomes += len(expected)
            selected = ((0, 2, b), (1, 1, a)) if law == 'C' else tuple(
                (coin, raw, port) for coin in (0, 1) for raw, port in ((1, a), (2, b)))
            for coin, raw, port in selected:
                inputs = tuple(i for i, _ in actual[coin, raw])
                if inputs != (index[doubles[port]],):
                    raise ArithmeticError('a double-occupation singleton effect is missing')
                witnesses[port] += 1
        expected_count = 5 if law == 'C' else 10
        if witnesses != Counter({port: expected_count for port in range(12)}):
            raise ArithmeticError('singleton records do not exclude every double occupation')
        laws[law] = {
            'complete_raw_outcomes_checked': raw_outcomes,
            'nonzero_effect_rank_counts': {str(k): v for k, v in sorted(effect_ranks.items())},
            'double_occupation_singleton_effects': sum(witnesses.values()),
            'singleton_witness_count_per_port': [witnesses[port] for port in range(12)],
            'singleton_raw_records': 'coin0/raw2 isolates double_b; coin1/raw1 isolates double_a'
            if law == 'C' else 'either coin: raw1/contrast+2 isolates double_a; raw2/contrast-2 isolates double_b',
        }
    return {'initial_total_two_carrier_dimension': 78,
            'double_occupation_basis_states_excluded_for_dimension_greater_than_one': 12,
            'remaining_binary_carrier_dimension': 66,
            'binary_support_derived_not_assumed': True,
            'reason': 'Each retained singleton effect has positive coefficient and rank1. For k>1 its scalar KL compression must vanish, which excludes its double-occupation basis vector.',
            'laws': laws}


def _integer_schedule_action(edges):
    """Execute emitted Gray-neighbor transpositions on every basis address.

    This is an exact schedule oracle, not a source-execution shortcut.  Its
    input is the emitted gate list, rather than a claimed target permutation.
    """
    images, locations = list(range(SIZE)), list(range(SIZE))
    for target, pattern in edges:
        if (type(target) is not int or not 0 <= target < 9
                or type(pattern) is not int or not 0 <= pattern < SIZE
                or pattern & (1 << target)):
            raise ArithmeticError('compiled Gray schedule contains an invalid controlled flip')
        left, right = pattern, pattern | (1 << target)
        a, b = locations[left], locations[right]
        images[a], images[b] = right, left
        locations[left], locations[right] = b, a
    return tuple(images)


def _compiled_layout_certificate():
    """Bind all actual emitted schedules to independent native/logical slots.

    A permutation word implementing its own supplied permutation is not
    enough: a wrong logical-column permutation is still unitary.  Here the
    encoder, native routes and every used inverse recovery schedule are
    compared with separately enumerated endpoint images and logical slots.
    """
    encoder = _encoder_permutation()
    forward = _integer_schedule_action(routing.permutation_edges(encoder))
    inverse = _integer_schedule_action(routing.permutation_edges(encoder, inverse=True))
    if any(inverse[forward[i]] != i for i in range(SIZE)):
        raise ArithmeticError('compiled encoder inverse schedule is inconsistent')
    for logical, matching in enumerate(MATCHINGS):
        for term, pair in enumerate(matching):
            state = tuple(int(p in pair) for p in range(12))
            if forward[logical+8*term] != circuit.encode(state):
                raise ArithmeticError('compiled encoder schedule changes a declared logical slot')
    native_states = native.states()
    source_addresses = {circuit.encode(state) for state in native_states}
    branches_checked = images_checked = 0
    for seam in SEAMS:
        a, b = seam
        route = _route_permutation(seam)
        routed = _integer_schedule_action(routing.permutation_edges(route))
        unrouted = _integer_schedule_action(routing.permutation_edges(route, inverse=True))
        if (any(unrouted[routed[i]] != i for i in range(SIZE))
                or any(routed[i] != i for i in range(SIZE) if i not in source_addresses)):
            raise ArithmeticError('compiled native route loses its inverse or changes padding')
        for state in native_states:
            work = circuit.decode(routed[circuit.encode(state)])
            if (work[0], work[1]) != (state[a], state[b]):
                raise ArithmeticError('compiled native route changes the selected endpoint law')
            for coin in (0, 1):
                total = state[a]+state[b]
                left, right = (total+coin)//2, (total+1-coin)//2
                direct, local = list(state), list(work)
                direct[a], direct[b] = left, right
                local[0], local[1] = left, right
                if unrouted[circuit.encode(tuple(local))] != circuit.encode(tuple(direct)):
                    raise ArithmeticError('compiled native route does not undo the local source repair')
        for law in ('C', 'E'):
            for coin in (0, 1):
                expected = [[[] for _ in range(DIMENSION)]
                            for _ in range(1 << circuit.RECORD_BITS[law])]
                for logical, matching in enumerate(MATCHINGS):
                    for pair in matching:
                        state = [int(p in pair) for p in range(12)]
                        total = state[a]+state[b]
                        left, right = (total+coin)//2, (total+1-coin)//2
                        current = state[a]-left
                        raw = 0 if current == 0 else 1 if current == 1 else 2
                        if law == 'E':
                            raw += 4*(total % 2)
                        state[a], state[b] = left, right
                        expected[raw][logical].append(circuit.encode(tuple(state)))
                for raw, rows in enumerate(expected):
                    if not rows[0]:
                        if any(rows):
                            raise ArithmeticError('independent native-image dimensions are inconsistent')
                        continue
                    permutation = _branch_permutation(law, seam, coin, raw)
                    decoded = _integer_schedule_action(routing.permutation_edges(permutation, inverse=True))
                    for logical, row in enumerate(rows):
                        for term, address in enumerate(row):
                            if decoded[address] != logical+8*term:
                                raise ArithmeticError('compiled recovery schedule changes a native image logical slot')
                            images_checked += 1
                    branches_checked += 1
    if branches_checked != 300 or images_checked != 4320:
        raise ArithmeticError('compiled all-seam recovery certificate omitted an image')
    return {'encoder_forward_and_inverse_checked_on_all512_addresses': True,
            'native_routes_checked_on_all78_states_and_padding': len(SEAMS),
            'inverse_recovery_schedules_checked': branches_checked,
            'independent_native_image_logical_slots_checked': images_checked,
            'all_emitted_integer_schedules_match_independent_native_logical_layout': True,
            'this_integer_certificate_replaces_numerical_source_execution': False}


def exact_certificate():
    """Exact combinatorial witness and the arbitrary-code sharp capacity proof.

    With the record retained the KL conditions are V*D_r V=p_r I, since
    different output labels are orthogonal.  On the full78 carrier, singleton
    C/E effects first exclude all twelve double occupations whenever k>1.
    On the resulting binary66 carrier, for each edge ab C and E both
    contain D_ab=n_a(1-n_b) and D_ba.  Thus compressed n_a-n_b is scalar.
    Connectivity and sum_a n_a=2I force every compressed n_a scalar.  Then
    compressed n_a n_b is scalar and rank at most one, hence zero if k>1.
    The code avoids all30 native-edge basis rows.  Each remaining n_a has
    rank6, and sum_a n_a=2I forces some nonzero scalar compression: k<=6.
    This argument concerns arbitrary complex subspaces and arbitrary CPTP
    recovery with full classical records and unrestricted reference systems.
    """
    full_support = _full_carrier_support_certificate()
    compiled_layout = _compiled_layout_certificate()
    all_pairs = set(combinations(range(12), 2))
    flattened = tuple(pair for row in MATCHINGS for pair in row)
    nonedges = all_pairs-set(SEAMS)
    if len(set(flattened)) != 36 or set(flattened) != nonedges:
        raise ArithmeticError('six matchings no longer partition the source nonedges')
    if any(sorted(p for pair in row for p in pair) != list(range(12)) for row in MATCHINGS):
        raise ArithmeticError('a declared logical support is not a perfect matching')
    neighbors = {v: {b if a == v else a for a, b in SEAMS if v in (a, b)} for v in range(12)}
    reached = {0}
    for _ in range(12):
        reached |= set().union(*(neighbors[v] for v in reached))
    if len(reached) != 12 or {len(neighbors[v]) for v in neighbors} != {5}:
        raise ArithmeticError('full native seam graph lost connectivity or degree5')
    basis = native.states(binary=True)
    index = {state: i for i, state in enumerate(basis)}
    code_indices = tuple({index[state] for state in row} for row in SUPPORTS)
    law_reports = {}
    for law in ('C', 'E'):
        counts = Counter()
        for seam in SEAMS:
            total = F(0)
            for coin, label, mapping in coins.coin_maps(law, seam):
                ins = {i for i, _ in mapping}
                overlaps = tuple(len(ins & row) for row in code_indices)
                if len(set(overlaps)) != 1:
                    raise ArithmeticError('source-enumerated code Gram is nonscalar')
                state = basis[next(iter(ins))]
                _, emitted = native.native(state, seam, coin)
                raw = {0: 0, 1: 1, -1: 2}[emitted[1]]
                if law == 'E':
                    raw += 4*((state[seam[0]]+state[seam[1]]) % 2)
                independent = branch_images(law, seam, coin, raw)
                actual = tuple(tuple(circuit.encode(basis[o]) for i, o in mapping if i in row)
                               for row in code_indices)
                if any(set(a) != set(b) for a, b in zip(actual, independent, strict=True)):
                    raise ArithmeticError('independent endpoint oracle differs from native maps')
                weight = F(overlaps[0], 12)
                counts[str(weight)] += 1
                total += weight
            if total != 1:
                raise ArithmeticError('source branches do not preserve complete probability')
        law_reports[law] = {'nonzero_record_branches': sum(counts.values()),
                            'probability_counts_across_thirty_seams': dict(sorted(counts.items()))}
    return {'logical_dimension': 6, 'native_total_two_carrier_dimension': 78,
            'native_binary_carrier_dimension': 66,
            'full_carrier_binary_support_derivation': full_support,
            'compiled_schedule_certificate': compiled_layout,
            'source_execution_register_dimension': SIZE,
            'matchings': [[list(pair) for pair in row] for row in MATCHINGS],
            'seams_checked': 30, 'native_graph_connected': True,
            'native_graph_degree': 5, 'nonedge_graph_degree': 6,
            'all_exact_branch_gram_identities': True,
            'maximum_exact_record_correctable_code_dimension': 6,
            'upper_bound_native_carrier_dimension': 78,
            'rank_fifteen_exact_record_correctable': False,
            'arbitrary_complex_codes_and_references_in_upper_bound': True,
            'record_conditioned_arbitrary_CPTP_recovery_in_upper_bound': True,
            'complete_unitary_recovery_exists_on_binary66_itself': True,
            'same66_recovery': 'Inverse endpoint swap followed by six disjoint matching-block Householder reflections; no additional quantum carrier in this algebraic extension.',
            'source_compilation_extension': 'The independently compiled permutation/preparation recovery uses the existing full512 source register and may use padding during intermediate controls; no fresh quantum ancilla is added per recovery.',
            'upper_bound_scope': 'An arbitrary fixed code on the full78 total-two carrier precedes state-independent selection with positive probability for every native seam; the seam, coin and C/E value are retained. For dimension greater than one, binary66 support is forced by the singleton effects. Each branch is corrected before the next attempt. Seam-dependent re-encoding before the native measurement is a different contract.',
            'lower_bound_scope': 'Six explicit matching superpositions with full source-word active recovery; coherent controls, code load/export and initialization remain specified services.',
            'laws': law_reports}


def execution_report():
    """Numerical source replay on two declared seams, with all code columns.

    The exact certificate covers every seam.  This separate floating audit
    executes the existing source words for C at (0,1) and E at (2,10),
    including general port routing, every raw record and all512 output rows.
    It does not call numerical replay a physical accuracy certificate.
    """
    compiled_layout = _compiled_layout_certificate()
    expected = ideal_code_matrix()
    executed = code_matrix()
    rows = {}
    for law, seam in (('C', (0, 1)), ('E', (2, 10))):
        branches = repair_branches(law, seam, executed)
        gram = sum((b['amplitude'].conj().T@b['amplitude'] for b in branches),
                   np.zeros((DIMENSION, DIMENSION), complex))
        branch_rows = []
        for branch in branches:
            probability = branch['code_branch_weight']
            target = math.sqrt(float(probability))*expected
            branch_rows.append({
                'coin': branch['coin'], 'raw': branch['record_code'],
                'abort': branch['abort'], 'probability_exact': str(probability),
                'operator_error': float(np.linalg.norm(branch['amplitude']-target, 2)),
                'source_recovery_resources': branch['recovery_resources'],
            })
        rows[law] = {
            'seam': list(seam), 'raw_branches_retained': len(branches),
            'reference_gram_error': float(np.linalg.norm(gram-np.eye(DIMENSION), 2)),
            'maximum_branch_operator_error': max(b['operator_error'] for b in branch_rows),
            'complete_retained_output_operator_error': float(np.linalg.norm(
                np.concatenate([b['amplitude']-math.sqrt(float(b['code_branch_weight']))*expected
                                for b in branches], axis=0), 2)),
            'source_repair_resources': repair_resources(law, seam),
            'branches': branch_rows,
        }
    _, forward = routing._executed_macro(False)
    _, inverse = routing._executed_macro(True)
    encoder_error = float(np.linalg.norm(executed-expected, 2))
    diagnostics = [encoder_error, *(row[key] for row in rows.values() for key in (
        'reference_gram_error', 'maximum_branch_operator_error',
        'complete_retained_output_operator_error'))]
    if any(not math.isfinite(value) or value > 1e-8 for value in diagnostics):
        raise ArithmeticError('complete matching-code source replay exceeds its numerical gate')
    return {'encoder_operator_error': encoder_error,
            'maximum_allowed_numerical_execution_error': 1e-8,
            'encoder_resources': encoding_resources(),
            'executed_code_columns': DIMENSION, 'output_rows_retained': SIZE,
            'all_thirty_seam_numerical_source_replay': False,
            'all_thirty_seam_exact_native_certificate': True,
            'compiled_schedule_certificate': compiled_layout,
            'forward_full512_source_macro_guard': dict(forward),
            'inverse_full512_source_macro_guard': dict(inverse),
            'code_or_padding_projection_applied': False,
            'physical_service_or_roundoff_certificate': False,
            'laws': rows}


def _ideal_primitive(kind, angle):
    """Analytical comparison only; never supplies an executed gate."""
    if kind == 'J':
        return np.diag([np.exp(1j*angle)]*3+[1.])
    identity = np.eye(2, dtype=complex)
    generators = {
        'X': np.asarray([[0, 1j], [1j, 0]]),
        'Y': np.asarray([[0, 1], [-1, 0]], complex),
        'Z': np.diag([1j, -1j]),
        'I': 1j*identity,
    }
    generators['H'] = (generators['X']+generators['Z'])/math.sqrt(2)
    return math.cos(angle)*identity+math.sin(angle)*generators[kind]


@lru_cache(maxsize=2048)
def _primitive_error(kind, angle):
    dimension = 4 if kind == 'J' else 2
    actual = compiler.native_pulse(kind, angle)[:dimension, :dimension]
    return float(np.linalg.norm(actual-_ideal_primitive(kind, angle), 2))


def _word_log_error(word):
    return math.fsum(math.log1p(_primitive_error(e['kind'], e['angle'])) for e in word)


def _literal_mixing_words():
    """Validate and freeze the actual emitted words before numerical caches."""
    rows = []
    for count in (1, 4, 5, 6):
        words = []
        for inverse in (False, True):
            word = compiler._validate(_mixing_word(count, inverse), 9)
            if any(target not in (3, 4, 5) for event in word for target in event['targets']):
                raise ArithmeticError('literal mixing word acts outside the declared syndrome wires')
            words.append(word)
        if words[1] != compiler.inverse(words[0]):
            raise ArithmeticError('literal inverse mixing word is not the reversed negative-angle forward word')
        rows.extend(tuple((event['kind'], event['angle'], tuple(event['targets']))
                          for event in word) for word in words)
    return tuple(rows)


def _mixing_synthesis_report(count, forward, inverse):
    """Separate angle synthesis from source realization of requested angles.

    An exactly unitary requested word U can be changed to a unitary T with
    T|0> equal to the declared uniform vector by a rotation on their span.
    Its operator change is at most twice the first-column discrepancy.  The
    same bound applies to U* versus T*.  This comparison supplies no gate
    to the execution path.  Values are numerical diagnostics, not intervals.
    """
    target = np.asarray([1/math.sqrt(count)]*count+[0.]*(8-count), complex)
    # Reduce the literal, already validated 9-bit targets to their actual
    # three-bit active block.  Recompiling an unbound reference word would
    # miss a wrong but unitary move onto the logical wires.
    matrices = []
    for word in (forward, inverse):
        requested = np.eye(8, dtype=complex)
        for kind, angle, targets in word:
            requested = compiler._apply(requested, _ideal_primitive(kind, angle),
                                        [wire-3 for wire in targets], 3)
        matrices.append(requested)
    requested, requested_inverse = matrices
    column_error = float(np.linalg.norm(requested[:, 0]-target))
    correction = 2*column_error
    residual = float(np.linalg.norm(requested.conj().T@requested-np.eye(8), 2))
    inverse_error = float(np.linalg.norm(requested_inverse-requested.conj().T, 2))
    inverse_column_error = float(np.linalg.norm(requested_inverse@target-np.eye(8)[:, 0]))
    if any(not math.isfinite(value) or value > 1e-10
           for value in (correction, residual, inverse_error, inverse_column_error)):
        raise ArithmeticError('uniform-support angle synthesis exceeds its numerical gate')
    return {'support_count': count,
            'requested_word_uniform_column_error': column_error,
            'uniform_target_unitary_completion_error_upper': correction,
            'requested_word_unitarity_residual': residual,
            'literal_syndrome_wire_support_validated': True,
            'literal_inverse_word_is_exact_reversed_negative_angle_word': True,
            'requested_inverse_adjoint_error': inverse_error,
            'requested_inverse_uniform_to_zero_error': inverse_column_error,
            'forward_and_inverse_use_the_same_target_completion_bound': True}


def _native_dilation_error(law, coin):
    """Compare every local source column to its independent integer dilation."""
    dimension = 1 << (3+circuit.RECORD_BITS[law])
    lower = list(range(32))
    for left, right in circuit.REPAIR_SWAPS['C', coin]:
        lower[left], lower[right] = lower[right], lower[left]
    images = []
    for address in range(dimension):
        final = (address & ~31)+lower[address & 31]
        if law == 'E' and (final >> coin) & 1:
            final ^= 32
        images.append(final)
    expected = np.zeros((dimension, dimension), complex)
    expected[images, np.arange(dimension)] = 1
    return float(np.linalg.norm(recorded._repair(law, coin)-expected, 2))


@lru_cache(maxsize=1)
def _source_composition_report(mixing_words):
    """Uniform numerical error budget from every all30 source schedule.

    If A_j differs from an exactly unitary U_j by delta_j, telescoping
    bounds the full product error by product_j(1+delta_j)-1.  Every emitted
    macro and conjugation is counted, with no ideal target used in execution.
    Each macro error includes both its stabilized target residual and the
    full-register polar change, so every invocation covers the raw executed
    source word as well as the stabilized numerical replay.
    The maximum over orthogonal classical control blocks costs a maximum,
    not a sum over their probabilities.  The calculated residuals and sums
    use floating arithmetic: this is a numerical composition diagnostic,
    not a directed-rounding, hardware or physical-service certificate.
    """
    exact_flip = np.eye(SIZE, dtype=complex)
    exact_flip[[0, 256]] = exact_flip[[256, 0]]
    macro_errors, stabilized_errors, projection_changes = {}, {}, {}
    for inverse in (False, True):
        actual, guard = routing._executed_macro(inverse)
        stabilized_errors[inverse] = float(np.linalg.norm(actual-exact_flip, 2))
        # compiler.stabilize_executed_unitary evaluates this change with the
        # spectral norm on all512 rows and columns.  The triangle inequality
        # charges it once per macro, not once per entire composed history.
        projection_changes[inverse] = guard['projection_change']
        macro_errors[inverse] = stabilized_errors[inverse]+projection_changes[inverse]
        if (guard['dimension'] != SIZE or not math.isfinite(macro_errors[inverse])
                or not 0 <= projection_changes[inverse] <= guard['maximum_projection_change']):
            raise ArithmeticError('full-register macro stabilization guard is invalid')
    exact_ix = np.asarray([[0, 1j], [1j, 0]])
    xerror = max(float(np.linalg.norm(compiler.native_pulse('X', sign*math.pi/2)[:2, :2]
                                      -sign*exact_ix, 2)) for sign in (1, -1))
    synthesis = {count: _mixing_synthesis_report(count, *mixing_words[2*i:2*i+2])
                 for i, count in enumerate((1, 4, 5, 6))}

    mixing_logs = {(count, inverse): (_word_log_error([
                       {'kind': kind, 'angle': angle, 'targets': list(targets)}
                       for kind, angle, targets in mixing_words[2*i+int(inverse)]])
                   +math.log1p(synthesis[count]['uniform_target_unitary_completion_error_upper']))
                   for i, count in enumerate((1, 4, 5, 6)) for inverse in (False, True)}

    def mixing_log(count, inverse=False):
        return mixing_logs[count, inverse]

    def permutation_log(images, inverse=False):
        edges = routing.permutation_edges(images, inverse=inverse)
        xcount = 2*sum(pattern.bit_count() for _, pattern in edges)
        return len(edges)*math.log1p(macro_errors[inverse])+xcount*math.log1p(xerror)

    encoder_log = permutation_log(_encoder_permutation())+mixing_log(6)
    decoder_log = permutation_log(_encoder_permutation(), True)+mixing_log(6, True)
    encoder_resources = encoding_resources()
    exact_coin = np.asarray([[1., -1.], [1., 1.]])/math.sqrt(2)
    coin_error = float(np.linalg.norm(compiler.native_pulse('Y', -math.pi/4)[:2, :2]-exact_coin, 2))
    coin_log = math.log1p(coin_error)
    rows = {}
    for law in ('C', 'E'):
        native_errors = [_native_dilation_error(law, coin) for coin in (0, 1)]
        native_log = math.log1p(max(native_errors))
        worst_log = 0.
        worst_resources = dict.fromkeys(routing.COUNT_KEYS, 0)
        worst_recovery_pulses = 0
        for seam in SEAMS:
            route = _route_permutation(seam)
            route_log = permutation_log(route)+permutation_log(route, True)
            largest_recovery_log = 0.
            largest_recovery_resources = dict.fromkeys(routing.COUNT_KEYS, 0)
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    count = len(branch_images(law, seam, coin, raw)[0])
                    if not count:
                        continue
                    current_log = (permutation_log(_branch_permutation(law, seam, coin, raw), True)
                                   +mixing_log(count, True)+encoder_log)
                    largest_recovery_log = max(largest_recovery_log, current_log)
                    resources = recovery_resources(law, seam, coin, raw)
                    for key in routing.COUNT_KEYS:
                        largest_recovery_resources[key] = max(largest_recovery_resources[key], resources[key])
            worst_log = max(worst_log, route_log+coin_log+native_log+largest_recovery_log)
            costs = _counts(repair_resources(law, seam), largest_recovery_resources)
            for key in routing.COUNT_KEYS:
                worst_resources[key] = max(worst_resources[key], costs[key])
            worst_recovery_pulses = max(worst_recovery_pulses, largest_recovery_resources['native_pulses'])
        rows[law] = {
            'maximum_corrected_attempt_operator_error': math.expm1(worst_log),
            'maximum_encoded_corrected_attempt_operator_error': math.expm1(encoder_log+worst_log),
            'maximum_corrected_attempt_log_error': worst_log,
            'maximum_corrected_attempt_resources': worst_resources,
            'maximum_recovery_native_pulses': worst_recovery_pulses,
            'complete_local_native_dilation_operator_errors_by_coin': native_errors,
            'all_thirty_seam_schedules_and_raw_recoveries_counted': True,
        }
    if any(not math.isfinite(row['maximum_encoded_corrected_attempt_operator_error'])
           or row['maximum_encoded_corrected_attempt_operator_error'] > 1e-6
           for row in rows.values()):
        raise ArithmeticError('all-seam source composition exceeds its numerical gate')
    return {
        'method': 'product(1+factor_operator_error)-1; maximum over retained classical control blocks',
        'forward_full512_macro_operator_error': macro_errors[False],
        'inverse_full512_macro_operator_error': macro_errors[True],
        'forward_full512_stabilized_macro_operator_error': stabilized_errors[False],
        'inverse_full512_stabilized_macro_operator_error': stabilized_errors[True],
        'forward_full512_macro_polar_change': projection_changes[False],
        'inverse_full512_macro_polar_change': projection_changes[True],
        'raw_and_stabilized_macro_conventions_both_covered': True,
        'macro_error_convention': 'Stabilized-to-exact operator residual plus raw-to-stabilized spectral norm change, charged at every emitted macro invocation.',
        'source_X_conjugation_operator_error': xerror,
        'fair_coin_preparation_operator_error': coin_error,
        'mixing_angle_synthesis_diagnostics': {str(k): v for k, v in synthesis.items()},
        'encoder_operator_error_bound': math.expm1(encoder_log),
        'encoder_log_error': encoder_log,
        'decoder_operator_error_bound': math.expm1(decoder_log),
        'decoder_log_error': decoder_log,
        'encoder_resources': encoder_resources,
        'numerical_source_composition_bound_not_interval_or_hardware_certificate': True,
        'unrestricted_reference_extension_has_same_operator_bound': True,
        'complete_raw_classical_control_blocks_retained': True,
        'laws': rows,
    }


def source_composition_report():
    """Return a fresh report; caller mutations cannot poison cached bounds."""
    compiled_layout = _compiled_layout_certificate()
    mixing_words = _literal_mixing_words()
    for i, count in enumerate((1, 4, 5, 6)):
        _mixing_synthesis_report(count, *mixing_words[2*i:2*i+2])
    # Word signatures key the numerical cache.  A changed nominal angle or
    # ordering cannot inherit a previously accepted synthesis diagnostic.
    report = deepcopy(_source_composition_report(mixing_words))
    report['compiled_schedule_certificate'] = compiled_layout
    return report


source_execution_report = execution_report
