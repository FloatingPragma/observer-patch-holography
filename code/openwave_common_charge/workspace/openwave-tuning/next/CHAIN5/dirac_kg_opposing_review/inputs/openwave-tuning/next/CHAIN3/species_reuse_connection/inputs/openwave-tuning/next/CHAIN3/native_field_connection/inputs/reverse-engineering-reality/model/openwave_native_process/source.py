"""Literal C/E repair, complete records and a finite causal event process.

The source is the unchanged native_program interpreter. This module supplies
an exact endpoint interface and event-addressed observer records. It accepts
no Hamiltonian, target spectrum, physical clock or hidden-state scheduler.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from model.native_program import interpreter as native
from model.federated_observer import source as federation
from model.observer import native as endpoint


@dataclass(frozen=True)
class LocalBranch:
    """Kraus amplitude sqrt(weight) times one endpoint partial injection.

    Endpoint address x has occupations (x&1, (x>>1)&1). Record is the
    five-tuple (kind, law, coin, raw, label); current is a decoded record.
    All 25 slots exist, including the zero hold and register complements.
    """
    record: tuple
    weight: F
    pairs: tuple[tuple[int, int], ...]
    current: int | None


def graph(level=0):
    """Original immutable port_pair federation: 240,960,3840 ports."""
    return federation.federation(level)


def fixed_command(seam, law='mixture', repeat=1):
    return {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(seam)},
            'law': law, 'repeat': repeat}


def admitted_program(level=0, total=1, commands=None, capacity=3,
                     overflow='export_oldest'):
    """Fair endpoint coins, C/E weight1/2, repair probability1, literal menu."""
    return native.program(federation={'level': level, 'total': total},
                          commands=commands, capacity=capacity, overflow=overflow)


def admitted_local_program(shell='H66', commands=None, capacity=3,
                           overflow='export_oldest'):
    return native.program(shell=shell, biased=False, commands=commands,
                          capacity=capacity, overflow=overflow)


@lru_cache(None)
def raw_local_branches(law='mixture'):
    """The actual endpoint rule and interpreter's complete raw record alphabet."""
    if law not in ('C', 'E', 'mixture'):
        raise ValueError('C, E or the fixed fair mixture required')
    result = [LocalBranch(('hold', None, None, None, None), F(0),
                          tuple((x, x) for x in range(4)), None)]
    for kind in ('C', 'E'):
        weight = F(1, 4) if law == 'mixture' else F(kind == law, 2)
        for coin in (0, 1):
            decoding = native.circuit.record_values(kind, coin)
            for raw in range(1 << native.circuit.RECORD_BITS[kind]):
                label = decoding.get(raw)
                current = (None if label is None else label if kind == 'C' else
                           (label-(2*coin-1)*(abs(label) % 2))//2)
                pairs = []
                for x in range(4):
                    loads = (x & 1, (x >> 1) & 1)
                    output, (_, actual_current) = endpoint.native(loads, (0, 1), coin)
                    actual_label = actual_current if kind == 'C' else loads[0]-loads[1]
                    if label is not None and label == actual_label:
                        pairs.append((x, output[0]+2*output[1]))
                result.append(LocalBranch(('repair', kind, coin, raw, label), weight,
                                          tuple(pairs), current))
    if len(result) != 25:
        raise ArithmeticError('complete native raw alphabet changed')
    return tuple(result)


def apply_branch(mask, seam, branch):
    """Unscaled partial map, with all other binary occupations untouched."""
    if type(mask) is not int or mask < 0:
        raise ValueError('nonnegative global occupation mask required')
    a, b = seam
    if type(a) is not int or type(b) is not int or not 0 <= a < b:
        raise ValueError('ordered distinct nonnegative endpoints required')
    x = ((mask >> a) & 1)+2*((mask >> b) & 1)
    image = dict(branch.pairs).get(x)
    if image is None:
        return None
    return (mask & ~((1 << a) | (1 << b))) | ((image & 1) << a) | ((image >> 1) << b)


def complete_local_events(occupied, seam, law='mixture'):
    """All25 classical raw records, including each probability-zero outcome.

    occupied is an integer mask. Zero probability and undefined partial-map
    image are separate fields. No successful outcome is normalized.
    """
    return tuple({'slot': k, 'record': branch.record, 'weight': branch.weight,
                  'image': apply_branch(occupied, seam, branch),
                  'probability': branch.weight if apply_branch(occupied, seam, branch) is not None else F(0),
                  'current': branch.current}
                 for k, branch in enumerate(raw_local_branches(law)))


def branch_column(left, right, seam, law='mixture'):
    """Complete matrix-unit instrument. Values are exact CP coefficients."""
    result = []
    for branch in raw_local_branches(law):
        a, b = apply_branch(left, seam, branch), apply_branch(right, seam, branch)
        result.append((branch.record, branch.weight if a is not None and b is not None else F(0),
                       None if a is None or b is None else (a, b)))
    return tuple(result)


def unread_column(left, right, seam, law='mixture'):
    result = defaultdict(F)
    for _, weight, pair in branch_column(left, right, seam, law):
        if weight:
            result[pair] += weight
    return {key: value for key, value in result.items() if value}


def causal_edges(seams):
    """Precedence for every pair of events that share a source endpoint.

    A supplied ordering fixes causal order at a shared port. Independent
    events have no edge. Event ordinals are addresses, not physical times.
    """
    return tuple((i, j) for i in range(len(seams)) for j in range(i+1, len(seams))
                 if set(seams[i]) & set(seams[j]))


def execute_event_tape(level, total, seams, initial, order=None):
    """Original interpreter execution of one linear extension of the event DAG.

    Every event has a distinct blank record slot. The terminal relabelling
    identifies its event address; it does not erase a record or a FIFO state.
    """
    seams = tuple(map(tuple, seams))
    order = tuple(range(len(seams))) if order is None else tuple(order)
    if sorted(order) != list(range(len(seams))):
        raise ValueError('one permutation of all event addresses required')
    positions = {event: position for position, event in enumerate(order)}
    if any(positions[i] >= positions[j] for i, j in causal_edges(seams)):
        raise ValueError('event order violates a shared-port causal dependency')
    program = admitted_program(level, total,
        [fixed_command(seams[i]) for i in order], capacity=max(1, len(seams)))
    result = native.execute(program, initial)
    return program, result


def addressed_blocks(history, order):
    """Complete terminal source/reference blocks indexed by event records."""
    result = {}
    for branch in history.branches:
        records = branch.archive+branch.memory
        if len(records) != len(order):
            raise ValueError('addressed finite tape requires all events completed')
        key = tuple((event, native.canonical(asdict(record)))
                    for event, record in sorted(zip(order, records, strict=True)))
        if key in result:
            raise ArithmeticError('two raw record histories were merged')
        result[key] = branch.density
    return result


def endpoint_certificate():
    data = []
    for law in ('C', 'E', 'mixture'):
        outgoing = [F(0)]*4
        incoming = [F(0)]*4
        for row in raw_local_branches(law):
            assert len(dict(row.pairs)) == len(row.pairs)
            assert len({b for _, b in row.pairs}) == len(row.pairs)
            for a, b in row.pairs:
                outgoing[a] += row.weight
                incoming[b] += row.weight
                assert a.bit_count() == b.bit_count()
                if row.weight:
                    assert (a & 1)-(b & 1) == row.current
                    assert ((b >> 1) & 1)-((a >> 1) & 1) == row.current
            data.append([law, list(row.record), str(row.weight), row.pairs, row.current])
        assert outgoing == incoming == [F(1)]*4
        for left, right in product(range(4), repeat=2):
            actual = unread_column(left, right, (0, 1), law)
            assert sum(v for (a, b), v in actual.items() if a == b) == int(left == right)
    # Disjoint events commute for every four-port matrix unit and every
    # pair of raw slots. This proves the arbitrary-spectator extension.
    count = 0
    for left, right in product(range(16), repeat=2):
        for first, second in product(raw_local_branches(), repeat=2):
            def compose(x, order):
                for seam, branch in order:
                    if x is None:
                        break
                    x = apply_branch(x, seam, branch)
                return x
            ab = (((0, 1), first), ((2, 3), second))
            ba = tuple(reversed(ab))
            assert (compose(left, ab), compose(right, ab)) == (compose(left, ba), compose(right, ba))
            count += 1
    return {'all_endpoint_matrix_units_per_law': 16, 'laws': ['C', 'E', 'mixture'],
        'complete_raw_slots_per_attempt': 25, 'local_maps': data,
        'exact_completeness_unitality_and_load_continuity': True,
        'disjoint_outcome_matrix_unit_comparisons': count,
        'arbitrary_quantum_reference': 'each endpoint Kraus operator tensor identity',
        'finite_theorem_attachment': 'actual endpoint tensor factors and complete Kraus family instantiate slot locality; disjoint-slot commutation gives linear-extension invariance'}


def _encoded_density(density):
    return [[list(key), str(value)] for key, value in sorted(density.items()) if value]


def compact_history(history):
    """Lossless zero-complement representation of the executed raw tape.

    All rows are visited. The declared alphabet and rectangular raw-slot
    indexing retain every zero history; nonzero blocks retain their entire
    source/reference matrix rather than only a selected probability.
    """
    nonzero = []
    digest_payload = []
    for index, branch in enumerate(history.branches):
        row = {'raw_history_index': index,
               'records': [asdict(r) for r in branch.archive+branch.memory],
               'terminal_record': asdict(branch.record),
               'density': _encoded_density(branch.density), 'halted': branch.halted}
        digest_payload.append(row)
        if branch.density:
            nonzero.append(row)
    return {'status': history.status, 'attempts_executed': history.attempted_steps,
            'complete_raw_blocks': len(history.branches),
            'exactly_zero_blocks': sum(not b.density for b in history.branches),
            'nonzero_blocks': nonzero,
            'zero_complement': 'Every unlisted raw_history_index in range(complete_raw_blocks) has exactly the zero source/reference matrix; its records are the literal Cartesian product in interpreter order.',
            'all_raw_histories_sha256': native.digest(digest_payload),
            'total_probability_exact': str(sum(native.trace(b.density) for b in history.branches)),
            'resources': history.resources}


def federation_tape(level):
    g = graph(level)
    a, b = g.inter_seams[0]
    c = min(v for v in g.neighbours[a] if v//12 == a//12)
    d = min(v for v in g.neighbours[b] if v//12 == b//12)
    seams = (tuple(sorted((a, c))), tuple(sorted((b, d))), (a, b))
    assert len({a, b, c, d}) == 4 and causal_edges(seams) == ((0, 2), (1, 2))
    initial = {(1 << x, r, 1 << y, s): F(1, 2)
               for x, r in ((a, 0), (b, 1)) for y, s in ((a, 0), (b, 1))}
    declaration, forward = execute_event_tape(level, 1, seams, initial)
    _, exchanged = execute_event_tape(level, 1, seams, initial, (1, 0, 2))
    assert forward.status == exchanged.status == 'COMPLETE'
    assert addressed_blocks(forward, (0, 1, 2)) == addressed_blocks(exchanged, (1, 0, 2))
    expected_reference = native.reference_marginal([initial])
    for history in (forward, exchanged):
        assert native.reference_marginal(b.density for b in history.branches) == expected_reference
        assert sum(native.trace(b.density) for b in history.branches) == 1
        assert all(i.bit_count() == j.bit_count() == 1 for block in history.branches
                   for i, _, j, _ in block.density)
    return {'level': level, 'cells': g.cells, 'ports': g.ports, 'seams': len(g.seams),
            'source_preparation': _encoded_density(initial), 'event_seams': seams,
            'causal_edges': causal_edges(seams), 'linear_extensions_executed': [[0, 1, 2], [1, 0, 2]],
            'program': declaration, 'complete_terminal_instrument_equal_after_event_relabelling': True,
            'source_total_exact': 1, 'reference_marginal_exact': [[0, 0, '1/2'], [1, 1, '1/2']],
            'output': compact_history(forward),
            'second_linear_extension_all_raw_sha256': compact_history(exchanged)['all_raw_histories_sha256']}


def local_shell_and_stop_controls():
    shells = []
    for shell in ('H66', 'H78'):
        size = int(shell[1:])
        declaration = admitted_local_program(shell, [fixed_command((0, 1))])
        initial = {(i, i, i, i): F(1, size) for i in range(size)}
        history = native.execute(declaration, initial)
        assert native.reference_marginal(b.density for b in history.branches) == native.reference_marginal([initial])
        assert sum(native.trace(b.density) for b in history.branches) == 1
        shells.append({'shell': shell, 'all_basis_columns_tagged_and_checked': size,
                       'complete_raw_blocks': len(history.branches),
                       'total_probability_exact': '1'})
    edge = graph(0).inter_seams[0]
    initial = {(1 << edge[0], 0, 1 << edge[0], 0): F(1)}
    command = fixed_command(edge, repeat=3)
    timeout = native.execute(admitted_program(0, 1, [command], 1, 'halt'), initial)
    stopped = native.execute(admitted_program(0, 1, [command], 3), initial, max_steps=1)
    assert timeout.status == 'FIFO_TIMEOUT' and timeout.attempted_steps == 2
    assert stopped.status == 'RESOURCE_STOP' and stopped.attempted_steps == 1
    assert all(b.halted for b in timeout.branches)
    for history in (timeout, stopped):
        assert sum(native.trace(b.density) for b in history.branches) == 1
        assert native.reference_marginal(b.density for b in history.branches) == {(0, 0): F(1)}
    control = admitted_program(0, 1, [fixed_command(edge)])
    normal = native.execute(control, initial)
    control['parameters']['repair_probability'] = '1/2'
    control['declared_resources'] = native._native_cost(control)
    changed = native.execute(control, initial)
    def arrival(history):
        return sum(v for block in history.branches for (i, r, j, s), v in block.density.items()
                   if i == j == 1 << edge[1] and r == s)
    assert arrival(normal) == F(1, 2) and arrival(changed) == F(1, 4)
    return {'all_local_shell_columns': shells,
        'FIFO_timeout': compact_history(timeout), 'host_resource_stop': compact_history(stopped),
        'altered_rule_control': {'changed_declared_repair_probability': '1/2',
            'same_preparation_seam_and_record_alphabet': True,
            'unchanged_rule_arrival_probability_exact': '1/2',
            'altered_rule_arrival_probability_exact': '1/4',
            'altered_complete_raw_blocks': len(changed.branches)}}


def run():
    return {'verdict': 'PASS_LITERAL_NATIVE_PROCESS_AND_CAUSAL_RECORDS',
        'architecture': {'source': 'unchanged model/native_program/interpreter.py',
            'endpoint_rule': 'unchanged model/observer/native.py:native',
            'port_weights': [1]*12, 'C_E_mixture_exact': ['1/2', '1/2'],
            'repair_probability_exact': '1',
            'scheduler': 'literal finite partial order; uniform-seam selection is a separately declared admitted experiment',
            'memory': 'one fresh record slot per event, complete terminal event-address relabelling',
            'source_state': 'global binary fixed-number sector; H66/H78 are local controls',
            'Z6_link_field': 'absent from this native carrier',
            'physical_time_or_energy_calibration': 'none',
            'target_Hamiltonian_accepted': False},
        'endpoint_instrument': endpoint_certificate(),
        'complete_federation_tapes': [federation_tape(level) for level in (0, 1)],
        'controls': local_shell_and_stop_controls(),
        'scope': 'Complete finite native events and observer records on the declared carriers. Event addresses impose no global physical clock. No physical particle dictionary, wave dispersion, GeV scale or continuum limit is inferred.'}
