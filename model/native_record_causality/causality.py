"""Exact disjoint-command commutation, including every native record slot.

This is an observation-map equivalence between fixed schedules, not a
statement that the actual FIFO order is unobservable or Lorentzian.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product

from ..native_program import interpreter as source
from ..native_causal_dataflow import dataflow


@dataclass(frozen=True)
class RawSlot:
    label: tuple
    pairs: tuple[tuple[int, int], ...]


def _shell(shell):
    if type(shell) is not str or shell not in ('H66', 'H78'):
        raise ValueError('the full native H66 or H78 carrier is required')
    return shell


def _edge(edge):
    if (type(edge) is not tuple or len(edge) != 2
            or any(type(x) is not int for x in edge)
            or edge not in source.local.SEAMS):
        raise ValueError('one canonical native seam is required')
    return edge


def raw_slots(shell, edge):
    """Every raw C/E address and the complete repair-hold slot.

    The partial permutation is separated from its source-controlled scalar
    weight. On a fixed seam that weight is input- and record-independent.
    """
    shell, edge = _shell(shell), _edge(edge)
    result = [RawSlot(('hold',), tuple((i, i) for i in range(int(shell[1:]))))]
    for law in ('C', 'E'):
        for coin in (0, 1):
            values = source.circuit.record_values(law, coin)
            for raw in range(1 << source.circuit.RECORD_BITS[law]):
                mapping = source._local_mapping(shell, law, edge, coin, values.get(raw))
                if len(mapping) != len(set(mapping.values())):
                    raise ArithmeticError('a raw source branch ceased to be a partial injection')
                result.append(RawSlot((law, coin, raw), tuple(sorted(mapping.items()))))
    if len(result) != 25 or len({row.label for row in result}) != 25:
        raise ArithmeticError('the complete native record alphabet changed')
    return tuple(result)


def compose(first, second):
    """Apply first, then second, with empty support kept as a zero operator."""
    right = dict(second)
    return tuple((i, right[j]) for i, j in first if j in right)


def commutation_certificate():
    rows = []
    for shell in ('H66', 'H78'):
        slots = {edge: raw_slots(shell, edge) for edge in source.local.SEAMS}
        pairs = tested = nonzero = 0
        for first, second in combinations(source.local.SEAMS, 2):
            if set(first) & set(second):
                continue
            pairs += 1
            for a, b in product(slots[first], slots[second]):
                forward, backward = compose(a.pairs, b.pairs), compose(b.pairs, a.pairs)
                if forward != backward:
                    raise ArithmeticError('disjoint native raw quantum operators do not commute')
                tested += 1
                nonzero += bool(forward)
        if pairs != 315 or tested != 196875:
            raise ArithmeticError('the exhaustive disjoint-seam census changed')
        rows.append({'carrier': shell, 'dimension': int(shell[1:]),
                     'disjoint_unordered_seam_pairs': pairs,
                     'raw_slots_per_command': 25, 'raw_operator_pairs': tested,
                     'nonzero_composed_operators': nonzero,
                     'zero_composed_operators_retained': tested-nonzero})
    return {'full_carriers': rows, 'partial_operator_commutators_exactly_zero': True,
            'all_input_coherences_and_untouched_references': True,
            'proof': 'endpoint domains and updates depend on disjoint load coordinates'}


def command(edge, law='mixture'):
    _edge(edge)
    if law not in ('C', 'E', 'mixture'):
        raise ValueError('a fixed native C/E law or tagged mixture is required')
    return {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(edge)},
            'law': law, 'repeat': 1}


def paired_programs(*, shell='H66', capacity=1, law='mixture',
                    edges=((0, 5), (1, 7)), repair_probability=F(3, 4), overflow='export_oldest'):
    _shell(shell)
    if (type(edges) is not tuple or len(edges) != 2
            or set(_edge(edges[0])) & set(_edge(edges[1]))):
        raise ValueError('two disjoint fixed native seam commands are required')
    if type(repair_probability) not in (int, F) or not 0 <= repair_probability <= 1:
        raise ValueError('an exact repair probability in [0,1] is required')
    first, second = (command(edge, law) for edge in edges)
    programs = [source.program(shell=shell, commands=commands, capacity=capacity, overflow=overflow)
                for commands in ([first, second], [second, first])]
    for value in programs:
        value['parameters']['repair_probability'] = str(F(repair_probability))
        value['declared_resources'] = source._native_cost(value)
        source.admit(value)
    return tuple(programs)


def split_tape(tape, capacity):
    if type(tape) is not tuple or type(capacity) is not int or capacity < 1:
        raise ValueError('a complete record tuple and positive FIFO capacity are required')
    cut = max(0, len(tape)-capacity)
    return tape[:cut], tape[cut:]


def permute_tape(tape, capacity, order):
    """Reversibly compare named command instances, retaining every raw slot."""
    split_tape(tape, capacity)
    if (type(order) is not tuple or any(type(i) is not int for i in order)
            or sorted(order) != list(range(len(tape)))):
        raise ValueError('a complete record-position permutation is required')
    return split_tape(tuple(tape[i] for i in order), capacity)


def comparison_blocks(history, capacity, *, swap=False):
    """Keep source/reference operators and every archive/FIFO output block."""
    if type(swap) is not bool:
        raise ValueError('the comparison order must be explicit')
    if history.status != 'COMPLETE' or history.attempted_steps != 2:
        raise ValueError('only two completely executed commands admit this comparison')
    if len(history.branches) != 625:
        raise ValueError('all 625 raw outcome blocks must be retained')
    blocks = {}
    for outcome in history.branches:
        tape = outcome.archive+outcome.memory
        if outcome.halted or len(tape) != 2:
            raise ValueError('terminal or missing native records cannot be discarded')
        if split_tape(tape, capacity) != (outcome.archive, outcome.memory):
            raise ArithmeticError('actual FIFO/export partition differs from the complete tape')
        key = permute_tape(tape, capacity, (1, 0) if swap else (0, 1))
        if key in blocks:
            raise ArithmeticError('raw outcome labels were merged')
        blocks[key] = outcome.density
    return blocks


def _index(shell, occupied):
    return source.local.basis(shell == 'H66').index(
        tuple(occupied.count(i) for i in range(12)))


def complete_pair(shell, capacity, *, imaginary=False, law='mixture'):
    """Actual source runs, with coherence tied to an untouched reference."""
    x, y = (_index(shell, pair) for pair in ((0, 1), (5, 7)))
    if imaginary:
        # Imaginary part of an entangled pure density; complex linearity is exact.
        operator = {(x, 0, y, 1): -F(1, 2), (y, 1, x, 0): F(1, 2)}
    else:
        operator = {(x, 0, x, 0): F(1, 2), (x, 0, y, 1): F(1, 2),
                    (y, 1, x, 0): F(1, 2), (y, 1, y, 1): F(1, 2)}
    values = paired_programs(shell=shell, capacity=capacity, law=law)
    left, right = (source.execute(value, operator) for value in values)
    first = comparison_blocks(left, capacity)
    second = comparison_blocks(right, capacity, swap=True)
    if first != second or len(first) != 625:
        raise ArithmeticError('complete source/reference output differs after reversible record comparison')
    if not imaginary:
        if sum(source.trace(x) for x in first.values()) != 1:
            raise ArithmeticError('complete probabilities do not normalize')
        if source.reference_marginal(first.values()) != {(0, 0): F(1, 2), (1, 1): F(1, 2)}:
            raise ArithmeticError('an untouched reference marginal changed')
    zero = sum(not value for value in first.values())
    return {'carrier': shell, 'FIFO_capacity': capacity, 'law': law,
            'input': 'imaginary coherence component' if imaginary else 'entangled pure density',
            'complete_blocks': len(first), 'zero_blocks_retained': zero,
            'archive_slots_per_output': max(0, 2-capacity),
            'same_after_reversible_record_permutation': True,
            'raw_chronological_outputs_equal': first == comparison_blocks(right, capacity),
            'attempts_per_program': 2, 'exported_records_not_erased': True}


def _occupancy_probability(history, shell, port):
    basis = source.local.basis(shell == 'H66')
    return sum((coefficient*int(basis[i][port] > 0)
                for outcome in history.branches
                for (i, a, j, b), coefficient in outcome.density.items()
                if i == j and a == b), F())


def boundaries():
    adaptive = dataflow.certificate()
    if adaptive['receiver_probabilities'] != ['1/4', '0']:
        raise ArithmeticError('the actual adaptive controller boundary changed')
    halted = paired_programs(capacity=1, repair_probability=F(1), overflow='halt')
    x = _index('H66', (0, 1));rho = {(x, 0, x, 0): F(1)}
    runs = [source.execute(value, rho) for value in halted]
    if any(run.status != 'FIFO_TIMEOUT' for run in runs):
        raise ArithmeticError('the intrinsic halt boundary was not retained')
    halt_readings = [_occupancy_probability(run, 'H66', 5) for run in runs]
    if halt_readings != [F(8, 33), F()]:
        raise ArithmeticError('the halted schedules no longer separate')
    overlapping = (command((0, 5)), command((0, 1)))
    x = _index('H66', (0, 2));rho = {(x, 0, x, 0): F(1)}
    runs = [source.execute(source.program(commands=list(commands)), rho)
            for commands in (overlapping, tuple(reversed(overlapping)))]
    overlap_readings = [_occupancy_probability(run, 'H66', 5) for run in runs]
    if overlap_readings != [F(8, 33), F(4, 33)]:
        raise ArithmeticError('the overlapping-seam control no longer separates')
    return {'adaptive_controller': {'receiver_probabilities': adaptive['receiver_probabilities'],
                'raw_histories_per_input': adaptive['complete_raw_branch_counts'],
                'real_controller_signal_cannot_be_quotiented_away': True},
            'intrinsic_FIFO_timeout': {'receiver_probabilities': list(map(str, halt_readings)),
                                      'terminal_outcomes_retained': True},
            'overlapping_fixed_seams': {'receiver_probabilities': list(map(str, overlap_readings))}}


def certificate():
    return {'schema': 'oph.native_record_causality.v1',
            'exact_operator_census': commutation_certificate(),
            'actual_interpreter_comparisons': [complete_pair(shell, capacity)
                for shell, capacity in product(('H66', 'H78'), (1, 2, 3))],
            'complex_linearity_controls': [complete_pair(shell, 1, imaginary=True)
                                          for shell in ('H66', 'H78')],
            'boundaries': boundaries(),
            'scope': {'same_native_endpoint_law_and_complete_record_alphabets': True,
                'comparison_is_reversible_record_permutation': True,
                'all_input_reference_channel_identity': True,
                'induced_diamond_distance_after_record_comparison': '0',
                'chronological_records_are_identified_without_readout_map': False,
                'adaptive_feedback_covered': False, 'early_halting_covered': False,
                'raw_single_FIFO_dependency_order_changed': False,
                'source_selected_physical_causal_quotient': False,
                'three_spatial_dimensions_or_Lorentz_limit_derived': False,
                'interacting_field_or_particle_mass_derived': False}}
