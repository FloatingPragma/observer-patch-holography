"""Native incident-clock fields in the conserved, fixed-number equilibrium.

Independent Bernoulli preparation is unnecessary. Each full screen supplies
its own conserved total and exact canonical restriction to the retained ball.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from itertools import combinations, product
import math
from pathlib import Path
import random
from typing import NamedTuple

from tools import strict_json
from . import core, dual, refinement as ref
from ..native_action import clock_histories as runtime, incident_refinement as clocks
from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*ref.SOURCE_PATHS,
    'model/observer_time_fields/canonical_sector.py'))))


def _sector(population, occupied):
    if (type(population) is not int or not 2 <= population <= 1000000
            or type(occupied) is not int or not 0 <= occupied <= population):
        raise ValueError('exact finite population and conserved occupation required')


def falling(number, degree):
    if (type(number) is not int or number < 0 or type(degree) is not int
            or not 0 <= degree <= 4):
        raise ValueError('nonnegative population and degree zero through four required')
    return math.prod(range(number-degree+1, number+1)) if degree <= number else 0


def initial_weight(population, occupied, initial):
    _sector(population, occupied)
    if (type(initial) is not tuple or len(initial) != 2
            or any(type(bit) is not int or bit not in (0, 1) for bit in initial)):
        raise ValueError('two original disclosed binary bits required')
    ones = sum(initial)
    return F(falling(occupied, ones)*falling(population-occupied, 2-ones),
             falling(population, 2))


def _expectation(summaries, initial, population, occupied, grid, mass):
    value, coefficients = F(0), 0
    remaining = occupied-sum(initial)
    for summary in summaries:
        for (mask, degree), count in summary:
            coefficients += count
            if (not mask & 1 or initial[0]) and (not mask & 2 or initial[1]):
                value += F(count, grid)*F(falling(remaining, degree), falling(population-2, degree))
    deficit = mass-F(coefficients, grid)
    if not 0 <= deficit <= mass or not 0 <= value <= mass:
        raise ArithmeticError('canonical expectation lost its exact enclosure')
    return value, min(mass, value+deficit)


class FieldLaw(NamedTuple):
    ports: tuple
    edges: tuple
    incident: tuple
    records: tuple
    lower_mass: tuple
    upper_mass: tuple
    timeout: F
    ticks: int
    cap: int
    bits: int
    population: int
    occupied: int
    feedback: bool

    def intervals(self):
        value = validate(self)
        return dict(zip(value.records, zip(value.lower_mass, value.upper_mass), strict=True))


def law(ports, edges, incident, *, population, occupied, ticks=4, cap=512, bits=52, feedback=False):
    """All success/timeout records for the exact full-sector marginal on this graph."""
    _sector(population, occupied)
    ref._integer(ticks, 2, 8); ref._integer(cap, 1, 2048); ref._integer(bits, 8, 52)
    dual.build(ports, edges, incident)
    if (type(feedback) is not bool or population < len(ports)
            or not {0, 1} <= set(ports) or (0, 1) not in incident
            or any(set(edge).intersection((0, 1)) and edge not in incident for edge in edges)):
        raise ValueError('declared full sector and complete read-port incident access required')
    first, final, ordinary, controlled = ref._coefficients(ports, edges, incident, ticks, cap, bits)
    q = 1-(1-F(len(incident), len(edges)))**cap
    completed = q**ticks
    records, lower, upper = [], [], []

    def retain(record, interval, weight):
        low, high = max(F(0), interval[0]), min(F(1), interval[1])
        if low > high:
            raise ArithmeticError('canonical field interval is empty')
        records.append(record); lower.append(weight*low); upper.append(weight*high)

    for initial in product((0, 1), repeat=2):
        weight = initial_weight(population, occupied, initial)
        def expectation(summary, mass):
            return _expectation(summary, initial, population, occupied, 1 << bits, mass) if weight else (F(0), F(0))
        f = expectation(first, q)
        f_completed = ref._scale(f, q**(ticks-1))
        last = expectation(final, completed)
        joint = expectation(ordinary, completed)
        active = expectation(controlled if feedback else ordinary, completed)
        success = {(1, 1): active, (1, 0): ref._subtract(f_completed, active),
            (0, 1): ref._subtract(last, joint),
            (0, 0): ref._add(ref._subtract(ref._subtract((completed, completed), f_completed), last), joint)}
        for a, b in product((0, 1), repeat=2):
            retain((initial, a, b, ticks, 'complete'), success[a, b], weight)
        retain((initial, None, None, 0, 'timeout'), (1-q, 1-q), weight)
        for done in range(1, ticks):
            for bit, probability in ((0, ref._subtract((q, q), f)), (1, f)):
                retain((initial, bit, None, done, 'timeout'),
                       ref._scale(probability, q**(done-1)*(1-q)), weight)
    result = FieldLaw(ports, edges, incident, tuple(records), tuple(lower), tuple(upper),
                      1-completed, ticks, cap, bits, population, occupied, feedback)
    if not 0 <= 1-sum(lower) <= 1 or sum(upper) < 1:
        raise ArithmeticError('canonical complete probability budget was lost')
    return result


def validate(value):
    if (type(value) is not FieldLaw or type(value.records) is not tuple
            or type(value.lower_mass) is not tuple or type(value.upper_mass) is not tuple
            or any(type(p) is not F for p in (*value.lower_mass, *value.upper_mass))
            or type(value.timeout) is not F):
        raise ValueError('an original canonical field law required')
    for record in value.records:
        if (type(record) is not tuple or len(record) != 5
                or type(record[0]) is not tuple or len(record[0]) != 2
                or any(type(bit) is not int or bit not in (0, 1) for bit in record[0])
                or any(bit is not None and (type(bit) is not int or bit not in (0, 1)) for bit in record[1:3])
                or type(record[3]) is not int or type(record[4]) is not str):
            raise ValueError('original canonical record types required')
    expected = law(value.ports, value.edges, value.incident, population=value.population,
        occupied=value.occupied, ticks=value.ticks, cap=value.cap, bits=value.bits, feedback=value.feedback)
    if value != expected:
        raise ValueError('canonical field law differs from its source calculation')
    return expected


def statistics(value, *, spatial=F(0)):
    if type(spatial) is not F or not 0 <= spatial <= 1:
        raise ValueError('an exact spatial allowance required')
    value = validate(value)
    moments = {subset: ref._moment(value, subset, spatial)
               for degree in range(1, 5) for subset in combinations(range(4), degree)}
    cumulant = F(0), F(0)
    for partition in core._partitions(list(range(4))):
        interval = F(1), F(1)
        for block in partition:
            other = moments[tuple(sorted(block))]
            interval = interval[0]*other[0], interval[1]*other[1]
        term = ref._scale(interval, F(math.factorial(len(partition)-1)))
        cumulant = ref._subtract(cumulant, term) if len(partition) % 2 == 0 else ref._add(cumulant, term)
    return {'means': [moments[i,] for i in range(4)], 'joint_four': moments[tuple(range(4))],
            'connected_four': cumulant}


def export(value, spatial):
    value = validate(value)
    lower = tuple(ref.dyadic(p) for p in value.lower_mass)
    stats = statistics(value, spatial=spatial)
    return {'feedback': value.feedback,
        'records': [{'initial_bits': list(record[0]), 'first_field': record[1],
            'final_field': record[2], 'scheduled_ticks': record[3], 'status': record[4],
            'probability_interval': [str(lo), str(ref.dyadic(hi, upper=True))]}
            for record, lo, hi in zip(value.records, lower, value.upper_mass, strict=True)],
        'capped_complete_readout_l1_allowance': str(ref.dyadic(1-sum(lower), upper=True)),
        'eventual_readout_missing_mass_allowance': str(ref.dyadic(1-sum(
            (p for r, p in zip(value.records, lower, strict=True) if r[-1] == 'complete'), F(0)), upper=True)),
        'interval_cap_timeout_allowance': str(ref.dyadic(value.timeout, upper=True)),
        'full_screen_eventual_means': [ref._interval(interval) for interval in stats['means']],
        'full_screen_eventual_joint_four': ref._interval(stats['joint_four']),
        'full_screen_eventual_connected_four': ref._interval(stats['connected_four'])}


def screen_distance(first, second, spatial):
    """Two normalized completed field laws, allowing distinct initial marginals."""
    first, second = validate(first), validate(second)
    if type(spatial) is not F or not 0 <= spatial <= 1:
        raise ValueError('an exact one-screen locality allowance required')
    attributes = ('ports', 'edges', 'incident', 'ticks', 'cap', 'bits', 'feedback')
    if any(getattr(first, name) != getattr(second, name) for name in attributes):
        raise ValueError('the same retained graph and observation protocol required')
    def successful(value):
        return {record: (lo, min(F(1), hi+value.timeout))
            for record, lo, hi in zip(value.records, value.lower_mass, value.upper_mass, strict=True)
            if record[-1] == 'complete'}
    a, b = successful(first), successful(second)
    lower = sum((max(F(0), a[key][0]-b[key][1], b[key][0]-a[key][1]) for key in a), F(0))/2
    upper = sum((max(abs(a[key][0]-b[key][1]), abs(a[key][1]-b[key][0])) for key in a), F(0))/2
    return max(F(0), lower-2*spatial), min(F(1), upper+2*spatial)


def equilibrium_certificate():
    """Check the local source and connected graphs needed by sector ergodicity."""
    states = tuple(product((0, 1), repeat=2))
    counts = [[0]*4 for _ in states]
    for i, state in enumerate(states):
        outputs = [native.native(state, (0, 1), coin)[0] for coin in (0, 1)]
        if sorted(outputs) != sorted((state, tuple(reversed(state)))):
            raise ArithmeticError('the native diagonal repair is not identity plus swap')
        for output in outputs:
            counts[states.index(output)][i] += 1
    if counts != [list(row) for row in zip(*counts)] or any(sum(row) != 2 for row in counts):
        raise ArithmeticError('native local reversibility was lost')
    graphs = []
    for level in runtime.LEVELS:
        graph = clocks.screen.screen(level)
        adjacency = graph.graph.neighbour_lists
        reached, queue = {0}, [0]
        for port in queue:
            for neighbour in adjacency[port]:
                neighbour = int(neighbour)
                if neighbour not in reached:
                    reached.add(neighbour); queue.append(neighbour)
        if len(reached) != graph.ports:
            raise ArithmeticError('a source screen is disconnected')
        graphs.append({'level': level, 'ports': int(graph.ports), 'seams': int(graph.seams),
            'reachable_ports': len(reached), 'conserved_occupied_ports': int(graph.ports)//3})
    return {'native_two_port_integer_counts': counts, 'screens': graphs,
        'unique_diagonal_stationary_law_within_each_number_sector': 'uniform on configurations with the declared conserved total',
        'sector_value_selected_by_repair': False,
        'finite_preparation_mixing_time_certified': False}


def draw_sector(rng, population, occupied):
    """Uniform subset under exact uniform integer draws; no floating probabilities.

    Seeded PRNG executions are reproducible examples of this sampler, not
    proof that a finite seed ensemble represents the full equilibrium law.
    """
    _sector(population, occupied)
    left = occupied
    state = bytearray(population)
    for port in range(population):
        state[port] = int(rng.randrange(population-port) < left)
        left -= state[port]
    if left or sum(state) != occupied:
        raise ArithmeticError('fixed-number sampler changed the conserved total')
    return bytes(state)


def _counts(counts):
    if type(counts) is not tuple or len(counts) != 2:
        raise ValueError('one original conserved count per full screen required')
    for population, occupied in zip(runtime.source().ports, counts, strict=True):
        _sector(population, occupied)


@lru_cache(maxsize=8)
def _initial(seed, ticks, cap, feedback, counts):
    data = runtime.source()
    rng = random.Random(seed)
    states = tuple(draw_sector(rng, population, occupied)
        for population, occupied in zip(data.ports, counts, strict=True))
    addresses = tuple(map(dict, data.matched_addresses))
    fields = tuple(state[table[0]] for state, table in zip(states, addresses, strict=True))
    return runtime.Checkpoint(data.source_digest, seed, ticks, cap, feedback, 0, 0,
        (0, 0), (0, 0), fields, fields, counts,
        tuple(hashlib.sha256(state).hexdigest() for state in states), states, rng.getstate(), ((), ()))


@lru_cache(maxsize=8)
def _compute(seed, ticks, cap, feedback, stop, counts):
    return runtime._continue(_initial(seed, ticks, cap, feedback, counts), stop)


class ObserverCheckpoint(NamedTuple):
    source: runtime.Checkpoint
    initial_bits: tuple
    counters: tuple
    first_fields: tuple


def _wrap(source):
    initial = _initial(source.seed, source.ticks, source.cap, source.feedback, source.initial_totals)
    addresses = tuple(map(dict, runtime.source().matched_addresses))
    disclosed = tuple(tuple(state[table[p]] for p in (0, 1))
        for state, table in zip(initial.states, addresses, strict=True))
    fields = tuple(ref._counter_view(bits, records)
        for bits, records in zip(disclosed, source.records, strict=True))
    return ObserverCheckpoint(source, disclosed, tuple(row[0] for row in fields), tuple(row[1] for row in fields))


def validate_checkpoint(value):
    if (type(value) is not ObserverCheckpoint or type(value.source) is not runtime.Checkpoint
            or not all(runtime._builtins(field) for field in value.source)
            or not all(runtime._builtins(field) for field in value[1:])):
        raise ValueError('original canonical source and observer checkpoint required')
    source = value.source
    runtime._parameters(source.seed, source.ticks, source.cap, source.feedback, source.slots)
    _counts(source.initial_totals)
    expected = _wrap(_compute(source.seed, source.ticks, source.cap, source.feedback,
                             source.slots, source.initial_totals))
    if (not all(runtime._exact_equal(a, b) for a, b in zip(source, expected.source, strict=True))
            or not all(runtime._exact_equal(a, b) for a, b in zip(value[1:], expected[1:], strict=True))):
        raise ValueError('canonical checkpoint differs from its actual source prefix')
    return expected


def execute(seed, *, counts=(5120, 81920), ticks=4, cap=runtime.MAX_CAP,
            feedback=False, stop_after_slots=None, checkpoint=None):
    stop = cap if stop_after_slots is None else stop_after_slots
    runtime._parameters(seed, ticks, cap, feedback, stop); _counts(counts)
    if checkpoint is None:
        return _wrap(_compute(seed, ticks, cap, feedback, stop, counts))
    value = validate_checkpoint(checkpoint)
    old = value.source
    if (old.seed, old.ticks, old.cap, old.feedback, old.initial_totals) != (seed, ticks, cap, feedback, counts):
        raise ValueError('continuation cannot replace its conserved sector or controller')
    if stop < old.slots:
        raise ValueError('continuation cannot rewind its source')
    continued = runtime._continue(old, stop)
    fields = tuple(ref._counter_view(value.counters[side],
        continued.records[side][len(old.records[side]):], value.first_fields[side]) for side in (0, 1))
    return ObserverCheckpoint(continued, value.initial_bits, tuple(row[0] for row in fields),
                              tuple(row[1] for row in fields))


def _payload(value):
    # Compare the actual JSON representation, retaining exact int/Boolean
    # distinctions while converting immutable source tuples to JSON arrays.
    source = json.loads(json.dumps(runtime._payload(value.source), allow_nan=False))
    return {'schema': 'oph.canonical_field_checkpoint.v1', 'source': source,
        'initial_bits': [list(row) for row in value.initial_bits],
        'counters': [list(row) for row in value.counters], 'first_fields': list(value.first_fields)}


def checkpoint_to_json(value):
    return json.dumps(_payload(validate_checkpoint(value)), sort_keys=True, separators=(',', ':'), allow_nan=False)


def checkpoint_from_json(text):
    if type(text) is not str or len(text) > 1100000:
        raise ValueError('a bounded canonical checkpoint string required')
    data = strict_json.loads(text)
    if type(data) is not dict or type(data.get('source')) is not dict:
        raise ValueError('complete canonical source checkpoint required')
    source = data['source']
    try:
        counts = tuple(source['initial_totals'])
        expected = execute(source['seed'], counts=counts, ticks=source['ticks'], cap=source['cap'],
            feedback=source['feedback'], stop_after_slots=source['slots'])
    except (KeyError, TypeError) as error:
        raise ValueError('canonical checkpoint parameters missing') from error
    compare_receipt(data, _payload(expected))
    return expected


def history(seed, *, feedback=False, cap=runtime.MAX_CAP):
    full = execute(seed, feedback=feedback, cap=cap)
    prefix = execute(seed, feedback=feedback, cap=cap, stop_after_slots=full.source.slots//2)
    encoded = checkpoint_to_json(prefix)
    resumed = execute(seed, feedback=feedback, cap=cap, checkpoint=checkpoint_from_json(encoded))
    if resumed != full:
        raise ArithmeticError('canonical source or observer memory changed on continuation')
    rows = []
    for side, level in enumerate(runtime.LEVELS):
        source = full.source
        addresses = dict(runtime.source().matched_addresses[side])
        actual = tuple(source.states[side][addresses[p]] for p in (0, 1))
        if actual != full.counters[side] or sum(source.states[side]) != source.initial_totals[side]:
            raise ArithmeticError('canonical current readback or total conservation failed')
        complete = source.scheduled_hits == source.ticks
        rows.append({'level': level, 'initial_total': source.initial_totals[side],
            'terminal_total': sum(source.states[side]), 'initial_bits': list(full.initial_bits[side]),
            'first_field': full.first_fields[side], 'terminal_counters': list(full.counters[side]),
            'final_field': full.counters[side][1] if complete else None,
            'status': 'complete' if complete else 'timeout',
            'native_attempts': source.uniform_attempts[side]+source.feedback_attempts[side],
            'scheduled_ticks': source.scheduled_hits, 'feedback_attempts': source.feedback_attempts[side],
            'incident_records': [list(record) for record in source.records[side]],
            'initial_state_sha256': source.initial_state_sha256[side],
            'terminal_state_sha256': hashlib.sha256(source.states[side]).hexdigest()})
    return {'seed': seed, 'feedback': feedback, 'external_union_slot_cap': cap,
        'external_union_slots_used': full.source.slots, 'screens': rows,
        'full_state_record_rng_and_observer_memory_continuation_verified': True,
        'seeded_examples_are_statistical_equilibrium_frequency_evidence': False,
        'independent_screen_preparations_assumed_to_have_identical_restrictions': False}


def run():
    certificate = equilibrium_certificate()
    ports, edges, incident, observed, exits, _ = clocks.geometry()
    spatial, _ = clocks.exit_bound(clocks.exit_counts(ports, edges, observed, exits, 48), 4, len(incident))
    rows, source_laws = [], []
    for info in certificate['screens']:
        population, occupied = info['ports'], info['conserved_occupied_ports']
        laws = [law(ports, edges, incident, population=population, occupied=occupied,
                    feedback=feedback) for feedback in (False, True)]
        source_laws.append(laws)
        stats = [statistics(value, spatial=spatial) for value in laws]
        response = ref._subtract(stats[1]['means'][3], stats[0]['means'][3])
        if not (stats[0]['connected_four'][0] > 0 and stats[1]['connected_four'][1] < 0 and response[0] > 0):
            raise ArithmeticError('canonical field/feedback witnesses are unresolved')
        rho = F(occupied, population)
        rows.append({**info, 'density': str(rho),
            'initial_distinct_port_covariance': str(-rho*(1-rho)/(population-1)),
            'full_screen_initial_restriction_used_exactly': True,
            'full_screen_to_its_own_ball_complete_output_allowance': str(ref.dyadic(spatial, upper=True)),
            'field_protocols': [export(value, spatial) for value in laws],
            'feedback_increase_of_final_field': ref._interval(response)})
    return {'schema': 'oph.observer_time_fields.canonical_sector.v1',
        'equilibrium_certificate': certificate, 'screens': rows,
        'two_screen_completed_field_record_TV_intervals': [
            {'feedback': bool(index), 'interval': ref._interval(screen_distance(
                source_laws[0][index], source_laws[1][index], spatial))} for index in (0, 1)],
        'full_source_histories': [history(seed, feedback=feedback)
            for feedback in (False, True) for seed in range(4)],
        'forced_external_timeout': history(1, feedback=True, cap=1),
        'full_source_external_cutoff_allowance': str(ref.dyadic(clocks.slot_timeout_bound(
            runtime.MAX_CAP, 4, 33, runtime.source().slots), upper=True)),
        'scope': {'independent_Bernoulli_product_preparation_required': False,
            'conserved_sector_density_selected_by_repair': False,
            'stationarity_after_finite_uncontrolled_run_asserted': False,
            'uniform_equilibrium_sampling_service_declared': True,
            'diagonal_record_law_not_general_coherent_state_uniqueness': True,
            'different_finite_sector_restrictions_silently_identified': False,
            'physical_lab_time_or_relativistic_field_equations_derived': False},
        'source_bindings': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path); mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Canonical-sector native clock and field laws reproduced.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
