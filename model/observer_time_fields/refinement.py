"""Complete finite-density field readouts in the actual incident source clock.

Squarefree occupation polynomials of degree at most two suffice for this
four-bit protocol. They are an exact backward dual, not a mean-field closure.
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
from typing import NamedTuple

from . import core, dual
from ..native_action import clock_histories, incident_refinement as clocks
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*core.SOURCE_PATHS, *dual.SOURCE_PATHS,
    *clock_histories.SOURCE_PATHS, *clocks.stopped_control.SOURCE_PATHS,
    'model/native_action/incident_refinement.py',
    'model/observer_time_fields/refinement.py'))))


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
    rho: F
    feedback: bool

    @property
    def mass_deficit(self):
        checked = validate(self)
        return 1-sum(checked.lower_mass, F(0))

    def intervals(self):
        checked = validate(self)
        return dict(zip(checked.records, zip(checked.lower_mass, checked.upper_mass), strict=True))


def _integer(value, lo, hi):
    if type(value) is not int or not lo <= value <= hi:
        raise ValueError('a bounded original integer is required')


def _summary(operator, vector):
    # Group the exact coefficients before evaluating the four preparations.
    # Mask records which initially disclosed sites occur in the monomial.
    result = {}
    for labels, coefficient in zip(dual.states(operator), vector, strict=True):
        mask = int(0 in labels)+2*int(1 in labels)
        remote = len(labels)-mask.bit_count()
        result[mask, remote] = result.get((mask, remote), 0)+coefficient
    return tuple(sorted(result.items()))


def _expectation(summaries, initial, rho, grid, mass):
    value, coefficient_mass = F(0), 0
    for summary in summaries:
        for (mask, degree), count in summary:
            coefficient_mass += count
            if (not mask & 1 or initial[0]) and (not mask & 2 or initial[1]):
                value += count*rho**degree/grid
    deficit = mass-F(coefficient_mass, grid)
    if not 0 <= deficit <= mass or not 0 <= value <= mass:
        raise ArithmeticError('backward source coefficients lost their lower-bound contract')
    return value, min(mass, value+deficit)


def _add(left, right):
    return left[0]+right[0], left[1]+right[1]


def _subtract(left, right):
    return left[0]-right[1], left[1]-right[0]


def _scale(value, factor):
    assert factor >= 0
    return value[0]*factor, value[1]*factor


@lru_cache(maxsize=4)
def _coefficients(ports, edges, incident, ticks, cap, bits):
    one, two = (dual.build(ports, edges, incident, degree=d) for d in (1, 2))
    single = dual.states(one)
    index = {labels[0]: i for i, labels in enumerate(single)}
    pairs = {labels: i for i, labels in enumerate(dual.states(two))}
    grid = 1 << bits

    def basis(port):
        return tuple(grid*int(label == (port,)) for label in single)

    def evolve(op, value, count=1):
        return dual.evolve(op, value, count, cap=cap, bits=bits)

    first = evolve(one, basis(0))
    future = evolve(one, basis(1), ticks-1)
    final = evolve(one, future)

    def joint(vector):
        diagonal, cross = [0]*len(single), [0]*len(pairs)
        diagonal[index[0]] = vector[index[0]]
        for port, value in zip(ports, vector, strict=True):
            if port != 0:
                cross[pairs[tuple(sorted((0, port)))]] = value
        return (_summary(one, evolve(one, tuple(diagonal))),
                _summary(two, evolve(two, tuple(cross))))

    ordinary = joint(future)
    controlled = list(future)
    # Native internal repair on the backward linear observable. A discarded
    # odd grid unit is charged by the final coefficient-mass deficit.
    controlled[index[0]] = controlled[index[1]] = (
        future[index[0]]+future[index[1]])//2
    return ((_summary(one, first),), (_summary(one, final),),
            ordinary, joint(tuple(controlled)))


def law(ports, edges, incident, *, ticks=4, cap=512, bits=52, rho=F(1, 3), feedback=False):
    """Every success/timeout record, with exact rational enclosing intervals.

    Records are (initial_two_bits, first_bit, last_bit, completed_ticks, status).
    Timeout before the first tick has neither bit; later timeout keeps the
    already observed first bit. Native feedback acts after that first read.
    """
    _integer(ticks, 2, 8); _integer(cap, 1, 2048); _integer(bits, 8, 52)
    if type(rho) is not F or not 0 < rho < 1 or type(feedback) is not bool:
        raise ValueError('an exact interior density and Boolean feedback policy are required')
    # Fraction exposes writable protected slots. Do not let returned laws
    # alias the default preparation or a caller-owned density object.
    rho = F(rho.numerator, rho.denominator)
    # Validate original geometry before the coefficient cache can equate keys.
    dual.build(ports, edges, incident, degree=1)
    if (not {0, 1} <= set(ports) or (0, 1) not in incident
            or any(set(edge).intersection((0, 1)) and edge not in incident for edge in edges)):
        raise ValueError('both read ports and their incident native feedback seam are required')
    first, final, ordinary, controlled = _coefficients(ports, edges, incident, ticks, cap, bits)
    grid = 1 << bits
    q = 1-(1-F(len(incident), len(edges)))**cap
    completed = q**ticks
    records, lower, upper = [], [], []

    def retain(record, interval, weight):
        lo, hi = max(F(0), interval[0]), min(F(1), interval[1])
        if lo > hi:
            raise ArithmeticError('inconsistent source readout interval')
        records.append(record); lower.append(weight*lo); upper.append(weight*hi)

    for initial in product((0, 1), repeat=2):
        weight = rho**sum(initial)*(1-rho)**(2-sum(initial))
        f = _expectation(first, initial, rho, grid, q)
        f_completed = _scale(f, q**(ticks-1))
        l0 = _expectation(final, initial, rho, grid, completed)
        j0 = _expectation(ordinary, initial, rho, grid, completed)
        j1 = _expectation(controlled if feedback else ordinary, initial, rho, grid, completed)
        success = {(1, 1): j1, (1, 0): _subtract(f_completed, j1),
                   (0, 1): _subtract(l0, j0),
                   (0, 0): _add(_subtract(_subtract((completed, completed), f_completed), l0), j0)}
        for a, b in product((0, 1), repeat=2):
            retain((initial, a, b, ticks, 'complete'), success[a, b], weight)
        retain((initial, None, None, 0, 'timeout'), (1-q, 1-q), weight)
        for done in range(1, ticks):
            factor = q**(done-1)*(1-q)
            for bit, probability in ((0, _subtract((q, q), f)), (1, f)):
                retain((initial, bit, None, done, 'timeout'), _scale(probability, factor), weight)
    result = FieldLaw(ports, edges, incident, tuple(records), tuple(lower), tuple(upper), 1-completed,
                      ticks, cap, bits, rho, feedback)
    if not 0 <= 1-sum(result.lower_mass) <= 1 or sum(result.upper_mass) < 1:
        raise ArithmeticError('the complete capped field law lost its enclosing mass budget')
    return result


def validate(value):
    """Compare immutable original fields with their source calculation, not a seal."""
    if type(value) is not FieldLaw:
        raise ValueError('an immutable source field law is required')
    if (type(value.records) is not tuple
            or type(value.lower_mass) is not tuple or type(value.upper_mass) is not tuple
            or any(type(p) is not F for p in (*value.lower_mass, *value.upper_mass))
            or type(value.timeout) is not F):
        raise ValueError('original exact field-law types are required')
    for record in value.records:
        if (type(record) is not tuple or len(record) != 5
                or type(record[0]) is not tuple or len(record[0]) != 2
                or any(type(bit) is not int or bit not in (0, 1) for bit in record[0])
                or any(bit is not None and (type(bit) is not int or bit not in (0, 1))
                       for bit in record[1:3])
                or type(record[3]) is not int or type(record[4]) is not str):
            raise ValueError('original exact record types are required')
    expected = law(value.ports, value.edges, value.incident, ticks=value.ticks,
                   cap=value.cap, bits=value.bits, rho=value.rho, feedback=value.feedback)
    if value != expected:
        raise ValueError('field law differs from the source calculation')
    return expected


def dyadic(value, *, upper=False, bits=80):
    """Outward compact export, without a decimal/float conversion."""
    grid = 1 << bits
    units = -((-value.numerator*grid)//value.denominator) if upper else (
        value.numerator*grid)//value.denominator
    return F(units, grid)


def _moment(law, subset, spatial=F(0)):
    # The complete eventual law dominates every retained successful lower
    # mass. All unresolved timeouts and arithmetic losses are charged once.
    success = [(record, lo, hi) for record, lo, hi in zip(
        law.records, law.lower_mass, law.upper_mass, strict=True) if record[-1] == 'complete']
    deficit = 1-sum((lo for _, lo, _ in success), F(0))
    low = sum((lo for (initial, first, last, _, _), lo, _ in success
               if all((*initial, first, last)[i] for i in subset)), F(0))
    high = min(F(1), low+deficit)
    return max(F(0), low-spatial), min(F(1), high+spatial)


def statistics(law, *, spatial=F(0)):
    """Moments/cumulant of the eventual four-bit law, then spatial error."""
    if type(spatial) is not F or not 0 <= spatial <= 1:
        raise ValueError('a source field law and rational spatial allowance are required')
    law = validate(law)
    moments = {subset: _moment(law, subset, spatial)
               for degree in range(1, 5) for subset in combinations(range(4), degree)}
    cumulant = F(0), F(0)
    for partition in core._partitions(list(range(4))):
        interval = F(1), F(1)
        for block in partition:
            other = moments[tuple(sorted(block))]
            interval = interval[0]*other[0], interval[1]*other[1]
        term = _scale(interval, F(math.factorial(len(partition)-1)))
        cumulant = _subtract(cumulant, term) if len(partition) % 2 == 0 else _add(cumulant, term)
    return {'means': [moments[i,] for i in range(4)], 'joint_four': moments[tuple(range(4))],
            'connected_four': cumulant,
            'completed_readout_missing_mass': 1-sum((p for r, p in zip(
                law.records, law.lower_mass, strict=True) if r[-1] == 'complete'), F(0))}


class ObserverCheckpoint(NamedTuple):
    source: clock_histories.Checkpoint
    initial_bits: tuple
    counters: tuple
    first_fields: tuple


def _counter_view(counters, records, first=None):
    fields = list(counters)
    for kind, _, scheduled, _, a, b, _, current, field0 in records:
        for port in (0, 1):
            fields[port] += current*(int(b == port)-int(a == port))
        if fields[0] != field0 or any(value not in (0, 1) for value in fields):
            raise ArithmeticError('authentic currents did not reconstruct both fields')
        if kind == 'scheduled' and scheduled == 1:
            first = fields[0]
    return tuple(fields), first


def _observer_checkpoint(source):
    source = clock_histories.validate_checkpoint(source)
    initial = clock_histories.execute(source.seed, ticks=source.ticks, cap=source.cap,
        feedback=source.feedback, stop_after_slots=0)
    addresses = tuple(map(dict, clock_histories.source().matched_addresses))
    disclosed = tuple(tuple(state[table[p]] for p in (0, 1)) for state, table in zip(
        initial.states, addresses, strict=True))
    projected = tuple(_counter_view(bits, records) for bits, records in zip(
        disclosed, source.records, strict=True))
    return ObserverCheckpoint(source, disclosed, tuple(row[0] for row in projected),
                              tuple(row[1] for row in projected))


def _checkpoint_payload(checkpoint):
    return {'schema': 'oph.incident_field_checkpoint.v1',
        'source_checkpoint': clock_histories.checkpoint_to_json(checkpoint.source),
        'initial_bits': [list(bits) for bits in checkpoint.initial_bits],
        'counters': [list(bits) for bits in checkpoint.counters],
        'first_fields': list(checkpoint.first_fields)}


def checkpoint_to_json(source_checkpoint):
    """Extend the canonical source checkpoint with the two-field observer memory."""
    return json.dumps(_checkpoint_payload(_observer_checkpoint(source_checkpoint)),
                      sort_keys=True, separators=(',', ':'), allow_nan=False)


def checkpoint_from_json(text):
    """Authenticate source ancestry and original observer-memory types before use.

    Source replay is a checkpoint audit, not an extra controller observation.
    Operational continuation below consumes only the saved counters and suffix.
    """
    if type(text) is not str or len(text) > 1100000:
        raise ValueError('a bounded serialized observer checkpoint is required')
    data = strict_json.loads(text)
    if type(data) is not dict or 'source_checkpoint' not in data:
        raise ValueError('the complete observer checkpoint is required')
    source = clock_histories.checkpoint_from_json(data['source_checkpoint'])
    expected = _observer_checkpoint(source)
    compare_receipt(data, _checkpoint_payload(expected))
    return expected


def history(seed, *, feedback=False, cap=clock_histories.MAX_CAP):
    """Project an actual complete-screen execution; continue its real checkpoint."""
    initial = clock_histories.execute(seed, ticks=4, cap=cap, feedback=feedback, stop_after_slots=0)
    full = clock_histories.execute(seed, ticks=4, cap=cap, feedback=feedback)
    prefix = clock_histories.execute(seed, ticks=4, cap=cap, feedback=feedback,
                                    stop_after_slots=full.slots//2)
    encoded = checkpoint_to_json(prefix)
    restored = checkpoint_from_json(encoded)
    continued = clock_histories.execute(seed, ticks=4, cap=cap, feedback=feedback,
        checkpoint=restored.source)
    if continued != full:
        raise ArithmeticError('field-history continuation changed the underlying source')
    source = clock_histories.source()
    rows = []
    for side, level in enumerate(clock_histories.LEVELS):
        addresses = dict(source.matched_addresses[side])
        disclosed = tuple(initial.states[side][addresses[p]] for p in (0, 1))
        suffix = continued.records[side][len(restored.source.records[side]):]
        fields, first = _counter_view(restored.counters[side], suffix, restored.first_fields[side])
        if (restored.initial_bits[side] != disclosed
                or (fields, first) != _counter_view(disclosed, full.records[side])):
            raise ArithmeticError('observer-memory continuation changed the field readout')
        actual = tuple(full.states[side][addresses[p]] for p in (0, 1))
        if tuple(fields) != actual:
            raise ArithmeticError('terminal source occupations disagree with their current records')
        complete = full.scheduled_hits == 4
        rows.append({'level': level, 'initial_bits': list(disclosed), 'first_field': first,
            'final_field': fields[1] if complete else None,
            'status': 'complete' if complete else 'timeout',
            'scheduled_ticks': full.scheduled_hits,
            'native_feedback_attempts': full.feedback_attempts[side],
            'feedback_signed_currents': [row[7] for row in full.records[side] if row[0] == 'feedback'],
            'total_incident_ticks': len(full.records[side]),
            'native_attempts': full.uniform_attempts[side]+full.feedback_attempts[side],
            'terminal_state_sha256': hashlib.sha256(full.states[side]).hexdigest()})
    return {'seed': seed, 'feedback': feedback, 'external_union_slot_cap': cap,
        'external_union_slots_used': full.slots, 'screens': rows,
        'both_fields_reconstructed_from_initial_bits_and_signed_currents': True,
        'serialized_continuation_preserves_full_states_records_and_rng': True,
        'serialized_two_field_memory_continues_from_suffix_only': True,
        'serialized_checkpoint_bytes': len(encoded.encode()),
        'source_digest': source.source_digest,
        'sampling_claim': 'fixed seeded executions, not statistical confidence bounds'}


def _interval(value):
    return [str(dyadic(value[0])), str(dyadic(value[1], upper=True))]


def export(result, spatial):
    result = validate(result)
    lowers = tuple(dyadic(p) for p in result.lower_mass)
    rows = [{'initial_bits': list(initial), 'first_field': first, 'final_field': last,
             'scheduled_ticks': done, 'status': status,
             'probability_interval': [str(lo), str(dyadic(hi, upper=True))]}
            for (initial, first, last, done, status), lo, hi in zip(
                result.records, lowers, result.upper_mass, strict=True)]
    stats = statistics(result, spatial=spatial)
    return {'feedback': result.feedback, 'records': rows,
        'capped_complete_readout_l1_allowance': str(dyadic(1-sum(lowers), upper=True)),
        'eventual_readout_missing_mass_allowance': str(dyadic(1-sum(
            (p for r, p in zip(result.records, lowers, strict=True) if r[-1] == 'complete'), F(0)), upper=True)),
        'interval_cap_timeout_allowance': str(dyadic(result.timeout, upper=True)),
        'either_full_screen_eventual_field_means': [_interval(row) for row in stats['means']],
        'either_full_screen_eventual_joint_four': _interval(stats['joint_four']),
        'either_full_screen_eventual_connected_four': _interval(stats['connected_four'])}


def run():
    ports, edges, incident, observed, exits, _ = clocks.geometry()
    carrier = dual.source_carrier_certificate(ports, edges, incident)
    lower = carrier['unnormalized_incident_transfer_lower']
    walk_counts = clocks.exit_counts(ports, edges, observed, exits, 48)
    spatial, _ = clocks.exit_bound(walk_counts, 4, len(incident))
    baseline = law(ports, edges, incident)
    controlled = law(ports, edges, incident, feedback=True)
    first, second = (statistics(row, spatial=spatial) for row in (baseline, controlled))
    response = _subtract(second['means'][3], first['means'][3])
    if not (first['connected_four'][0] > 0 and second['connected_four'][1] < 0
            and response[0] > 0):
        raise ArithmeticError('the declared connected-field and feedback witnesses are unresolved')
    histories = [history(seed, feedback=feedback) for feedback in (False, True) for seed in range(8)]
    cutoff = clocks.slot_timeout_bound(clock_histories.MAX_CAP, 4, 33, clock_histories.source().slots)
    return {'schema': 'oph.observer_time_fields.refinement.v1',
        'source': {'levels': [3, 5], 'global_ports': list(clock_histories.source().ports),
            'matched_ports': len(ports), 'retained_seams': len(edges), 'incident_seams': len(incident),
            'clock_region_ports': 12, 'current_record_support_ports': len(observed),
            'preparation': 'independent Bernoulli 1/3 on every binary source port',
            'readout': ['initial n0', 'initial n1', 'first scheduled tick n0', 'fourth scheduled tick n1'],
            'feedback': 'first field one requests one additional native seam(0,1) repair',
            'complete_capped_records_per_protocol': len(baseline.records),
            'one_label_states': len(ports), 'two_label_states': math.comb(len(ports), 2),
            'retained_interval_attempt_cap': 512, 'integer_grid_bits': 52,
            'export_grid_bits': 80, 'maximum_native_attempts_per_capped_ball_history': 2049,
            'maximum_native_feedback_attempts': 1,
            'full_source_joint_configuration_space_enumerated': False},
        'positive_effective_field_connection': {
            'stationary_boundary': 'Omega_rho(x)=sqrt(rho**N(x)*(1-rho)**(V-N(x)))',
            'transfer': 'T=(I-A)^(-1/2) B (I-A)^(-1/2)',
            'Hamiltonian': '-log(T) in scheduled incident-attempt units',
            'all_local_field_and_internal_feedback_words_preserved': True,
            'binary_carrier_certificate': carrier,
            'strict_transfer_lower_common_ball': str(F(lower, len(edges))),
            'strict_transfer_lower_full_screens': [str(F(lower, clocks.screen.screen(level).seams))
                                                  for level in (3, 5)],
            'uniform_global_Hamiltonian_norm_or_all_level_limit_proved': False},
        'one_full_screen_to_common_ball_complete_output_allowance': str(dyadic(spatial, upper=True)),
        'two_screen_complete_output_allowance': str(dyadic(2*spatial, upper=True)),
        'field_protocols': [export(row, spatial) for row in (baseline, controlled)],
        'either_full_screen_feedback_increase_of_final_field': _interval(response),
        'full_source_histories': histories,
        'forced_external_timeout': history(1, cap=1),
        'full_source_external_cutoff_allowance': str(dyadic(cutoff, upper=True)),
        'scope': {'finite_density_full_readout_not_a_moment_closure': True,
            'same_native_source_preparation_clock_and_feedback_for_field_refinement': True,
            'C_and_E_share_this_diagonal_field_law': True,
            'uncoarsened_current_words_or_terminal_quantum_law_enumerated': False,
            'physical_duration_calibration_or_relativistic_continuum_proved': False},
        'source_bindings': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', type=Path); modes.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Finite-density incident-clock field refinement replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
