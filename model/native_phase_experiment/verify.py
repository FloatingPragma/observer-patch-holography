"""Validate the complete finite experiment profile and exact ideal targets."""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import product
import math

from . import experiment
from ..observer.source_discrimination import circuit
from ..observer.scalar_response.experiment import compare_receipt


TOLERANCE = Fraction(str(experiment.TOLERANCE))


def _exact_close(actual, expected, label, tolerance=TOLERANCE):
    if abs(actual - Fraction(expected)) > tolerance:
        raise ValueError('independent comparison failed at ' + label)


def _close(actual, expected, label, tolerance=experiment.TOLERANCE):
    if type(actual) not in (int, float) or not math.isfinite(actual):
        raise ValueError('independent comparison failed at ' + label)
    _exact_close(Fraction(actual), expected, label, Fraction(str(tolerance)))


def _require(condition, label):
    if not condition:
        raise ValueError('independent structure mismatch at ' + label)


def _keys(actual, expected, label):
    _require(type(actual) is dict and actual.keys() == set(expected), label)


def _list(actual, length, label):
    _require(type(actual) is list and len(actual) == length, label)


def _probability(actual, expected, label, tolerance=experiment.TOLERANCE):
    _require(type(actual) is float and math.isfinite(actual) and actual >= 0,
             label + ' finite probability')
    _close(actual, expected, label, tolerance)


def _components(actual, carrier, q):
    """Reconstruct slot labels from the independent carrier, not the producer."""
    slots = defaultdict(list)
    index = 0
    for a, b in carrier['seams']:
        choices = [('point', tuple(i for i, bit in enumerate(state) if bit))
                   for state in carrier['states'] if state[a] == state[b]]
        choices += [('pair', (a, b, spectator)) for spectator in range(12)
                    if spectator not in (a, b)]
        _require(len(choices) == 56, 'independent preparation basis')
        for key in choices:
            slots[key].append(index)
            index += 1
    _require(index == 1680, 'independent preparation slot count')
    _require(len(slots) == 366, 'independent preparation component count')
    _list(actual, len(slots), 'preparation components')
    for item, ((kind, ports), indices) in zip(actual, sorted(slots.items()), strict=True):
        _keys(item, ('kind', 'ports', 'weight', 'selector_slots'), 'component fields')
        compare_receipt({'kind': kind, 'ports': list(ports), 'selector_slots': indices},
                        {key: item[key] for key in ('kind', 'ports', 'selector_slots')})
        _probability(item['weight'], q * len(indices) / 1680, 'component weight')


def _density(actual, entries, q, label):
    """Compare the exact squared Frobenius error with the unchanged norm bound."""
    _keys(actual, ('real', 'imag', 'outside_total_two_probability'), 'density fields')
    _probability(actual['outside_total_two_probability'], 0, label + ' leakage')
    target = {(row, column): q * Fraction(value) for row, column, value in entries}
    for field in ('real', 'imag'):
        _list(actual[field], 66, 'prepared density rows')
        for row in actual[field]:
            _list(row, 66, 'prepared density columns')
            _require(all(type(value) is float and math.isfinite(value) for value in row),
                     'prepared density finite entries')
    squared_error = Fraction(0)
    for row, column in product(range(66), repeat=2):
        real_error = Fraction(actual['real'][row][column]) - target.get((row, column), 0)
        imaginary = Fraction(actual['imag'][row][column])
        squared_error += real_error * real_error + imaginary * imaginary
        if squared_error > TOLERANCE * TOLERANCE:
            raise ValueError('independent comparison failed at ' + label + ' full density')


def _resources(runs, selector):
    bits = selector['word_bits'] * selector['attempts']
    preparation = max(len(runs[kind].word) for kind in ('point', 'pair'))
    phase = max(len(runs[setting + '_phase'].word) for setting in ('reference', 'twisted'))
    pulses = bits * len(runs['selector'].word)
    return {'maximum_selector_source_pulses': pulses,
            'maximum_selector_measurements': bits,
            'terminal_computational_measurements': 2,
            'maximum_trial_source_pulses': {
                law: pulses + preparation + phase + len(runs[law].word)
                + len(runs['analyzer'].word) for law in ('C', 'E')},
            'proper_code_load_export_events_per_source_pulse': 2}


def verify(result, expected):
    """Check every profile field using the trusted exact oracle and source words."""
    _keys(result, ('schema', 'target', 'pointer_records', 'selector', 'preparation',
                   'source_execution', 'source_sha256', 'protected_obstruction_reuse',
                   'trials', 'resource_envelope'), 'receipt fields')
    _require(result['schema'] == 'oph.native_phase_experiment.v1', 'schema')
    _require(expected['compact'] is False and expected['checks']['all_pass'] is True,
             'complete independent expectation')
    _require(result['source_sha256'] == experiment._dependencies(), 'source closure')
    compare_receipt({'occupation_bits': 12, 'binary_total': 2, 'phase_port': 0,
                     'analyzer_ports': [0, 1], 'seams': expected['carrier']['seams']}, result['target'])
    pointers = result['pointer_records']
    _list(pointers, 36, 'pointer record count')
    for pointer, (first, second) in zip(pointers, product(range(6), repeat=2), strict=True):
        compare_receipt({'first': first, 'second': second,
                         'abort': first >= 2 or second >= 2,
                         'score': int((first, second) == (1, 0))
                         - int((first, second) == (0, 1))}, pointer)
    selection = result['selector']
    _keys(selection, ('bits_per_draw', 'maximum_draws', 'accepted_slots', 'retry_slots',
                      'accepted_probability', 'source_failure_probability', 'timeout_probability',
                      'normalization_error', 'bit_record_probabilities', 'slot_probabilities',
                      'archive'), 'selector fields')
    selector_expected = expected['finite_selector']
    for actual, wanted in (('bits_per_draw', 'word_bits'), ('maximum_draws', 'attempts'),
                           ('accepted_slots', 'accepted_word_count'), ('retry_slots', 'rejected_word_count')):
        _require(type(selection[actual]) is int and selection[actual] == selector_expected[wanted],
                 'selector ' + actual)
    q = Fraction(selector_expected['accepted_probability'])
    _probability(selection['accepted_probability'], q, 'selector accepted weight')
    _probability(selection['timeout_probability'], Fraction(selector_expected['timeout_probability']),
                 'selector timeout', 2e-20)
    _probability(selection['source_failure_probability'], 0, 'selector source failure')
    _list(selection['bit_record_probabilities'], 6, 'selector bit records')
    for label, value in enumerate(selection['bit_record_probabilities']):
        _probability(value, Fraction(1, 2) if label < 2 else 0, 'selector bit record')
    _list(selection['slot_probabilities'], 1680, 'selector slots')
    for index, value in enumerate(selection['slot_probabilities']):
        _probability(value, q / 1680, 'selector slot ' + str(index))
    slot_mass = sum(map(Fraction, selection['slot_probabilities']), Fraction(0))
    _exact_close(slot_mass, q, 'complete selector accepted weight')
    selector_balance = (slot_mass + Fraction(selection['source_failure_probability'])
                        + Fraction(selection['timeout_probability']))
    _exact_close(selector_balance, 1, 'complete selector probability')
    _close(selection['normalization_error'], abs(selector_balance - 1), 'selector balance diagnostic')
    _require(selection['archive'] == 'literal bit records, retry indices, first source failure or timeout',
             'selector archive description')
    preparation = result['preparation']
    _keys(preparation, ('components', 'source_failure_probability', 'success_probability', 'densities'),
          'preparation fields')
    _components(preparation['components'], expected['carrier'], q)
    _probability(preparation['success_probability'], q, 'preparation unconditioned success')
    _probability(preparation['source_failure_probability'], 0, 'preparation source failure')
    _keys(preparation['densities'], ('reference', 'twisted'), 'preparation settings')
    for setting, name in (('reference', 'rho'), ('twisted', 'sigma')):
        _density(preparation['densities'][setting], expected['preparation'][name]['entries'], q, setting)
    runs = experiment.source_runs()
    _keys(result['source_execution'], runs, 'source execution names')
    for name, run in runs.items():
        compare_receipt(experiment._run_summary(run), result['source_execution'][name])
    compare_receipt(_resources(runs, selector_expected), result['resource_envelope'])
    lineage = result['protected_obstruction_reuse']
    compare_receipt(experiment.classify_phase_loss(experiment.ROOT), lineage)
    for law in ('C', 'E'):
        _require(lineage[law]['survival'] == expected['one_repair'][law]['delta_factor'], 'phase quotient')
        _require(lineage[law]['classification']['almostSure'] is True
                 and lineage[law]['classification']['closedTrap'] is False, 'phase first loss')
    _require(lineage['hold_control']['almostSure'] is False
             and lineage['hold_control']['closedTrap'] is True, 'hold control')
    _list(result['trials'], 6, 'six trials')
    for trial in result['trials']:
        _keys(trial, ('setting', 'law', 'records', 'score', 'measured_probability',
                      'phase_source_failure_probability', 'repair_source_failure_probability',
                      'analyzer_source_failure_probability', 'complete_probability_error',
                      'regional_factorization_error'), 'trial fields')
    trials = {(trial['setting'], trial['law']): trial for trial in result['trials']}
    _require(set(trials) == {(setting, law) for setting in ('reference', 'twisted')
                            for law in ('before', 'C', 'E')}, 'unique trials')
    max_probability_error = max_readout_error = max_pointer_error = 0.
    checked_raw = checked_pointers = 0
    for (setting, law), trial in trials.items():
        name = 'rho' if setting == 'reference' else 'sigma'
        exact = expected['preparation'] if law == 'before' else expected['one_repair'][law]
        _close(trial['score'], q * Fraction(exact[name + '_readout']), setting + ' ' + law + ' score')
        if law == 'before':
            expected_rows = {((0, 1), None, None): {
                name + '_probability': '1', name + '_readout_weight': exact[name + '_readout'],
                name + '_pointer_probabilities': exact[name + '_pointer_probabilities']}}
        else:
            expected_rows = {(tuple(row['seam']), row['coin'], row['raw']): row
                             for row in exact['raw_branches']}
        records = trial['records']
        _require(type(records) is list, 'raw records')
        for row in records:
            _keys(row, ('seam', 'coin', 'raw', 'raw_abort', 'outcomes'), 'raw fields')
            _list(row['seam'], 2, 'raw seam')
            _require(all(type(port) is int for port in row['seam']), 'raw seam integer labels')
            _require((row['coin'] is None and row['raw'] is None) if law == 'before'
                     else type(row['coin']) is int and type(row['raw']) is int, 'raw integer labels')
        keys = [(tuple(row['seam']), row['coin'], row['raw']) for row in records]
        _require(len(keys) == len(set(keys)) and set(keys) == set(expected_rows), 'complete raw grid')
        measured = scored = Fraction(0)
        for key, row in zip(keys, records, strict=True):
            wanted = expected_rows[key]
            _require(type(row['raw_abort']) is bool, 'raw failure type')
            _require(row['raw_abort'] == (False if law == 'before' else
                     row['raw'] not in circuit.record_values(law, row['coin'])), 'raw failure label')
            outcomes = row['outcomes']
            _list(outcomes, 36, 'complete pointer grid')
            wanted_pointers = wanted[name + '_pointer_probabilities']
            _require(len(wanted_pointers) == 36, 'independent pointer grid')
            probability = score = Fraction(0)
            for value, wanted_pointer, pointer in zip(outcomes, wanted_pointers, pointers, strict=True):
                _require(type(value) is float and math.isfinite(value) and value >= 0,
                         'pointer probability')
                number = Fraction(value)
                probability_target = q * Fraction(wanted_pointer)
                _exact_close(number, probability_target, 'individual pointer ' + str(key))
                max_pointer_error = max(max_pointer_error, float(abs(number - probability_target)))
                probability += number
                score += number * pointer['score']
                checked_pointers += 1
            wanted_probability = q * Fraction(wanted[name + '_probability'])
            wanted_contribution = q * Fraction(wanted[name + '_readout_weight'])
            _exact_close(probability, wanted_probability, 'raw probability ' + str(key))
            _exact_close(score, wanted_contribution, 'raw readout ' + str(key))
            max_probability_error = max(max_probability_error, float(abs(probability - wanted_probability)))
            max_readout_error = max(max_readout_error, float(abs(score - wanted_contribution)))
            measured += probability
            scored += score
            checked_raw += law != 'before'
        _probability(trial['measured_probability'], q, 'measured target')
        _close(trial['measured_probability'], measured, 'measured total')
        _close(trial['score'], scored, 'scored total')
        failure_values = [selection['source_failure_probability'], selection['timeout_probability'],
                          preparation['source_failure_probability']]
        for field in ('phase_source_failure_probability', 'repair_source_failure_probability',
                      'analyzer_source_failure_probability'):
            _probability(trial[field], 0, field)
            failure_values.append(trial[field])
        failures = sum(map(Fraction, failure_values), Fraction(0))
        _exact_close(measured + failures, 1, 'complete probability')
        _close(trial['complete_probability_error'], abs(measured + failures - 1), 'balance diagnostic')
        _close(trial['regional_factorization_error'], 0, 'regional marginal')
    return {'all_pass': True, 'independent_checks': expected['checks']['count'],
            'raw_records_compared': checked_raw, 'pointer_records_compared': checked_pointers,
            'maximum_raw_probability_error': max_probability_error,
            'maximum_raw_readout_error': max_readout_error,
            'maximum_pointer_probability_error': max_pointer_error}
