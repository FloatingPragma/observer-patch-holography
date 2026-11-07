"""Compare complete source execution against separately reconstructed counts."""
from __future__ import annotations

from fractions import Fraction
from itertools import product
import math

import numpy as np

from . import experiment
from ..observer.source_discrimination import circuit
from ..observer.scalar_response.experiment import compare_receipt


def _close(actual, expected, label, tolerance=experiment.TOLERANCE):
    if (type(actual) not in (int, float) or not math.isfinite(actual)
            or abs(Fraction(actual) - Fraction(expected)) > Fraction(str(tolerance))):
        raise ValueError('independent comparison failed at ' + label)


def _require(condition, label):
    if not condition:
        raise ValueError('independent structure mismatch at ' + label)


def verify(result, expected):
    """Check every retained raw slot and the unconditioned measured contrast."""
    _require(result['schema'] == 'oph.native_phase_experiment.v1', 'schema')
    _require(expected['compact'] is False and expected['checks']['all_pass'] is True,
             'complete independent expectation')
    _require(result['source_sha256'] == experiment._dependencies(), 'source closure')
    target = result['target']
    _require(target == {'occupation_bits': 12, 'binary_total': 2, 'phase_port': 0,
                       'analyzer_ports': [0, 1], 'seams': expected['carrier']['seams']}, 'target')
    pointers = result['pointer_records']
    _require(len(pointers) == 36, 'pointer record count')
    for pointer, (first, second) in zip(pointers, product(range(6), repeat=2), strict=True):
        _require(pointer == {'first': first, 'second': second,
                             'abort': first >= 2 or second >= 2,
                             'score': int((first, second) == (1, 0))
                             - int((first, second) == (0, 1))}, 'pointer contract')
    selection = result['selector']
    selector_expected = expected['finite_selector']
    for actual, wanted in (('bits_per_draw', 'word_bits'), ('maximum_draws', 'attempts'),
                           ('accepted_slots', 'accepted_word_count'), ('retry_slots', 'rejected_word_count')):
        _require(type(selection[actual]) is int and selection[actual] == selector_expected[wanted],
                 'selector ' + actual)
    q = Fraction(selector_expected['accepted_probability'])
    _close(selection['accepted_probability'], q, 'selector accepted weight')
    _close(selection['timeout_probability'], Fraction(selector_expected['timeout_probability']),
           'selector timeout', 2e-20)
    _require(len(selection['slot_probabilities']) == 1680, 'selector slots')
    for index, value in enumerate(selection['slot_probabilities']):
        _close(value, q / 1680, 'selector slot ' + str(index))
    _require(len(result['preparation']['components']) == 366, 'preparation components')
    for name, run in experiment.source_runs().items():
        compare_receipt(experiment._run_summary(run), result['source_execution'][name])
    lineage = result['protected_obstruction_reuse']
    for law in ('C', 'E'):
        _require(lineage[law]['survival'] == expected['one_repair'][law]['delta_factor'], 'phase quotient')
        _require(lineage[law]['classification']['almostSure'] is True
                 and lineage[law]['classification']['closedTrap'] is False, 'phase first loss')
    _require(lineage['hold_control']['almostSure'] is False
             and lineage['hold_control']['closedTrap'] is True, 'hold control')
    for setting, name in (('reference', 'rho'), ('twisted', 'sigma')):
        actual = result['preparation']['densities'][setting]
        _close(actual['outside_total_two_probability'], 0., setting + ' preparation leakage')
        wanted = np.zeros((66, 66))
        for row, column, value in expected['preparation'][name]['entries']:
            wanted[row, column] = float(Fraction(value))
        density = np.asarray(actual['real']) + 1j * np.asarray(actual['imag'])
        _require(density.shape == wanted.shape and np.isfinite(density).all(), 'prepared density')
        _close(float(np.linalg.norm(density - float(q) * wanted)), 0., setting + ' full density')
    _require(len(result['trials']) == 6, 'six trials')
    trials = {(trial['setting'], trial['law']): trial for trial in result['trials']}
    _require(set(trials) == {(setting, law) for setting in ('reference', 'twisted')
                            for law in ('before', 'C', 'E')}, 'unique trials')
    max_probability_error = max_readout_error = max_pointer_error = 0.
    checked_raw = 0
    checked_pointers = 0
    for (setting, law), trial in trials.items():
        name = 'rho' if setting == 'reference' else 'sigma'
        exact = expected['preparation'] if law == 'before' else expected['one_repair'][law]
        wanted_score = q * Fraction(exact[name + '_readout'])
        _close(trial['score'], wanted_score, setting + ' ' + law + ' score')
        if law == 'before':
            expected_rows = {((0, 1), None, None): {
                name + '_probability': '1', name + '_readout_weight': exact[name + '_readout'],
                name + '_pointer_probabilities': exact[name + '_pointer_probabilities']}}
        else:
            expected_rows = {(tuple(row['seam']), row['coin'], row['raw']): row
                             for row in exact['raw_branches']}
        records = trial['records']
        keys = [(tuple(row['seam']), row['coin'], row['raw']) for row in records]
        _require(len(keys) == len(set(keys)) and set(keys) == set(expected_rows), 'complete raw grid')
        measured = scored = 0.
        for key, row in zip(keys, records, strict=True):
            wanted = expected_rows[key]
            _require(type(row['raw_abort']) is bool, 'raw failure type')
            if law == 'before':
                _require(row['raw_abort'] is False, 'initial readout raw status')
            else:
                _require(row['raw_abort'] == (row['raw'] not in circuit.record_values(law, row['coin'])),
                         'raw failure label')
            outcomes = row['outcomes']
            _require(len(outcomes) == 36, 'complete pointer grid')
            wanted_pointers = wanted[name + '_pointer_probabilities']
            _require(len(wanted_pointers) == 36, 'independent pointer grid')
            for value, wanted_pointer in zip(outcomes, wanted_pointers, strict=True):
                _require(type(value) is float and math.isfinite(value) and value >= 0,
                         'pointer probability')
                probability_target = q * Fraction(wanted_pointer)
                _close(value, probability_target, 'individual pointer ' + str(key))
                max_pointer_error = max(max_pointer_error, abs(value - float(probability_target)))
                checked_pointers += 1
            probability = math.fsum(outcomes)
            score = math.fsum(value * pointer['score']
                             for value, pointer in zip(outcomes, pointers, strict=True))
            wanted_probability = q * Fraction(wanted[name + '_probability'])
            wanted_contribution = q * Fraction(wanted[name + '_readout_weight'])
            _close(probability, wanted_probability, 'raw probability ' + str(key))
            _close(score, wanted_contribution, 'raw readout ' + str(key))
            max_probability_error = max(max_probability_error, abs(probability - float(wanted_probability)))
            max_readout_error = max(max_readout_error, abs(score - float(wanted_contribution)))
            measured += probability
            scored += score
            checked_raw += law != 'before'
        _close(trial['measured_probability'], measured, 'measured total')
        _close(trial['score'], scored, 'scored total')
        failures = (selection['source_failure_probability'] + selection['timeout_probability']
                    + result['preparation']['source_failure_probability']
                    + trial['phase_source_failure_probability']
                    + trial['repair_source_failure_probability']
                    + trial['analyzer_source_failure_probability'])
        for value in (selection['source_failure_probability'], selection['timeout_probability'],
                      result['preparation']['source_failure_probability'],
                      trial['phase_source_failure_probability'], trial['repair_source_failure_probability'],
                      trial['analyzer_source_failure_probability']):
            _require(type(value) is float and math.isfinite(value) and value >= 0, 'failure probability')
        _close(measured + failures, 1., 'complete probability')
        _close(trial['complete_probability_error'], abs(measured + failures - 1.), 'balance diagnostic')
        _close(trial['regional_factorization_error'], 0., 'regional marginal')
    return {'all_pass': True, 'independent_checks': expected['checks']['count'],
            'raw_records_compared': checked_raw,
            'pointer_records_compared': checked_pointers,
            'maximum_raw_probability_error': max_probability_error,
            'maximum_raw_readout_error': max_readout_error,
            'maximum_pointer_probability_error': max_pointer_error}
