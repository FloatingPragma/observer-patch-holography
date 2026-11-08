"""Deterministic source-bound certificate for a common native experiment clock."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.native_rooted_refinement import packet as rooted
from model.native_exchange_clock import packet as exchange_owner
from model.native_fair_local_preparation import packet as preparation_owner
from model.native_fair_local_preparation import preparation
from model.native_fair_local_readout import packet as readout_owner
from model.native_phase_experiment import provenance
from . import clock as c, history as h, composition as j

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_fair_local_clock/'+name for name in (
    '__init__.py', 'clock.py', 'history.py', 'composition.py', 'packet.py', 'test_clock.py',
    'test_composition.py', 'PROOF.md', 'COMPOSITION.md', 'README.md')) + (
    'model/native_phase_experiment/README.md', 'tools/strict_json.py', 'model/requirements.txt')


def source_bindings():
    paths = (set(SOURCE_FILES) | set(rooted.source_bindings()) |
             set(exchange_owner.source_bindings()) | set(preparation_owner.source_bindings()) |
             set(readout_owner.source_bindings()) |
             set(provenance.source_hashes()))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) in (tuple, list):
        return [encode(x) for x in value]
    if type(value) is dict and all(type(k) is str for k in value):
        return {k: encode(v) for k, v in value.items()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise TypeError('exact JSON values required')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def composed_certificate():
    def compact(value):
        value = dict(value); plan = value.pop('plan')
        value.pop('detector_source_plan', None)
        value['source_programme'] = {
            'programme_sha256': plan['program_digest'], 'level': plan['level'],
            'declared_global_total': plan['total'], 'partition': plan['partition'],
            'epochs': plan['epochs'], 'attempts': plan['total_attempts'],
            'fixed_attempts': plan['fixed_attempts'], 'uniform_attempts': plan['uniform_attempts'],
            'resources': plan['program']['declared_resources']}
        return value
    return {'complete_fair_E_unitality': j.unitality_certificate(),
            'canonical_reference_then_uniform_readout': compact(j.canonical_protocol()),
            'sparse_two_source_reference_then_exact_tree':
                [compact(j.sparse_protocol(level)) for level in (0, 2, 6, 14, 200)]}


def produce():
    before = source_bindings()
    schedules = []
    for level in (0, 2, 9):
        plan = c.prepared_readout_plan(level, 240*4**level//3)
        row = {'level': level, 'programme_sha256': plan['program_digest'],
               'counted_partition': plan['partition'], 'epochs': plan['epochs'],
               'total_attempts': plan['total_attempts'],
               'AST_nodes': plan['program']['declared_resources']['AST_nodes'],
               'fixed_preparation_attempts': plan['fixed_attempts'],
               'coverage_budget_given_erasure_one_tenth':
                   c.readout_coverage_budget(plan, 1, 11, F(1, 10)),
               'preparation_reference_error_after_mixing_archive_marginal': F(1, 2**32),
               'expanded_host_completion_asserted': False}
        schedules.append(row)
    refused = False
    try:
        c.prepared_readout_plan(10, 240*4**10//3)
    except ValueError as error:
        refused = ('resource/domain boundary' in str(error) and
                   8*630*4**10 > c.native.MAX_REPEAT and
                   12*630*4**10 > c.native.MAX_ATTEMPTS)
    if not refused:
        raise ArithmeticError('the common fixed-duration programme exceeded its original budget')
    windows = []; attached = []
    preparation_word = tuple((edge, 'E') for edge in preparation.edges(32))
    for level in (6, 10, 14):
        window = c.routing.routed_window(level, 4)
        plan = c.window_plan(window, 240*4**level//3, (F(0), F(1, 630)),
                              (preparation_word,), ('E',))
        observed = tuple(e for e in window['internal_seams'] if set(e) <= set(window['root_ports']))
        epsilon = c.canonical_error(plan, window, F(1, 3))
        row = c.attach_window(plan, window, 0, observed, window['root_ports'], initial_error=epsilon,
            initial_reference='Bernoulli(1/3) on the canonical window tensor the same independent R')
        windows.append({'level': level, 'window_ports': len(window['canonical_ports']),
                        'window_seams': len(window['internal_seams']), 'attachment': row,
                        'fixed_preparation_attempts': plan['fixed_attempts'],
                        'native_programme_sha256': plan['program_digest'],
                        'expanded_high_level_history_executed': False})
        attached.append((plan, window, row))
    if attached[-1][2]['budget']['native_to_common_finite_window_upper'] >= F(1, 80000):
        raise ArithmeticError('the complete prepared-window comparison exceeds its displayed gate')
    pair = c.compare_attachments(*attached[1], *attached[2])
    result = {'schema': 'oph.native_fair_local_clock.v1', 'source_bindings': before,
        'shared_preparation_and_readout_schedules': schedules,
        'same_twelve_time_units_level10_admission_rejected': refused,
        'authenticated_finite_window_comparisons': windows,
        'level10_level14_complete_local_output_upper': pair,
        'actual_complete_raw_instrument_control': h.complete_fixed_control(),
        'source_admission_and_service_controls': h.source_domain_controls(),
        'composed_preparation_and_readout': composed_certificate(),
        'scope': {'one_actual_native_attempt_counter': True,
            'fixed_preparation_and_readout_attempts_charged': True,
            'complete_selected_local_CQ_output_with_untouched_reference': True,
            'source_retains_all_original_records': True,
            'observation_traces_unselected_archive_values_and_exact_ticks': True,
            'preparation_archive_independence_asserted': False,
            'coverage_erasure_premise_inferred_from_arbitrary_input': False,
            'canonical_spatial_product_specialization_has_independent_external_reference': True,
            'canonical_readout_composition_allows_correlated_reference': True,
            'sparse_readout_composition_requires_canonical_input': False,
            'preparation_mixing_archive_factored_as_independent': False,
            'literal_growing_coherent_load_produced': False,
            'coherent_service_impossibility_asserted': False,
            'source_pulse_to_repair_attempt_rate_selected': False,
            'all_radius_routed_continuation_or_unbounded_host_proved': False,
            'physical_time_density_relativistic_action_or_particles_selected': False,
            'canonical_physical_premise_discharged': False}}
    if before != source_bindings():
        raise ArithmeticError('source bytes changed during replay')
    return encode(result)


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('native clock receipt differs from exact typed replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Native fair local clock receipt written.')
    else:
        verify()
        print('Native fair local clock receipt verified.')
