"""Source-bound replay of complete local history capacity allowances."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.native_local_response import packet as owning
from model.native_spatial_scaling import geometry
from . import limit as l

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_local_history_limit/'+name for name in (
    '__init__.py', 'limit.py', 'packet.py', 'test_limit.py', 'PROOF.md', 'README.md')) + (
    'model/native_capacity_restriction/restriction.py',
    'model/native_capacity_restriction/PROOF.md',
    'model/native_capacity_restriction/__init__.py',
    'model/native_causal_dataflow/__init__.py',
    'model/native_record_causality/__init__.py',
    'model/native_spatial_scaling/__init__.py',
    'model/native_spatial_scaling/geometry.py',
    'model/observer_time_fields/canonical_sector.py',
    'model/observer_time_fields/dual.py',
    'model/observer_time_fields/refinement.py',
    'model/weighted_source_family/__init__.py',
    'model/federated_observer/SCALING.md',
    'model/native_exchange_clock/PROOF.md')


def source_bindings():
    paths = set(SOURCE_FILES) | set(owning.source_bindings())
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {k: encode(v) for k, v in value.items()}
    if type(value) in (tuple, list):
        return [encode(v) for v in value]
    if value is None or type(value) in (int, str, bool):
        return value
    raise TypeError('exact JSON values required')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def produce():
    before = source_bindings()
    attachment = owning.source_attachment()
    partition = (F(0), F(1, 126), F(1, 63))
    signature = None
    rows = []
    for level in (0, 1, 2):
        halo = l.halo_signature(level)
        key = (halo['internal_edges'], halo['attachment_ports'])
        if signature is None:
            signature = key
        if key != signature or not halo['endpoint_orientation_preserved']:
            raise ArithmeticError('actual source halos do not have the declared record isomorphism')
        census = geometry.census(level)
        v, m = census['ports'], census['seams']
        attempts = l.aligned_clock(m, partition)[-1]
        program = l.cone.source.program(federation={'level': level, 'total': v//3}, biased=False,
            commands=[{'op': 'repair', 'selector': {'kind': 'uniform'},
                       'law': 'mixture', 'repeat': attempts}], capacity=1, overflow='export_oldest')
        l.cone.audit_program(program)
        state = l.canonical_state_bound(v, v//3, len(halo['window']), F(1, 3))
        selection = l.poissonization_bound(m, len(halo['internal_edges']), partition)
        # Exact selector-only cone law for this source radius-one halo.
        exact = l.bounds.halo_global_error(m, attempts)['cone_escape_probability']
        row = {'level': level, 'ports': v, 'seams': m, 'global_occupation': v//3,
               'window_ports': len(halo['window']), 'window_seams': len(halo['internal_edges']),
               'source_address_isomorphism': halo, 'clock_partition': partition,
               'native_attempts': attempts, 'program_digest': l.cone.source.digest(program),
               'initial_state_upper': state, 'exact_native_halo_escape': exact,
               'complete_window_instrument_upper': selection,
               'native_to_common_window_upper': min(F(1), state+exact+selection),
               'selector_control': l.selector_law_control(m, 33, partition),
               'expanded_native_history_execution_claimed': False}
        rows.append(row)
    formal = geometry.census(200)
    # Census-only examples are conditional on a matching window. They do
    # not assert that the floating graph constructor or program admits it.
    census_budget = l.local_budget(formal['ports'], formal['ports']//3,
        formal['seams'], 15, 33, 12, 6, 1, partition, F(1, 3))
    result = {'schema': 'oph.native_local_history_limit.v1', 'source_bindings': before,
        'native_endpoint_attachment': attachment,
        'actual_matching_halo_controls': rows,
        'pairwise_complete_observation_bounds': [
            {'levels': [a, b], 'half_trace_upper': min(F(1),
             rows[a]['native_to_common_window_upper']+rows[b]['native_to_common_window_upper'])}
            for a, b in ((0, 1), (0, 2), (1, 2))],
        'conditional_level200_matching_window_bound': {
            'census': formal, 'budget': census_budget,
            'formal_kernel_attempts': l.aligned_clock(formal['seams'], partition)[-1],
            'native_attempt_cap': l.cone.source.MAX_ATTEMPTS,
            'all_level_geometry_or_native_admission_asserted': False},
        'complete_Poisson_window_overflow_upper': l.poisson_tail(F(33, 63), 6),
        'scope': {'complete_local_raw_records_and_quantum_disturbance': True,
            'visible_native_coins_retained': True, 'arbitrary_reference_channel_bound': True,
            'canonical_state_specialization_requires_independent_reference': True,
            'within_bin_event_order_retained': True,
            'exact_global_tick_addresses_retained': False,
            'remote_archive_erased_from_native_source': False,
            'uniform_fixed_attempt_scheduler_replaced_in_native_source': False,
            'complete_global_history_equivalence': False,
            'native_refinement_unique_local_limit_proved': False,
            'physical_capacity_density_clock_or_dimension_selected': False,
            'interacting_chiral_field_or_particle_identified': False,
            'new_canonical_physical_premise_discharged': False}}
    if before != source_bindings():
        raise ArithmeticError('source bytes changed during replay')
    return encode(result)


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('complete local history receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Complete local history receipt written.')
    else:
        verify()
        print('Complete local history receipt verified.')
