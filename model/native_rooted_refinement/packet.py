"""Replay exact topology, finite routing certificates and local-record bounds."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.native_local_history_limit import packet as owning
from . import refinement as r, routing, attachment as a

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_rooted_refinement/'+name for name in (
    '__init__.py', 'refinement.py', 'routing.py', 'attachment.py', 'packet.py',
    'test_refinement.py', 'PROOF.md', 'README.md')) + (
    'tools/strict_json.py', 'model/requirements.txt')


def source_bindings():
    paths = set(SOURCE_FILES) | set(owning.source_bindings()) | set(r.affine.SOURCE_PATHS)
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


def produce():
    before = source_bindings()
    topology_rows = []
    for radius in (0, 1, 2, 4, 8, 16, 32):
        level = r.minimum_level(radius)
        signature = r.infinite_ball(radius)
        if r.signature(level+1, radius) != signature or r.signature(level+3, radius) != signature:
            raise ArithmeticError('the exact addressed corner signatures differ')
        topology_rows.append({'cell_radius': radius, 'sufficient_level': level,
            'cells': len(signature),
            'signature_sha256': hashlib.sha256(canonical(encode(signature)).encode()).hexdigest()})
    windows = []
    for radius in (1, 2, 4, 8):
        first = None
        for level in (r.minimum_level(radius), 6, 10, 14):
            window = routing.routed_window(level, radius)
            if first is None:
                first = window
            elif not routing.same_window(first, window):
                raise ArithmeticError('certified finite routed windows differ')
            certificates = window['routing_certificates']
            windows.append({'level': level, 'cell_radius': radius,
                'window': window, 'summary': {
                    'level': level, 'cell_radius': radius,
                    'ports': len(window['canonical_ports']),
                    'internal_seams': len(window['internal_seams']),
                    'certified_cell_routing_rows': len(certificates),
                    'collision_rows': sum(row['collision_repaired'] for row in certificates),
                    'stable_tie_rows': sum(row['certified_score_gap'] == 0 for row in certificates),
                    'score_variation_allowance': routing.SCORE_ERROR,
                    'canonical_window_sha256': a.window_digest(window),
                    'source_address_sha256': hashlib.sha256(
                        canonical(encode(window['source_ports'])).encode()).hexdigest(),
                    'source_address_formula': '12*(original_face*4^level+child_offset)+port',
                    'all_source_endpoints_order_preserved': True,
                    'global_constructor_or_native_history_executed': False}})
    partition = (F(), F(1, 1260), F(1, 630))
    histories = [a.history_budget(row['window'], F(1, 3), partition)
                 for row in windows if row['cell_radius'] == 4]
    if histories[-1]['budget']['native_to_common_finite_window_upper'] >= F(1, 10**6):
        raise ArithmeticError('the displayed finite full-record bound exceeds its gate')
    output = {'schema': 'oph.native_rooted_refinement.v1', 'source_bindings': before,
        'exact_corner_topology': topology_rows,
        'literal_source_prefix_attachment': [a.source_prefix(level) for level in range(4)],
        'certified_matching_finite_windows': [row['summary'] for row in windows],
        'canonical_complete_local_history_bounds': histories,
        'pairwise_level10_level14_upper': min(F(1), sum(
            row['budget']['native_to_common_finite_window_upper']
            for row in histories if row['level'] in (10, 14))),
        'source_domain_controls': a.host_and_arithmetic_controls(),
        'scope': {
            'exact_oriented_stencil_all_radius_corner_full_sequence_proved': True,
            'literal_floating_orientation_all_level_agreement_proved': False,
            'finite_local_source_routing_certificates': True,
            'canonical_static_state_family_density_declared': True,
            'fine_native_records_and_reference_transport_preserved': True,
            'infinite_routed_full_sequence_limit_proved': False,
            'constructive_infinite_routed_subsequence_produced': False,
            'all_level_floating_constructor_totality_asserted': False,
            'high_level_expanded_native_history_execution_asserted': False,
            'arithmetic_continuation_selected': False,
            'physical_capacity_density_clock_or_dimension_selected': False,
            'canonical_physical_identification_premises_discharged': False}}
    if before != source_bindings():
        raise ArithmeticError('source bytes changed during replay')
    return encode(output)


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('rooted refinement receipt differs from source replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Rooted refinement receipt written.')
    else:
        verify()
        print('Rooted refinement receipt verified.')
