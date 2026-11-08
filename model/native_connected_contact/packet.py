"""Replay source-bound connected dressing and threshold certificates."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.native_contact_response import packet as contact_packet
from model.native_pair_response import packet as pair_packet
from . import dressing, threshold, records

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple(sorted({
    'tools/strict_json.py', 'model/requirements.txt',
    *(f'model/native_connected_contact/{name}' for name in (
        '__init__.py', 'dressing.py', 'threshold.py', 'records.py', 'packet.py',
        'test_dressing.py', 'test_threshold.py', 'test_records_packet.py', 'PROOF.md', 'README.md'))}))


def source_bindings():
    names = set(SOURCE_FILES)
    for owner in (contact_packet, pair_packet):
        names.update(owner.source_bindings())
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(names) if not name.startswith('docs/')}


def produce():
    before = source_bindings(); response = dressing.certificate()
    for row in response['finite_source_controls']:
        row['normalized_Abel_interval'] = contact_packet.outward(row['normalized_Abel_interval'])
    result = {'schema': 'oph.native_connected_contact.v1', 'source_bindings': before,
              'source': {'occupation': 2, 'port_weights': 'all one', 'repair_probability': '1',
                         'selector': 'independent uniform seam', 'reference': 'I/binom(V,2)',
                         'attempt_transfer': 'I-G/(2M)', 'ports': '240*4^level', 'seams': '630*4^level'},
              'connected_comparison_threshold': threshold.certificate(),
              'controlled_exterior_dressing': response,
              'complete_native_Abel_readout': records.certificate(),
              'scope': {'supplied_spatial_or_Higgs_action': False,
                        'positive_energy_dressed_particle_pole_proved': False,
                        'all_possible_bound_sectors_or_resonances_excluded': False,
                        'connected_threshold_collapse_from_vanishing_lower_bound': False,
                        'bare_cell_binding_threshold_valid_on_connected_family': False,
                        'physical_particle_mass_or_Lorentzian_QFT_derived': False,
                        'native_caps_changed': False}}
    if source_bindings() != before:
        raise ArithmeticError('connected contact sources changed during replay')
    return contact_packet.rational_tree(result)


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('connected contact receipt differs from exact source replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true'); group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce())); print('Connected native contact receipt written.')
    else:
        verify(); print('Connected native contact receipt verified.')
