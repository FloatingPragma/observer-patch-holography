"""Replay the native contact response, its complete record map and leakage."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path

from tools import strict_json
from model.native_correlation_sectors import packet as sector_packet
from model.native_correlation_preparation import packet as preparation_packet
from model.native_exchange_clock import packet as clock_packet
from model.native_spatial_scaling import packet as spatial_packet
from . import contact as c, leakage, records

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple(sorted({
    'model/native_pair_response/PROOF.md', 'model/native_pair_response/response.py',
    'model/native_effective_fields/construction.py', 'model/native_effective_fields/README.md',
    'tools/strict_json.py', 'model/requirements.txt',
    *(f'model/native_contact_response/{name}' for name in (
        '__init__.py', 'contact.py', 'leakage.py', 'records.py', 'packet.py',
        'test_contact.py', 'test_leakage.py', 'test_records_packet.py', 'PROOF.md', 'README.md'))}))


def source_bindings():
    names = set(SOURCE_FILES)
    for owner in (sector_packet, preparation_packet, clock_packet, spatial_packet):
        names.update(owner.source_bindings())
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(names) if not name.startswith('docs/')}


def rational_tree(value):
    if type(value) is F:
        return str(value)
    if type(value) in (tuple, list):
        return [rational_tree(x) for x in value]
    if type(value) is dict:
        return {k: rational_tree(v) for k, v in value.items()}
    return value


def outward(interval, scale=10**15):
    lower, upper = interval
    return (F((lower*scale).numerator//(lower*scale).denominator, scale),
            F(-((-upper*scale).numerator//(-upper*scale).denominator), scale))


def produce():
    before = source_bindings()
    spectral = c.spectral_data(); pole = c.pole_certificate()
    for key in ('pole_interval', 'normalized_residue_interval',
                'native_attempt_transfer_interval', 'raw_pole_weight_interval'):
        pole[key] = outward(pole[key])
    controls = []
    for name, edges in (('path4', ((0, 1), (1, 2), (2, 3))),
                        ('complete4', tuple(combinations(range(4), 2)))):
        graph = c.source.source.Family((1,)*4, edges, (1,)*len(edges))
        for z in (F(-1), F(), F(9, 2)):
            controls.append({'graph': name, 'z': z, **c.contact_identity(graph, z)})
    result = {'schema': 'oph.native_contact_response.v1', 'source_bindings': before,
              'source': {'ports': 12, 'seams': 30, 'occupation': 2,
                         'port_weights': [1]*12, 'repair_probability': '1',
                         'reference': 'I66/66', 'attempt_transfer': 'I-G2/60'},
              'exact_response': {
                  'probe_seams': c.PROBE, 'primitive_dimension': 54,
                  'observed_cyclic_dimension': 12,
                  'numerator_descending': [int(x) for x in spectral['numerator'].all_coeffs()],
                  'denominator_descending': [int(x) for x in spectral['denominator'].all_coeffs()],
                  'primitive_characteristic_factors': [
                      {'coefficients': [int(x) for x in p.all_coeffs()], 'multiplicity': m}
                      for p, m in c.FACTORS],
                  'normalized_exchange_moments': spectral['moments'], 'lowest_pole': pole},
              'contact_resolvent_controls': controls,
              'complete_native_record_attachment': records.certificate(),
              'federation_leakage': leakage.certificate(),
              'endpoint_counterexample': {
                  'exchange_parameter': F(1, 5760),
                  'isolated_mode_heat_error_lower': leakage.endpoint_lower(F(1, 5760)),
                  'native_bound_requires_t_equals_n_over_2M': True},
              'scope': {'contact_coefficient_fitted': False,
                        'finite_comparison_threshold_is_infinite_volume_scattering_threshold': False,
                        'correlated_mode_is_a_localized_physical_particle': False,
                        'connected_federation_has_the_isolated_cell_eigenvectors': False,
                        'coherent_transfer_vacuum_or_unitary_clock_prepared': False,
                        'physical_particle_masses_or_relativistic_continuum_derived': False,
                        'native_caps_changed': False}}
    if source_bindings() != before:
        raise ArithmeticError('native contact sources changed during replay')
    return rational_tree(result)


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('native contact receipt differs from exact source replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true'); group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce())); print('Native contact response receipt written.')
    else:
        verify(); print('Native contact response receipt verified.')
