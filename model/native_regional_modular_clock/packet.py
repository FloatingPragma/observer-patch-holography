"""Replay the exact regional reference and authenticated finite collars."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

from model.native_rooted_refinement import packet as rooted
from model.primitive_observer import transfer
from tools import strict_json
from . import certificate, collar_certificate

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
OWNED_FILES = tuple('model/native_regional_modular_clock/' + name for name in (
    '__init__.py', 'regional.py', 'collar.py', 'certificate.py', 'collar_certificate.py',
    'packet.py', 'test_regional.py', 'test_collar.py', 'test_packet.py',
    'PROOF.md', 'COLLAR.md', 'README.md'))
ADDITIONAL_FILES = (
    'model/federated_observer/SCALING.md',
    'model/native_spatial_scaling/PROOF.md',
    'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean')


def source_bindings():
    paths = set(rooted.source_bindings()) | set(transfer.SOURCE_PATHS)
    paths.update(OWNED_FILES)
    paths.update(ADDITIONAL_FILES)
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in sorted(paths)}


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) in (tuple, list):
        return [encode(item) for item in value]
    if type(value) is dict and all(type(key) is str for key in value):
        return {key: encode(item) for key, item in value.items()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise TypeError('exact typed scientific values required')


def canonical(value):
    return json.dumps(encode(value), sort_keys=True, indent=2, allow_nan=False) + '\n'


def produce():
    before = source_bindings()
    native = certificate.certificate()
    windows = [collar_certificate.row(14, radius) for radius in (2, 4, 8)]
    bound = F(windows[1]['bounds']['regional_modular_to_larger_exchange_flow_upper'])
    if bound >= F(1068, 10**9):
        raise ArithmeticError('the authenticated radius-four collar lost its stated bound')
    value = {
        'schema': 'oph.native_regional_modular_clock.v1',
        'source_bindings': before,
        'same_reference_regional_certificate': native,
        'authenticated_finite_source_collars': windows,
        'limiting_contract': {
            'clock': 'rho_A^(i*2M*t) O rho_A^(-i*2M*t); Heisenberg exchange time -t',
            'observable_algebra': 'regional observables preserving local occupation number',
            'local_generator': '[G_A+(1-2p)N_boundary]/[1-p(1-p)]',
            'capacity_condition': 'v^2/M -> 0 for bounded degree and M proportional to V',
            'density_condition': 'v*abs(q/V-p) -> 0; exact or rounded fixed density suffices',
            'collar_condition': 'crossing_count*T_r(4*degree*J*abs(t)) -> 0',
            'routing_condition': 'one authenticated compatible rooted port graph or justified subsequence'},
        'scope': {
            'exact_same_transfer_reference_restriction': True,
            'rescaled_regional_modular_exchange_limit': True,
            'finite_interior_collar_bound': True,
            'fixed_region_modular_overlap_compatibility': False,
            'global_transfer_reference_scalable_preparation_supplied': False,
            'all_level_numerical_port_routing_proved': False,
            'controlled_unitary_equals_unread_repair': False,
            'physical_lorentz_boost_or_responsive_metric': False,
            'physical_G1_to_G6_discharged': False},
        'custody': {
            'fresh_owned_regional_and_collar_replay': True,
            'inherited_source_files_bound': True,
            'all_inherited_producers_replayed': False,
            'full_python_import_closure_claimed': False}}
    if source_bindings() != before:
        raise ArithmeticError('source files changed during the owning replay')
    return value


def check_payload(retained, fresh):
    if canonical(retained) != canonical(fresh):
        raise ArithmeticError('the retained regional receipt differs from a fresh owning replay')
    return True


def verify(path=OUTPUT):
    fresh = produce()
    check_payload(strict_json.load(Path(path)), fresh)
    return {'verified': True, 'native_transfer_entries': 4356,
            'native_regions': 30, 'authenticated_collars': 3,
            'source_files': len(fresh['source_bindings']),
            'physical_G1_to_G6_discharged': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        print(json.dumps(verify(), indent=2))
    else:
        value = produce()
        if args.write:
            OUTPUT.write_text(canonical(value), encoding='utf-8')
            check_payload(strict_json.load(OUTPUT), value)
        print(json.dumps({'source_files': len(value['source_bindings']),
                          'native_regions': 30, 'authenticated_collars': 3,
                          'physical_G1_to_G6_discharged': False}, indent=2))


if __name__ == '__main__':
    main()
