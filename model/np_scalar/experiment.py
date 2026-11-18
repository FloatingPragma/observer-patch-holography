"""Frozen-source replay of the finite native scalar construction.

Verification is read-only. Recording requires an explicit unused output path
and an exact production-source freeze. Scientific gates run under Python -O.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

from tools.pin_scope import bound, bound_pins
from .carrier import load_source
from .compiler import require

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def source_pins():
    paths = [p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py')]
    paths += ['model/np_repair.py', 'model/np_scalar/protocol.json',
              'model/polynomial_sv/carrier.py', 'model/polynomial_sv/law.py',
              'model/polynomial_sv/sectors.py',
              'model/np_export/inputs/CELL.npz', 'tools/pin_scope.py', 'code/pin_scope.py']
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
            for p in sorted(bound(paths))}


def load_json(path):
    def finite_float(token):
        value = float(token)
        require(math.isfinite(value), 'Nonfinite JSON number: '+token)
        return value
    def pairs(items):
        out = {}
        for k, v in items:
            require(k not in out, 'Duplicate JSON key: '+k)
            out[k] = v
        return out
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs,
                      parse_float=finite_float,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('Nonfinite JSON: '+x)))


def compare(saved, actual, path='result'):
    if type(actual) is dict:
        require(type(saved) is dict and saved.keys() == actual.keys(), 'Different keys: '+path)
        for k in actual:
            compare(saved[k], actual[k], path+'.'+k)
    elif type(actual) is list:
        require(type(saved) is list and len(saved) == len(actual), 'Different sequence: '+path)
        if path.endswith(('probabilities', 'complete_basis_law')):
            for law in (saved, actual):
                require(all(type(p) in (float, int) and math.isfinite(p) and 0 <= p <= 1
                            for p in law) and abs(sum(law)-1) < 1e-12,
                        'Invalid probability law: '+path)
        for i, (a, b) in enumerate(zip(saved, actual)):
            compare(a, b, path+f'[{i}]')
    elif type(actual) is float:
        require(type(saved) in (float, int) and math.isfinite(saved) and math.isfinite(actual),
                'Invalid numerical value: '+path)
        # Structural zeros are exact; nonzero roundoff is compared separately.
        require(saved == 0 if actual == 0 else abs(saved-actual) <= 2e-11*(1+abs(actual)),
                'Numerical replay differs: '+path)
        if path.endswith(('error', 'residual')) and actual < 1e-12:
            require(0 <= saved < 1e-12, 'Saved residual exceeds its gate: '+path)
    else:
        require(type(saved) is type(actual) and saved == actual, 'Changed value: '+path)


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k,v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def validate_protocol():
    """Reject a frozen declaration that differs from this fixed implementation."""
    protocol = load_json(HERE/'protocol.json')
    require(protocol['schema'] == 'oph.native_base_protocol.v1', 'Unknown protocol schema')
    expected = {
        'scalar': {'mass_coefficient': 1, 'half_link_coefficient': '1/sqrt(2)',
                   'penalty': 1024, 'loop_sign_seed': 193019,
                   'default_ramp_substeps_per_star': 16,
                   'extended_ramp_substeps_per_star': 256,
                   'total_ramp_action_per_star': 16},
        'contrast': {'coupling': 8, 'competing_coupling': 0, 'cell_load': 2,
                     'sample_actions': [0, 0.125, 0.5, 1, 2, 4, 8]},
        'apparatus': {'pulse_strength': 64, 'pulse_action': 'pi/4',
                      'queries': [0, 4, 8], 'source_action_per_count': '1/8',
                      'phases': ['0', 'pi/2', 'pi', '3*pi/2'],
                      'targets': ['single_seam', 'circulation', 'five_seam', 'complex_seed']},
        'architecture': {'code_bits': 432, 'source_modes': 150, 'fuel_modes': 12,
                         'parity_reference_modes': 12, 'additional_local_scratch_bits': 12,
                         'additional_local_pointer_bits': 2,
                         'total_coherent_bits_including_scratch_and_pointers': 446,
                         'controller': 'supplied finite local classical program',
                         'record_policy': 'retain all loop, diagnostic and failure records; no archive reuse',
                         'whole_star_code_gates': 'supplied full unitary extensions, identity on complement',
                         'owner_initialization': 'supplied vacuum on each owner'},
    }
    for name, values in expected.items():
        compare(protocol[name], values, 'protocol.'+name)
    return protocol


def scientific_result(extended=False):
    protocol = validate_protocol()
    from . import replay
    from .routed import run as routed
    from .instruments import run as instruments
    from .gauge_capacity import (compiler_capacity, current_algebra_dimensions,
                                 fock_capacity, phase_response_witness)
    from model.np_repair import run as contrast, primary_join
    d = load_source()
    with np.load(ROOT/'model/np_export/inputs/CELL.npz', allow_pickle=False) as original:
        require(np.array_equal(d['edges'], original['edges']), 'Imported carrier edges differ')
        require(np.array_equal(d['wedges'], original['wedges']), 'Imported carrier wedges differ')
    result = replay()
    result['routed_preparation'] = routed()
    result['instruments'] = instruments(extended=extended)
    result['capacity'] = compiler_capacity()
    result['seven_species_capacity'] = fock_capacity(7)
    result['raw_currents'] = current_algebra_dimensions()
    phase = phase_response_witness()['records']
    result['phase_records'] = {k: {key: [int(p) for p in v] if isinstance(v, tuple) else int(v)
                                 for key,v in r.items()} for k,r in phase.items()}
    result['contrast_alternative'] = contrast()
    result['contrast_primary_join'] = primary_join()
    return _plain(dict(schema='oph.native_base_replay.v1', extended=extended,
                       protocol=protocol, scientific=result, physical_rows_validated=[]))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--record', type=Path)
    modes.add_argument('--verify', type=Path)
    parser.add_argument('--freeze', type=Path, default=HERE/'FREEZE.json')
    args = parser.parse_args(argv)
    freeze = load_json(args.freeze)
    require(bound_pins(freeze['source_sha256']) == source_pins(), 'Production source freeze differs')
    require(freeze['schema'] == 'oph.native_base_freeze.v1', 'Unknown freeze schema')
    validate_protocol()
    if args.record:
        require(not args.record.exists(), 'Refusing to replace recorded evidence')
    result = scientific_result(args.extended)
    result['source_sha256'] = source_pins()
    # A new receipt has the same domain gates as a replayed receipt.
    compare(result, result)
    path = args.verify or HERE/('extended.json' if args.extended else 'results.json')
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        with args.record.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, sort_keys=True, separators=(',', ':'), allow_nan=False)
            f.write('\n')
        print(json.dumps({'recorded': str(args.record), 'extended': args.extended}))
    else:
        saved = load_json(path)
        # Schema, values, types, structural zeros, failed controls and complete
        # outcome arrays are all compared; no unknown fields are ignored.
        compare(saved, result)
        print(json.dumps({'verified': True, 'extended': args.extended,
                          'construction': 'OPH-NSC-1', 'physical_rows_validated': []}))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
