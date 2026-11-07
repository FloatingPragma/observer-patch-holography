"""Replay the source-bound finite clock-identification comparison."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import comparison as c


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_clock_comparison/'+name for name in (
    '__init__.py', 'comparison.py', 'packet.py', 'test_comparison.py', 'PROOF.md', 'README.md'))+(
    'model/native_clock_observables/results.json', 'tools/strict_json.py')


def source_bindings():
    result = c.response.source_bindings()
    result.update(c.observable_packet.source_bindings())
    result.update({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_FILES})
    return dict(sorted(result.items()))


def produce():
    bindings = source_bindings()
    data, probe = c.checked_source()
    detector = c.detector_contract()
    value = {'schema': 'oph.native_clock_comparison.v1',
             'source': c.response.source(),
             'finite_time_comparison': c.finite_time_certificate(data, probe),
             'executed_detector': detector,
             'complete_output_comparison': c.operational_certificate(detector),
             'bare_record_clock_blindness': c.blindness_certificate(data),
             'scope': {
                 'transfer_generator_is_on_the_scalar_record_Hilbert_space': True,
                 'working_density_identified_with_transfer_ground': False,
                 'identity_reference_preserving_clock_identification_rejected': True,
                 'a_transfer_clock_service_is_provided': False,
                 'arbitrary_parameter_native_clock_compilation_proved': False,
                 'clock_specification_readout_is_forbidden': False,
                 'spatial_or_relativistic_clock_identification_proved': False,
                 'different_reference_or_general_encoding_ruled_out': False,
                 'particle_mass_interpretation': False},
             'source_bindings': bindings}
    if source_bindings() != bindings:
        raise ArithmeticError('clock comparison source files changed during replay')
    return value


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if strict_json.load(Path(path)) != produce():
        raise ValueError('clock comparison receipt differs from the owning replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Native clock comparison receipt written.')
    else:
        verify()
        print('Native clock comparison receipt verified.')
