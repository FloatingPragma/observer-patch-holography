"""Cold reproducible complete-spectrum and localization certificate."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import native, spectrum

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def bindings():
    names = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/one_particle_search/PROOF.md',
        'source_receipts': ['model/particle_species_detector/results.json'],
        'independent_checks': ['model/one_particle_search/test_search.py'],
    })
    names.add('model/one_particle_search/README.md')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


_IMPORTED_BINDINGS = bindings()


def build():
    before = bindings()
    if before != _IMPORTED_BINDINGS:
        raise ValueError('source changed after import; replay in a fresh interpreter')
    value = {'schema': 'oph.one_particle_search.v1', 'source_bindings': before,
             'exact_complete_spectrum': spectrum.exact_certificate(),
             'complete_native_detector': native.certificate(),
             'finite_stationary_localization': spectrum.localization_certificate(),
             'infinite_lattice_certificate': spectrum.laurent_certificate(),
             'independent_source_columns': spectrum.source_replay(),
             'independent_native_local_columns': native.local_execution(),
             'scope': {'source_geometry': [8, 8, 8], 'source_speed': '1/8',
                       'source_number_sector': 1, 'all_contact_U_allowed': True,
                       'original_energy_origin_retained': True,
                       'finite_spectral_label_premise_discharged': True,
                       'localized_stationary_N1_core_constructed': False,
                       'observed_species_identifications': 0,
                       'physical_rest_masses_derived': 0,
                       'many_body_defects_excluded': False,
                       'full_interacting_QFT_constructed': False,
                       'entire_completion_gate_accepted': False,
                       'prior_receipt_producers_all_replayed': False}}
    if before != bindings():
        raise ValueError('source changed during complete-spectrum replay')
    return value


def report():
    pins = bindings()
    with tempfile.TemporaryDirectory(prefix='oph-one-particle-') as temporary:
        path = Path(temporary)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '1'
        run = subprocess.run([sys.executable, '-W', 'error', '-X',
                              'pycache_prefix='+str(Path(temporary)/'pycache'),
                              '-m', 'model.one_particle_search.packet', '--worker', str(path)],
                             cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if run.returncode:
            raise RuntimeError('cold one-particle search failed: '+run.stderr[-4000:])
        value = strict_json.load(path)
    if pins != bindings() or value['source_bindings'] != pins:
        raise ValueError('source changed across cold replay')
    return value


def validate(value, replay=None):
    actual = report() if replay is None else replay
    if canonical(value) != canonical(actual):
        raise ValueError('spectrum, localization, detector or scientific scope differs from replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    modes.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    modes.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    value = build() if args.worker else report()
    if args.check:
        validate(strict_json.load(args.check), value)
    else:
        (args.write or args.worker).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    if not args.worker:
        print('PASS complete 1024-mode detector, finite stationary bounds and infinite N1 exclusion')


if __name__ == '__main__':
    main()
