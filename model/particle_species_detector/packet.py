"""Cold replay of the complete model-particle detector and physical comparison."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import comparison, detector, native_fft

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/particle_species_detector/README.md',
        'source_receipts': ['model/relativistic_sector_audit/results.json'],
        'independent_checks': [],
    })
    paths.update(('model/particle_species_detector/physical_targets.json',
                  'paper/tex_fragments/PARTICLE_SPECIES_DETECTOR.tex'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def validate(value):
    modes = value['identified_model_modes']
    if len(modes) != 16 or len({row['flavor'] for row in modes}) != 8:
        raise ValueError('all sixteen signed modes and eight flavors required')
    actual = value['complete_ideal_comparison']
    probabilities = np.asarray(actual['output_probabilities_all_1024_rows'])
    if (probabilities.shape != (1024, 16) or not np.isfinite(probabilities).all() or
            np.any(probabilities < 0) or np.max(abs(probabilities.sum(axis=0)-1)) > 2e-11):
        raise ValueError('complete one-particle output law changed')
    addresses = [row['detector_address'] for row in modes]
    if np.max(abs(probabilities[addresses, :]-np.eye(16))) > 2e-11:
        raise ValueError('detector no longer identifies its declared source modes')
    native = value['complete_native_word_certificate']
    if F(native['preparation_plus_detection_complete_source_error']) >= F(831, 10**11):
        raise ValueError('advertised native error exceeded')
    if native['full_Fock_native_state_array_executed'] is not False:
        raise ValueError('the Gaussian comparison cannot be relabeled full native execution')
    samples = value['retained_ideal_model_counts']
    if (len(samples['counts_all_1024_single_particle_addresses']) != 1024 or
            sum(samples['counts_all_1024_single_particle_addresses']) != samples['total_trials']):
        raise ValueError('complete count dictionary changed')
    if value['retrospective_physical_comparison']['new_physical_species_identifications'] != 0:
        raise ValueError('physical identification cannot follow from this detector alone')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; use a fresh interpreter')
    native = native_fft.certificate()
    value = {'schema': 'oph.particle_species_detector.v1', 'all_pass': True,
             'exact_source_identification': detector.exact_species_certificate(),
             'identified_model_modes': detector.mode_catalog(),
             'complete_native_word_certificate': native,
             'formal_native_pulses_per_preparation_and_detection':
                 1+2*native['native_pulses_per_Fourier_transform']+native['occupation_pointer_native_pulses'],
             'formal_source_first_failure_slots':
                 4+2*native['first_export_abort_slots_per_Fourier_transform']+native['occupation_pointer_first_export_abort_slots'],
             'correct_native_label_probability_lower': str(1-F(831, 10**11)),
             'correct_label_probability_scope': 'specified vacuum preparation, exact original-H model-law continuation and complete native detector; every wrong/abort/outside output included',
             'complete_ideal_comparison': detector.simulate_detector(),
             'retained_ideal_model_counts': detector.sample_counts(),
             'finite_periodic_refinement': detector.refinement_table(),
             'retrospective_physical_comparison': comparison.certificate(),
             'scope': {'operational_model_Weyl_flavors_identified': 8,
                       'signed_model_modes_identified': 16,
                       'physical_species_identifications_added': 0,
                       'massive_relativistic_species_constructed': False,
                       'positive_energy_interacting_vacuum_constructed': False,
                       'protected_cube_repair_invariance_transferred': False,
                       'original_field_law_changed': False,
                       'occupation_or_abort_postselection_used': False,
                       'field_time_evolution_compiled_by_this_package': False,
                       'SI_unit_fit_performed': False},
             'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during detector construction')
    validate(value)
    return value


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; use a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-model-particle-detector-') as temporary:
        path = Path(temporary)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-W', 'error', '-X',
                                 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.particle_species_detector.packet', '--worker', str(path)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError('cold model-particle detector failed: '+result.stderr[-4000:])
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold detector replay')
    validate(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    modes.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    modes.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    value = _produce() if args.worker else report()
    if args.check:
        expected = strict_json.load(args.check)
        compare_receipt(expected, value)
        validate(expected)
    elif args.write or args.worker:
        (args.write or args.worker).write_text(canonical(value))
    if not args.worker:
        print('Verified eight operational Weyl flavors, sixteen signed modes, complete native error and scoped physical comparison.')


if __name__ == '__main__':
    main()
