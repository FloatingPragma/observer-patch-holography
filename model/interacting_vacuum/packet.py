"""Cold exact replay of the original-cube global ground and charged measures."""
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
from . import charged, global_ground, ground


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/interacting_vacuum/README.md',
        'source_receipts': ['model/cube_eta_protection/results.json'],
        'independent_checks': [],
    })
    paths.add('paper/tex_fragments/INTERACTING_VACUUM.tex')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def validate(value):
    if value['global_ground']['global_ground_dimension'] != 5:
        raise ValueError('the entire five-state ground multiplet must be retained')
    if value['global_ground']['all_N_not4_energy_gap_strict_lower'] != '7/1024':
        raise ValueError('charged-sector bound changed')
    if len(value['global_ground']['Gram_blocks']) != 4:
        raise ValueError('all four nontrivial number sectors must be certified')
    if value['charged_spectral_measure']['down_addition']['residual_weight'] != '4/5':
        raise ValueError('unresolved interacting spectral weight was lost')
    if value['same_process_repair']['protection_kind'] != 'passive outside-service identity of the declared full-Fock completion':
        raise ValueError('outside-service scope changed')
    if (value['scope']['physical_species_identifications_added'] != 0 or
            value['scope']['relativistic_masses_constructed'] is not False or
            value['scope']['native_ground_preparation_executed'] is not False):
        raise ValueError('finite ground theorem was promoted beyond its checked scope')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; use a fresh interpreter')
    result = {'schema': 'oph.interacting_vacuum.v1', 'all_pass': True,
              'global_ground': global_ground.certificate(),
              'exact_ground_kernel_and_reference': ground.certificate(),
              'charged_spectral_measure': charged.spectral_certificate(),
              'same_process_repair': charged.repair_certificate(),
              'scope': {'unchanged_original_Weyl_contact_action': True,
                        'all65536_Fock_states_in_ground_comparison': True,
                        'finite_positive_energy_interacting_ground_constructed': True,
                        'ground_uniqueness_claimed': False,
                        'native_ground_preparation_executed': False,
                        'native_ground_detector_executed': False,
                        'all_charged_poles_resolved': False,
                        'thermodynamic_interacting_vacuum_constructed': False,
                        'relativistic_masses_constructed': False,
                        'physical_species_identifications_added': 0,
                        'Einstein_premises_completed': False},
              'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during exact ground replay')
    validate(result)
    return result


def report():
    if source_bindings() != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; use a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-interacting-vacuum-') as temporary:
        path = Path(temporary)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '4'
        worker = subprocess.run([sys.executable, '-W', 'error', '-X',
                                 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.interacting_vacuum.packet', '--worker', str(path)],
                                cwd=ROOT, env=environment, capture_output=True, text=True)
        if worker.returncode:
            raise RuntimeError('fresh exact ground worker failed: '+worker.stderr)
        result = strict_json.loads(path.read_text())
    if source_bindings() != _IMPORTED_BINDINGS or result['source_bindings'] != _IMPORTED_BINDINGS:
        raise ValueError('source changed across cold replay')
    validate(result)
    return result


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    modes.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        args.worker.write_text(canonical(_produce()))
        return
    expected = strict_json.loads(args.check.read_text()) if args.check else None
    result = report()
    if args.write: OUTPUT.write_text(canonical(result))
    if expected is not None and canonical(expected) != canonical(result):
        raise ValueError('exact receipt or source bindings changed')
    print(json.dumps({'all_pass': True, 'global_ground_energy': '-sqrt3/4',
                      'global_ground_dimension': 5, 'positive_charged_gap_lower': '7/1024',
                      'source_bindings': len(result['source_bindings'])}, indent=2))


if __name__ == '__main__': main()
