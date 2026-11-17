"""Cold replay of the controlled one-dimensional relativistic charge sector."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import diagnostics, dispersion, source

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/hubbard_relativistic_sector/README.md',
        'source_receipts': ['model/mass_generation_audit/results.json'],
        'independent_checks': [],
    })
    paths.add('paper/tex_fragments/HUBBARD_RELATIVISTIC_SECTOR.tex')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def validate(value):
    if value['schema'] != 'oph.hubbard_relativistic_sector.v1' or value['all_pass'] is not True:
        raise ValueError('invalid relativistic-sector receipt')
    if value['source_duality']['complete_Fock_dimension'] != 256:
        raise ValueError('complete source witness changed')
    data = value['controlled_relativistic_dispersion']
    if len(data['finite_parameter_dispersion']) != 15 or len(data['uniform_bounds']) != 6:
        raise ValueError('finite-parameter controls missing')
    bound = data['uniform_bounds'][-1]['absolute_mass_shell_residual_over_M_squared']
    if F(bound['upper_exact']) >= F(17,1000):
        raise ValueError('advertised uniform mass-shell bound exceeded')
    if data['anchored_trajectory']['trajectory_selected_by_unchanged_source'] is not False:
        raise ValueError('chosen trajectory relabeled source-selected')
    mass = data['anchored_trajectory']['chosen_trajectory_rest_mass']
    if not F(1375,100) < F(mass['lower_exact']) <= F(mass['upper_exact']) < F(1376,100):
        raise ValueError('anchored model mass changed')
    scopes = value['scope']
    if (scopes['spatial_dimension'] != 1 or scopes['centered_canonical_reference_required'] is not True or
            any(scopes[name] is not False for name in ('running_family_derived_from_fixed_bare_point',
                     'native_ground_preparation_executed', 'native_local_CAR_pole_proved',
                     'three_dimensional_massive_species_identified', 'physical_SI_mass_calculated'))):
        raise ValueError('unsupported relativistic-sector scope promotion')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; use a fresh interpreter')
    value = {'schema': 'oph.hubbard_relativistic_sector.v1', 'all_pass': True,
             'source_duality': source.certificate(),
             'controlled_relativistic_dispersion': dispersion.certificate(),
             'separate_fixed_source_diagnostics': diagnostics.fixed_source(),
             'scope': {'spatial_dimension': 1, 'centered_canonical_reference_required': True,
                       'thermodynamic_exact_dispersion_imported': True,
                       'uniform_relativistic_mass_shell_limit_proved': True,
                       'running_family_derived_from_fixed_bare_point': False,
                       'native_ground_preparation_executed': False,
                       'native_local_CAR_pole_proved': False,
                       'three_dimensional_massive_species_identified': False,
                       'physical_SI_mass_calculated': False},
             'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during relativistic-sector construction')
    validate(value)
    return value


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; use a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-hubbard-relativistic-') as temporary:
        path = Path(temporary)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-W', 'error', '-X',
                                 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.hubbard_relativistic_sector.packet', '--worker', str(path)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError('cold relativistic-sector replay failed: '+result.stderr[-4000:])
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold relativistic-sector replay')
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
        print('Verified controlled one-dimensional relativistic dispersion and chosen-trajectory rest mass.')


if __name__ == '__main__':
    main()
