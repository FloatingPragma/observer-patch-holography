"""Cold replay of the exact original-law Hubbard charge-gap join."""
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
from . import chain, gap, scaling

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/mass_generation_audit/README.md',
        'source_receipts': ['model/field_symmetry/results.json',
                            'model/common_source/mass_response/results.json'],
        'independent_checks': [],
    })
    paths.add('paper/tex_fragments/MASS_GENERATION_AUDIT.tex')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def validate(value):
    if value['schema'] != 'oph.mass_generation_audit.v1' or value['all_pass'] is not True:
        raise ValueError('invalid mass-generation audit')
    original = value['unchanged_source_join']
    if (original['complete_Fock_witness_dimension'] != 256 or original['contact_changed'] is not False or
            original['half_filling_derived_as_mu_zero_global_ground'] is not False):
        raise ValueError('source or reference scope changed')
    charge = value['interaction_generated_charge_gap']
    interval = charge['exact_charge_gap_interval']
    lo, hi = F(interval['lower_exact']), F(interval['upper_exact'])
    if not 0 < lo <= hi < 1 or hi-lo >= F(167, 10**27):
        raise ValueError('exact charge-gap enclosure changed')
    if F(charge['original_H_removal_energy_interval']['upper_exact']) >= 0:
        raise ValueError('canonical reference was confused with original vacuum')
    if value['continuum_boundary']['target_M_determined_by_fixed_source'] is not False:
        raise ValueError('a chosen running threshold cannot be called a derived mass')
    scope = value['scope']
    if (scope['thermodynamic_theorem_imported'] is not True or
            any(scope[name] is not False for name in ('massive_three_dimensional_species_identified',
                                                      'native_ground_preparation_executed',
                                                      'complete_repair_ground_invariance_proved',
                                                      'empirical_mass_fit_used'))):
        raise ValueError('unsupported scope promotion')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; use a fresh interpreter')
    value = {'schema': 'oph.mass_generation_audit.v1', 'all_pass': True,
             'unchanged_source_join': chain.certificate(),
             'interaction_generated_charge_gap': gap.certificate(),
             'continuum_boundary': scaling.certificate(),
             'scope': {'thermodynamic_theorem_imported': True,
                       'exact_original_law_join_proved': True,
                       'positive_canonical_charge_gap_certified': True,
                       'massive_three_dimensional_species_identified': False,
                       'native_ground_preparation_executed': False,
                       'complete_repair_ground_invariance_proved': False,
                       'empirical_mass_fit_used': False},
             'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during charge-gap audit')
    validate(value)
    return value


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; use a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-mass-generation-') as temporary:
        path = Path(temporary)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-W', 'error', '-X',
                                 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.mass_generation_audit.packet', '--worker', str(path)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError('cold mass-generation audit failed: '+result.stderr[-4000:])
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold charge-gap replay')
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
        print('Verified unchanged-law interaction-generated charge gap and fixed-U continuum boundary.')


if __name__ == '__main__':
    main()
