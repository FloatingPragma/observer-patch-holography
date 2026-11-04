"""Cold replay of complete-source finite modular-inclusion certificates."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..plaquette_reference.preparation import error_theorem
from . import inclusion


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PRIOR = 'model/certified_reference/results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/finite_modular_inclusion/README.md',
        'source_receipts': [PRIOR], 'independent_checks': [],
    })
    paths.add('Lean/QFT/FiniteModularInclusion.lean')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prepared_reference_bound(prior):
    if prior.get('schema') != 'oph.certified_reference.v1' or prior.get('all_pass') is not True:
        raise ValueError('the literal source requires its proved complete-output certificate')
    for path, digest in prior['source_sha256'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('literal source certificate custody changed: '+path)
    error = F(prior['preparation_certificate']['complete_preparation_output_error'])
    if not 0 <= error < F(17, 10**13):
        raise ValueError('the retained full-source preparation bound changed')
    design = error_theorem()
    density_error = (2+error)*error
    minimum = design['full_minimum_weight']-density_error
    if minimum <= 0:
        raise ValueError('the actual source success reference is not certified faithful')
    delta = (design['density_trace_error']+density_error)/minimum
    if delta != F(prior['modular_certificate']['full_modular_generator_error']):
        raise ValueError('the full-source modular bound disagrees with its independent recomputation')
    lower = inclusion.leakage_bound(generator_error=delta)
    if not lower > F(1, 25000):
        raise ValueError('the literal source error erased the finite obstruction')
    return {'inherited_source_certificate': PRIOR,
            'inherited_source_hashes_checked': True,
            'inherited_native_control_proof_reexecuted': False,
            'full_success_reference_minimum_weight': str(minimum),
            'complete_preparation_error': str(error),
            'full_modular_generator_error': str(delta),
            'distance_from_every_proper_regional_algebra_lower': str(lower),
            'common_simplified_distance_lower': '1/25000',
            'success_block_is_subnormalized_with_all_abort_outputs_retained': True}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def validate(value):
    if value.get('schema') != 'oph.finite_modular_inclusion.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid finite modular-inclusion certificate')
    result = value['finite_source']
    expected = {tuple(s) for n in (1, 2, 3)
                for s in inclusion.combinations(range(4), n)}
    rows = result['all_proper_site_regions']
    if len(rows) != 14 or {tuple(row['sites']) for row in rows} != expected:
        raise ValueError('all fourteen proper site regions are required')
    if result['full_Fock_dimension'] != 256:
        raise ValueError('the complete field carrier is required')
    for key in ('geometric_null_translation_constructed',
                'continuum_half_sided_inclusion_excluded',
                'Einstein_or_physical_particle_gate_accepted'):
        if result[key] is not False:
            raise ValueError('unsupported geometric or physical promotion')


def _produce():
    pins = source_bindings()
    value = {'schema': 'oph.finite_modular_inclusion.v1', 'all_pass': True,
             'issues': [1145], 'gates': ['CP-09'],
             'finite_source': inclusion.certificate(),
             'literal_reference': prepared_reference_bound(strict_json.loads((ROOT/PRIOR).read_text())),
             'source_bindings': pins,
             'proof_scope': 'Finite-dimensional inverse invariance and exact complete-CAR boundary leakage. No continuum half-sided or Einstein exclusion.'}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during the exact regional replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-finite-inclusion-') as temp:
        path = Path(temp)/'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(temp)/'pycache'),
                        '-m', 'model.finite_modular_inclusion.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.loads(path.read_text())
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across the cold regional replay')
    validate(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--worker', type=Path)
    args = parser.parse_args()
    if args.worker:
        args.worker.write_text(canonical(_produce()))
        return
    value = report()
    if args.write:
        args.write.write_text(canonical(value))
    else:
        retained = strict_json.loads(args.check.read_text())
        validate(retained)
        if canonical(value) != canonical(retained):
            raise ValueError('retained finite inclusion result differs from its replay')
    print(json.dumps({'all_pass': True, 'proper_regions': 14,
                      'literal_reference_distance_lower': '1/25000',
                      'continuum_exclusion': False}))


if __name__ == '__main__':
    main()
