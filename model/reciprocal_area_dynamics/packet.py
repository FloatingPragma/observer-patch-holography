"""Cold, source-bound replay of the finite reciprocal area dynamics."""
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
from . import bounds, records


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PRIOR = 'model/finite_area_response/results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/reciprocal_area_dynamics/README.md',
        'source_receipts': [PRIOR], 'independent_checks': [],
    })
    paths.add('model/reciprocal_area_dynamics/DERIVATION.md')
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def inherited_custody():
    prior = strict_json.loads((ROOT/PRIOR).read_text())
    if prior.get('schema') != 'oph.finite_area_response.v1' or prior.get('all_pass') is not True:
        raise ValueError('the exact finite area-response source receipt is required')
    for path, digest in prior['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('area-response source custody changed: '+path)
    return {'source_receipt': PRIOR, 'inherited_source_bindings_verified': True,
            'inherited_native_pulse_experiments_reexecuted': False}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def validate(value):
    if value.get('schema') != 'oph.reciprocal_area_dynamics.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid reciprocal area-dynamics receipt')
    result = value['response']
    if (len(result['spectral_enclosures']) != 16
            or len(result['orthogonal_memory_Gram_pivots']) != 14
            or result['minimal_scalar_memory_dimension'] != 14):
        raise ValueError('complete finite spectrum and certified memory directions required')
    for key in ('conditional_thermal_manifold_invariant_under_actual_law',
                'thermal_area_coefficient_is_universal_Newton_normalization',
                'regional_G4_or_Einstein_dynamics_proved'):
        if result[key] is not False:
            raise ValueError('unsupported stationary-manifold or Einstein promotion')
    for law in ('C', 'E'):
        if value['records']['laws'][law]['complete_records'] != 512:
            raise ValueError('all source record branches must be retained')
    if value['records']['pulse_controller_work_and_fresh_bath_preparation_derived'] is not False:
        raise ValueError('the supplied bath and controller services cannot be removed')


def _produce():
    pins = source_bindings()
    value = {'schema': 'oph.reciprocal_area_dynamics.v1', 'all_pass': True,
             'issues': [1144, 1145], 'gates': ['CP-08', 'CP-09'],
             'response': bounds.certificate(), 'records': records.certificate(),
             'inherited_custody': inherited_custody(), 'source_bindings': pins,
             'proof_scope': 'Exact finite area/energy memory, conditional thermal area stationarity and complete-record entropy production under the unchanged reciprocal law. No regional G4, spatial refinement or Einstein dynamics.'}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during reciprocal dynamics replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-reciprocal-area-') as temporary:
        path = Path(temporary)/'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(temporary)/'pycache'),
                        '-m', 'model.reciprocal_area_dynamics.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.loads(path.read_text())
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across the cold reciprocal dynamics replay')
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
            raise ValueError('reciprocal area dynamics receipt differs from replay')
    print(json.dumps({'all_pass': True, 'memory_dimension': 14,
                      'conditional_thermal_area_coefficient': [.578, .579],
                      'Einstein_dynamics': False}))


if __name__ == '__main__':
    main()
