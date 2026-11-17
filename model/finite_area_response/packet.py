"""Fresh exact finite area-response replay with inherited source custody."""
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
from . import response


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PRIOR = 'model/native_equilibration/reciprocal_vacuum_results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/finite_area_response/README.md',
        'source_receipts': [PRIOR], 'independent_checks': [],
    })
    paths.add('model/native_equilibration/reciprocal_vacuum_results.npz')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def inherited_custody():
    prior = strict_json.loads((ROOT/PRIOR).read_text())
    if prior.get('schema') != 'oph.reciprocal_vacuum.v1':
        raise ValueError('the specified reciprocal source receipt is required')
    for path, digest in prior['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('reciprocal source custody changed: '+path)
    artifact = ROOT/'model/native_equilibration'/prior['artifact']['filename']
    if hashlib.sha256(artifact.read_bytes()).hexdigest() != prior['artifact']['sha256']:
        raise ValueError('reciprocal source array custody changed')
    return {'source_receipt': PRIOR, 'inherited_hashes_checked': True,
            'inherited_native_pulse_experiments_reexecuted': False}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def validate(value):
    if value.get('schema') != 'oph.finite_area_response.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid finite area-response certificate')
    result = value['response']
    if len(result['complete_spectra']) != 2 or any(len(x) != 8 for x in result['complete_spectra']):
        raise ValueError('every state in both geometry blocks must be retained')
    for key in ('complete_interacting_QFT_or_Einstein_geometry_established',
                'same_reference_as_Weyl_particle_source'):
        if result[key] is not False:
            raise ValueError('unsupported field/geometry promotion')
    if result['fixed_law_state_response']['unchanged_Hamiltonian_bath_gaps_readout_and_transition'] is not True:
        raise ValueError('the primary response must use the unchanged law')
    if result['equilibrium_Hamiltonian_family_response']['requires_matching_bath_gaps_and_energy_readout_family'] is not True:
        raise ValueError('the equilibrium susceptibility must retain its changed family')


def _produce():
    pins = source_bindings()
    value = {'schema': 'oph.finite_area_response.v1', 'all_pass': True,
             'issues': [1145], 'gates': ['CP-09'],
             'response': response.certificate(), 'inherited_custody': inherited_custody(),
             'source_bindings': pins,
             'proof_scope': 'Exact finite thermal area response under the unchanged reciprocal law, and a separately specified noncommuting equilibrium susceptibility. No common interacting continuum or Einstein derivation.'}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during finite area replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-area-response-') as temp:
        path = Path(temp)/'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(temp)/'pycache'),
                        '-m', 'model.finite_area_response.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.loads(path.read_text())
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across the cold finite area replay')
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
            raise ValueError('finite area-response receipt differs from its replay')
    print(json.dumps({'all_pass': True, 'fixed_law_area_slope': [-0.000628, -0.000627],
                      'physical_Einstein_completion': False}))


if __name__ == '__main__':
    main()
