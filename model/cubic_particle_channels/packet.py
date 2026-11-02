"""Cold replay of the actual cubic particle channels and their exact filters."""
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
from . import group, channels

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PRIORS = ('model/centered_particle_spectrum/results.json',
          'model/centered_neutral_spectrum/results.json')
FRAGMENT = 'paper/tex_fragments/CUBIC_PARTICLE_CHANNELS.tex'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/cubic_particle_channels/README.md',
        'source_receipts': list(PRIORS), 'independent_checks': []})
    paths.add(FRAGMENT)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_custody():
    result = {}
    for name in PRIORS:
        row = strict_json.loads((ROOT/name).read_text())
        bindings = row.get('source_bindings', row.get('source_sha256'))
        if not isinstance(bindings, dict) or not bindings:
            raise ValueError('missing prior source bindings')
        changed = sorted(p for p, digest in bindings.items()
                         if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != digest)
        inputs=row['inherited_inputs']
        if 'historical_receipt_source_custody_valid' in inputs:
            inherited_ok=inputs['historical_receipt_source_custody_valid']
            inherited_mismatches=inputs['historical_receipt_source_mismatches']
        else:
            inherited_ok=all(r['historical_source_custody_valid'] for r in inputs['receipts'])
            inherited_mismatches=sorted({p for r in inputs['receipts'] for p in r['mismatches']})
        result[name] = {'current_source_mismatches': changed,
                        'current_source_verified': not changed,
                        'historical_all_pass': row.get('all_pass'),
                        'inherited_dependency_custody_pass': inherited_ok,
                        'inherited_source_mismatches': inherited_mismatches,
                        'historical_expensive_experiments_reexecuted': False}
    return result


def theorem():
    return {
        'source': 'unchanged centered open cube, c=1/8, U=1, N=8 ground',
        'spatial_action': 'rotation about cube center, original Weyl spinor SU2 lift',
        'ground_cubic_irrep': 'A1',
        'binary_group_order': 48,
        'hidden_spin_vector_is_spatial_vector': False,
        'hidden_spin_vector_irreps': {'A2': 1, 'E': 1},
        'each_charged_line': {'G1': 1, 'G2': 1, 'H': 1},
        'first_neutral_triplet': {'A1': 1, 'E': 1},
        'first_density_line': {'T1': 1},
        'classification_scope': 'all prior certified local probe-cyclic lines; exact dark candidates retained',
        'continuous_rotation_symmetry_proved': False,
        'physical_particle_species_identified': False,
        'relativistic_rest_masses_identified': False,
        'new_eigensolve_claimed': False,
        'native_symmetry_filter_execution_claimed': False,
    }


def validate(value):
    if value.get('schema') != 'oph.cubic_particle_channels.v1' or value.get('new_checks_pass') is not True:
        raise ValueError('invalid cubic-channel receipt')
    for key, expected in (('theorem', theorem()), ('group', group.certificate()),
                          ('structure', channels.structure()), ('spectrum', channels.spectrum()), ('prior_custody', prior_custody())):
        if canonical(value.get(key)) != canonical(expected):
            raise ValueError(key+' payload, scope, or JSON type changed')
    inherited = value['prior_custody']
    passed = all(r['current_source_verified'] and r['historical_all_pass'] is True
                 and r['inherited_dependency_custody_pass'] is True for r in inherited.values())
    if value.get('dependency_custody_pass') is not passed or value.get('all_pass') is not passed:
        raise ValueError('inherited dependency status concealed')


def _produce():
    pins = source_bindings()
    prior = prior_custody()
    complete = all(r['current_source_verified'] and r['historical_all_pass'] is True
                   and r['inherited_dependency_custody_pass'] is True for r in prior.values())
    value = {'schema': 'oph.cubic_particle_channels.v1', 'new_checks_pass': True,
             'dependency_custody_pass': complete, 'all_pass': complete,
             'theorem': theorem(), 'group': group.certificate(), 'structure': channels.structure(), 'spectrum': channels.spectrum(),
             'prior_custody': prior, 'source_bindings': pins,
             'proof_scope': 'Exact binary cubic symmetry, commutants and central filters; certified spectral intervals inherited from FA85/88. Source covariance is independently tested. Historical eigensolver certificates are not replayed.'}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-cubic-channel-') as temp:
        path = Path(temp)/'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X', 'pycache_prefix='+str(Path(temp)/'cache'),
                        '-m', 'model.cubic_particle_channels.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.loads(path.read_text())
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold replay')
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
        if canonical(retained) != canonical(value):
            raise ValueError('cubic-channel cold replay changed')
    print(json.dumps({k: value[k] for k in ('new_checks_pass', 'dependency_custody_pass', 'all_pass')}))


if __name__ == '__main__':
    main()
