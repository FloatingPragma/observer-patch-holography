"""Cold replay of the recorded-repair spin limit and its exact controls."""
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
from . import algebra, bounds

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PRIORS = ('model/strong_coupling_refinement/results.json',
          'model/defect_field_bridge/full_domain_results.json')
FRAGMENT = 'paper/tex_fragments/PROTECTED_SPIN_LIMIT.tex'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/protected_spin_limit/README.md',
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
        result[name] = {'current_source_mismatches': changed,
                        'current_source_verified': not changed,
                        'historical_all_pass': row.get('all_pass'),
                        'historical_expensive_experiments_reexecuted': False}
    return result


def theorem():
    return {
        'protected_carrier': 'one fermion per site, two internal spin states',
        'repair_source': 'complete two-flavour full_domain C/E field extension',
        'recorded_repair_identity_on_protected_carrier': True,
        'finite_history_arbitrary_reference_and_record_adaptive_policy': True,
        'joint_volume_coupling_trajectory': 'c=J*E^4,U=J*E^8,J>0,E>=2',
        'growing_volume_interacting_spin_dynamics_constructed': True,
        'spin_lattice_spacing': 'fixed',
        'source_derived_ground_state_exists': True,
        'ground_state_uniqueness_or_algorithm_claimed': False,
        'torus_extension': 'L=2*W^4,W>=3 odd,E=6*W^6,real signed nearest-neighbour hopping',
        'finite_source_first_spacing_upper': '5*J/W^2+J/W^6',
        'first_spacing_counts_ground_multiplicity': True,
        'gap_above_entire_ground_space_bounded': False,
        'infinite_GNS_gapless_particle_proved': False,
        'interacting_relativistic_continuum_constructed': False,
        'physical_particle_species_identified': False,
        'Einstein_dynamics_derived': False,
        'large_carrier_native_pulse_execution_claimed': False,
        'spatial_continuum_or_clock_calibration_assumed': False,
    }


def validate(value):
    if value.get('schema') != 'oph.protected_spin_limit.v1' or value.get('new_checks_pass') is not True:
        raise ValueError('invalid protected-spin receipt')
    for key, expected in (('theorem', theorem()), ('bounds', bounds.certificate()),
                          ('controls', algebra.certificate()), ('prior_custody', prior_custody())):
        if canonical(value.get(key)) != canonical(expected):
            raise ValueError(key+' payload, scope, or JSON type changed')
    inherited = value['prior_custody']
    # The older full-domain record has no all_pass key. Its actual columns are
    # rechecked here; no missing historical status is invented. The strong
    # coupling record explicitly carries the unresolved inherited custody.
    passed = all(r['current_source_verified'] for r in inherited.values()) and (
        inherited[PRIORS[0]]['historical_all_pass'] is True)
    if value.get('dependency_custody_pass') is not passed or value.get('all_pass') is not passed:
        raise ValueError('inherited dependency status concealed')


def _produce():
    pins = source_bindings()
    prior = prior_custody()
    complete = all(r['current_source_verified'] for r in prior.values()) and (
        prior[PRIORS[0]]['historical_all_pass'] is True)
    value = {'schema': 'oph.protected_spin_limit.v1', 'new_checks_pass': True,
             'dependency_custody_pass': complete, 'all_pass': complete,
             'theorem': theorem(), 'bounds': bounds.certificate(), 'controls': algebra.certificate(),
             'prior_custody': prior, 'source_bindings': pins,
             'proof_scope': 'Analytic full-instrument, joint-limit, ground-state and finite-twist proofs; exact integer/rational controls. Numerical full-Fock histories are independent tests. Historical expensive experiments are not replayed.'}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-protected-spin-') as temp:
        path = Path(temp)/'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X', 'pycache_prefix='+str(Path(temp)/'cache'),
                        '-m', 'model.protected_spin_limit.packet', '--worker', str(path)],
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
            raise ValueError('protected-spin cold replay changed')
    print(json.dumps({k: value[k] for k in ('new_checks_pass', 'dependency_custody_pass', 'all_pass')}))


if __name__ == '__main__':
    main()
