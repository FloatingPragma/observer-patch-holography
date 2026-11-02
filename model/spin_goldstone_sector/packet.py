"""Cold replay of the ordered gapless spin sector and its exact controls."""
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
from . import algebra, bounds, infrared, experiment

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PRIORS = ('model/protected_spin_limit/results.json',)
FRAGMENT = 'paper/tex_fragments/SPIN_GOLDSTONE_SECTOR.tex'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/spin_goldstone_sector/README.md',
        'source_receipts': list(PRIORS), 'independent_checks': []})
    paths.add(FRAGMENT)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_custody():
    result = {}
    for name in PRIORS:
        row = strict_json.loads((ROOT/name).read_text())
        if row.get('new_checks_pass') is not True:
            raise ValueError('missing prior scientific checks')
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
        'carrier': 'one fermion per site, two internal spin states',
        'source_branch': 'c=J*E^4,U=J*E^8,J>0; fixed-spacing cubic AF Heisenberg limit',
        'spatial_boundary_extension': 'even periodic cubes, side at least four',
        'staggered_bragg_weight_strict_lower': '1/1000',
        'selected_magnetization_lower': '1/1000',
        'ordered_ground_sector_constructed': True,
        'pure_ordered_ground_with_positive_gapless_spectrum_exists': True,
        'field_selected_ground_explicit_spectral_window': True,
        'field_selected_ground_itself_proved_pure': False,
        'source_initial_states': 'selected spin states injected through P; no field added to source dynamics',
        'all_exact_source_ground_subsequences_proved_ordered': False,
        'internal_spin_identified_with_spatial_spin': False,
        'large_volume_native_state_preparation_executed': False,
        'sharp_one_particle_pole_proved': False,
        'relativistic_dispersion_or_rest_mass_proved': False,
        'physical_particle_species_identified': False,
        'spatial_continuum_or_Einstein_dynamics_proved': False,
        'primary_imports': {
            'DLS1978': 'finite even-torus Duhamel bound: equations16,21,42,43; Theorems4.1-4.2 and AF Lemma6.1; DOI10.1007/BF01106729',
            'KLS1988': 'energy-sum-rule strategy equations1-4, DOI10.1007/BF01023854; finite IR factor rederived here from DLS',
            'Koma2020': 'arXiv1712.09018v4 Theorems2.1,2.2: maximal magnetization and pure gapless ground',
            'Koma2021': 'arXiv2105.04970v2 equations3.4,4.1,5.3-5.12: susceptibility and filtered thermodynamic limit',
        },
    }


def direct_checks():
    ir=infrared.certificate()
    if ir['integral_below_two_fifths'] is not True:
        raise ArithmeticError('infrared upper bound does not establish order')
    return {'infrared':ir,'algebra':algebra.certificate(),'bounds':bounds.certificate(),
            'order':bounds.order_margins(),'finite_controls':experiment.check()}


def validate(value):
    if value.get('schema') != 'oph.spin_goldstone_sector.v1' or value.get('new_checks_pass') is not True:
        raise ValueError('invalid spin-goldstone receipt')
    for key, expected in (('theorem', theorem()), ('direct_checks', direct_checks()),
                          ('prior_custody', prior_custody())):
        if canonical(value.get(key)) != canonical(expected):
            raise ValueError(key+' payload, scope, or JSON type changed')
    inherited = value['prior_custody']
    # Fresh algebra and integral checks cannot repair inherited receipt custody.
    passed = all(r['current_source_verified'] for r in inherited.values()) and (
        inherited[PRIORS[0]]['historical_all_pass'] is True)
    if value.get('dependency_custody_pass') is not passed or value.get('all_pass') is not passed:
        raise ValueError('inherited dependency status concealed')
    if canonical(value.get('source_bindings')) != canonical(source_bindings()):
        raise ValueError('current source bindings changed')


def _produce():
    pins = source_bindings()
    prior = prior_custody()
    complete = all(r['current_source_verified'] for r in prior.values()) and (
        prior[PRIORS[0]]['historical_all_pass'] is True)
    value = {'schema': 'oph.spin_goldstone_sector.v1', 'new_checks_pass': True,
             'dependency_custody_pass': complete, 'all_pass': complete,
             'theorem': theorem(), 'direct_checks': direct_checks(),
             'prior_custody': prior, 'source_bindings': pins,
             'proof_scope': 'Analytic source-state transfer, finite-to-infinite order and positive-energy spectral-window proofs; exact rational integral and algebra controls. Primary mathematical theorems are imported with explicit hypotheses. Finite diagonalization is a normalization control, not infinite-volume evidence. Historical experiments are not replayed.'}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-spin-goldstone-') as temp:
        path = Path(temp)/'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X', 'pycache_prefix='+str(Path(temp)/'cache'),
                        '-m', 'model.spin_goldstone_sector.packet', '--worker', str(path)],
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
            raise ValueError('spin-goldstone cold replay changed')
    print(json.dumps({k: value[k] for k in ('new_checks_pass', 'dependency_custody_pass', 'all_pass')}))


if __name__ == '__main__':
    main()
