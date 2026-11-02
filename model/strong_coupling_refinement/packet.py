"""Cold current-source replay, with explicit inherited receipt status."""
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
PRIORS = ('continuum_free_limit', 'centered_source_ground', 'field_symmetry')
FRAGMENT = 'paper/tex_fragments/STRONG_COUPLING_REFINEMENT.tex'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/strong_coupling_refinement/README.md',
        'source_receipts': [f'model/{name}/results.json' for name in PRIORS],
        'independent_checks': [],
    })
    paths.add(FRAGMENT)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_custody():
    result = {}
    for name in PRIORS:
        path = ROOT/f'model/{name}/results.json'
        row = strict_json.loads(path.read_text())
        if not row.get('source_bindings'):
            raise ValueError('prior has no source manifest')
        changed = sorted(p for p, digest in row['source_bindings'].items()
                         if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != digest)
        result[name] = {'historical_all_pass': row.get('all_pass') is True,
                        'current_source_mismatches': changed,
                        'current_source_verified': not changed,
                        'historical_expensive_experiments_reexecuted': False}
    return result


def theorem():
    return {
        'source': 'centered K=T+UD-(U/2)N, real signed t_xy=+/-c/2',
        'finite_graph_separation': 'U>=4 tau, tau=c |E|',
        'full_output_spin_dynamics_bound': 'min(2,4 tau/U+3 |t| tau^3/U^2)',
        'charged_sector_gap_lower': 'U/2-tau',
        'fixed_graph_interacting_spin_limit_proved': True,
        'growing_U_free_class': 'U=o((log c)^(3/2)), c=v/a->infinity',
        'new_unbounded_U_compactness_argument': True,
        'both_charge_branches_retained': True,
        'outside_cyclicity_promotes_to_full_hierarchy': True,
        'free_limit_valleys': 8,
        'free_limit_energy': 'v |p|',
        'free_limit_rest_masses': 'all zero',
        'old_c1_over8_U1_cube_in_new_global_expansion_domain': False,
        'interacting_relativistic_continuum_constructed': False,
        'all_growing_coupling_trajectories_excluded': False,
        'global_bound_failure_is_a_model_no_go': False,
        'physical_particle_species_identified': False,
        'Einstein_dynamics_derived': False,
        'unknown_ground_or_effective_generator_oracle_used': False,
    }


def validate(value):
    if value.get('schema') != 'oph.strong_coupling_refinement.v1' or value.get('new_checks_pass') is not True:
        raise ValueError('invalid strong-coupling receipt')
    # JSON comparison distinguishes false from 0 and exact integers from floats.
    if (canonical(value.get('theorem')) != canonical(theorem()) or
            canonical(value.get('bounds')) != canonical(bounds.certificate())):
        raise ValueError('theorem scope, JSON type, or exact bound changed')
    controls = value['controls']
    if canonical(controls) != canonical(algebra.certificate()):
        raise ValueError('exact source control payload or JSON type changed')
    if controls.get('basis_truncation_before_hopping') is not False:
        raise ValueError('intermediate doublon-hole states were deleted')
    if [r['half_filled_dimension'] for r in controls['exact_graph_audits']] != [6,20,70,12870]:
        raise ValueError('complete finite carrier coverage changed')
    for row in controls['exact_graph_audits']:
        if not all(row[k] is True for k in ('PTP_zero','one_hop_denominator_exactly_U','second_order_exchange_exact')):
            raise ValueError('finite spin reduction premise missing')
    inherited = value['prior_custody']
    if canonical(inherited) != canonical(prior_custody()):
        raise ValueError('prior source coverage, current custody, or JSON type changed')
    for row in inherited.values():
        if row['current_source_verified'] is not (not row['current_source_mismatches']):
            raise ValueError('inherited source mismatch hidden')
        if row['historical_expensive_experiments_reexecuted'] is not False:
            raise ValueError('historical execution falsely attributed')
    passed = all(r['current_source_verified'] and r['historical_all_pass'] for r in inherited.values())
    if value.get('dependency_custody_pass') is not passed or value.get('all_pass') is not passed:
        raise ValueError('integrated receipt status hides a pending prerequisite')


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def _produce():
    pins = source_bindings()
    prior = prior_custody()
    complete = all(r['current_source_verified'] and r['historical_all_pass'] for r in prior.values())
    value = {'schema': 'oph.strong_coupling_refinement.v1',
             'new_checks_pass': True, 'dependency_custody_pass': complete, 'all_pass': complete,
             'theorem': theorem(), 'bounds': bounds.certificate(), 'controls': algebra.certificate(),
             'prior_custody': prior, 'source_bindings': pins,
             'proof_scope': 'Analytic graph and continuum proofs with exact integer/rational controls; inherited expensive runs are not replayed. Current source custody is distinct from historical receipt custody.'}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-strong-coupling-') as temp:
        path = Path(temp)/'worker.json'
        env = dict(os.environ)
        env['OPENBLAS_NUM_THREADS'] = env['OMP_NUM_THREADS'] = '1'
        subprocess.run([sys.executable, '-W', 'error', '-X', 'pycache_prefix='+str(Path(temp)/'cache'),
                        '-m', 'model.strong_coupling_refinement.packet', '--worker', str(path)],
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
            raise ValueError('strong-coupling cold replay changed')
    print(json.dumps({k: value[k] for k in ('new_checks_pass','dependency_custody_pass','all_pass')}))


if __name__ == '__main__':
    main()
