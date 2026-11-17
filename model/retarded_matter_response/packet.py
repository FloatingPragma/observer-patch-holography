"""Cold replay of exact source response, memory and certified static values."""
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
from . import response, bounds, static, spectrum


ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/retarded_matter_response/README.md',
        'source_receipts':['model/interacting_vacuum/results.json',
                           'model/thermodynamic_ground/results.json'],
        'independent_checks':[],
    })
    paths.update(['paper/tex_fragments/RETARDED_MATTER_RESPONSE.tex',
                  'model/observer_particles/geometric_matter_join.py',
                  'model/relative_observer_geometry/intrinsic_geometry.py',
                  'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean'])
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_certificates():
    """Authenticate the proved full-ground gap and infinite response interface."""
    for package in ('interacting_vacuum','thermodynamic_ground'):
        prior=strict_json.loads((ROOT/'model'/package/'results.json').read_text())
        if prior.get('all_pass') is not True:raise ValueError('prior proof did not pass')
        for path,digest in prior['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior proof source changed: '+path)
        if package=='interacting_vacuum':
            value=prior['global_ground']
            if (value['global_ground_energy']!='-sqrt(3)/4' or
                    value['all_N_not4_energy_gap_strict_lower']!='7/1024' or
                    value['full_Fock_dimension']!=65536):
                raise ValueError('prior global ground or charged-gap domain changed')
        else:
            value=prior['retarded_response']
            if (value['same_ground_subsequence_as_existence_theorem'] is not True or
                    value['real_axis_limit_or_poles_or_gap_proved'] is not False):
                raise ValueError('thermodynamic response interface changed')
    return {'full_global_ground_gap_authenticated':'7/1024',
            'all_prior_source_hashes_verified':True,
            'prior_full_Gram_solver_reexecuted':False,
            'infinite_response_interface':'same-ground-subsequence convergence uniformly on compact subsets of Imz>0; no finite-cube state convergence rate or inherited gap'}


def validate(value):
    if value.get('schema')!='oph.retarded_matter_response.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid retarded matter response receipt')
    core=value['source_response']
    if (core['same_source_state_dependent_response_proved'] is not True or
            core['nonzero_dynamic_memory_proved'] is not True or
            core['original_H_without_mu_or_scalar_shift'] is not True):
        raise ValueError('exact source response premise missing')
    for key in ('common_metric_or_Einstein_premise_completed',
                'native_ground_preparation_or_response_measurement_executed'):
        if core[key] is not False:raise ValueError('response promoted beyond its proof')
    if core['physical_particle_assignments_added']!=0:raise ValueError('unproved physical assignment')
    if (value['frequency_bounds']['gap_extended_to_excited_reference_or_infinite_volume'] is not False or
            value['static_response']['floating_point_solve_admits_the_result'] is not False):
        raise ValueError('finite domain or exact admission boundary lost')
    resolved=value['resolved_local_spectrum']
    if (resolved['resolved_distinct_positive_lines']!=78 or
            resolved['unresolved_local_probe_weight']!='0' or
            resolved['entire560_state_spectrum_or_degeneracies_resolved'] is not False or
            resolved['native_preparation_detector_or_finite_shot_resolution_executed'] is not False or
            resolved['infinite_volume_poles_or_relativistic_masses_identified'] is not False):
        raise ValueError('local resolved spectrum or its domain changed')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    value={'schema':'oph.retarded_matter_response.v1','all_pass':True,
           'prior_proofs':prior_certificates(),
           'source_response':response.certificate(),
           'frequency_bounds':bounds.certificate(),
           'static_response':static.certificate(),
           'resolved_local_spectrum':spectrum.certificate(),
           'geometry_interface':{
               'existing_constructor':'model/observer_particles/geometric_matter_join.py:source_tetrad',
               'existing_density_metric':'model/relative_observer_geometry/intrinsic_geometry.py:wave_density_metric',
               'available_output':'source-dependent retarded matrix and exact Schur inverse including memory',
               'required_next_join':'controlled local spatial/low-energy principal symbol and full action density shared by matter channels',
               'Hartree_diagonal_is_a_metric':False,
               'fundamental_metric_register_update_required':False,
               'EinsteinRegisterData_inhabitant_from_this_package':False},
           'validation_boundary':'Exact source identities, all78 local spectral lines and rational error bounds are replayed. Numerical linear solves and spectral brackets are proposals; integer residuals, exact annihilation/rank and interval residues admit the results. Prior source hashes authenticate the global gap; its costly Gram solve is not repeated. No infinite response values, metric or native detector history are computed.',
           'source_bindings':pins}
    if source_bindings()!=pins:raise ValueError('source changed during response replay')
    validate(value)
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-retarded-matter-') as temp:
        path=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        worker=subprocess.run([sys.executable,'-W','error','-X',
                               'pycache_prefix='+str(Path(temp)/'pycache'),'-m',
                               'model.retarded_matter_response.packet','--worker',str(path)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if worker.returncode:raise RuntimeError('fresh response worker failed: '+worker.stderr)
        result=strict_json.loads(path.read_text())
    if source_bindings()!=pins or result['source_bindings']!=pins:
        raise ValueError('source changed across fresh response replay')
    validate(result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args()
    expected=None
    if args.check:
        expected=strict_json.loads(args.check.read_text());validate(expected)
    value=_produce() if args.worker else report()
    if args.check:
        if canonical(expected)!=canonical(value):raise ValueError('response receipt or source bindings changed')
    elif args.write or args.worker:(args.write or args.worker).write_text(canonical(value))
    if not args.worker:
        print(json.dumps({'all_pass':True,'source_response_channels':4,
                          'resolved_local_lines':78,
                          'ground_down_memory_Gram':'I/4',
                          'source_bindings':len(value['source_bindings'])}))


if __name__=='__main__':main()
