"""Replay exact finite premises of the thermodynamic-ground existence proof."""
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
from . import filled_sea, local, volume, retarded, energy

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/thermodynamic_ground/README.md',
        'source_receipts':['model/field_symmetry/results.json',
                           'model/interacting_vacuum/results.json'],
        'independent_checks':[],
    })
    paths.update(['model/full_fock_field/README.md',
                  'model/full_fock_field/LOCAL_SIMULATION.md',
                  'paper/tex_fragments/THERMODYNAMIC_GROUND.tex'])
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


_IMPORTED_BINDINGS=source_bindings()


def validate(value):
    if value.get('schema')!='oph.thermodynamic_ground.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid thermodynamic ground receipt')
    theorem=value['existence_theorem']
    for key in ('stationary_under_original_thermodynamic_dynamics',
                'gauge_invariant_state_exists','translation_invariant_ground_state_exists_by_Folner_average',
                'interacting_derivation_represented_faithfully'):
        if theorem[key] is not True:raise ValueError('proved ground premise missing')
    for key in ('Hamiltonian_scalar_shift_or_chemical_potential_added',
                'ground_state_computed_or_prepared','chosen_subsequence_or_convergence_rate',
                'state_functional_local_faithfulness_asserted','uniqueness_or_purity_asserted',
                'uniform_gap_or_clustering_asserted','complete_repair_stationarity_asserted',
                'spacing_zero_or_Lorentz_continuum_asserted',
                'nontrivial_scattering_or_vacuum_connected_four_point_asserted'):
        if theorem[key] is not False:raise ValueError('ground existence promoted beyond its proof')
    if (value['density']['number_density_lower']!='1/4' or
            value['density']['number_density_upper']!='1' or
            value['density']['particle_hole']['half_filling_concluded'] is not False):
        raise ValueError('original-mu density scope changed')
    if (value['local_source']['complete_one_site_CAR_basis_checked']!=16 or
            value['free_filled_sea_control']['interacting_energy_variance']!='2/3' or
            theorem['physical_particle_identifications_added']!=0):
        raise ValueError('complete local source or negative control changed')
    if (value['retarded_response']['same_ground_subsequence_as_existence_theorem'] is not True or
            value['retarded_response']['real_axis_limit_or_poles_or_gap_proved'] is not False or
            value['retarded_response']['infinite_response_values_computed'] is not False):
        raise ValueError('upper-half-plane response was promoted to real poles or execution')
    if (value['energy_density']['energy_density_exists'] is not True or
            value['energy_density']['state_convergence_rate_or_particle_rest_mass'] is not False or
            value['energy_density']['finite_ground_solver_run_for_growing_boxes'] is not False):
        raise ValueError('energy-density bracket scope changed')


def _produce():
    pins=source_bindings()
    if pins!=_IMPORTED_BINDINGS:raise ValueError('source changed since import')
    result={'schema':'oph.thermodynamic_ground.v1','all_pass':True,
            'existence_theorem':local.existence_certificate(),
            'local_source':local.source_certificate(),
            'density':volume.density_certificate(),
            'free_filled_sea_control':filled_sea.certificate(),
            'retarded_response':retarded.certificate(),
            'energy_density':energy.certificate(),
            'validation_boundary':'Exact finite algebra, rational boundary bounds and complete proof are replayed. No infinite ground density matrix or subsequence is computed. Existing receipt hashes bind prior proofs; the expensive earlier global-cube Gram replay is not rerun here.',
            'source_bindings':pins}
    if source_bindings()!=pins:raise ValueError('source changed during replay')
    validate(result)
    return result


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def report():
    if source_bindings()!=_IMPORTED_BINDINGS:raise ValueError('source changed since import')
    with tempfile.TemporaryDirectory(prefix='oph-thermodynamic-ground-') as temp:
        path=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        worker=subprocess.run([sys.executable,'-W','error','-X',
                               'pycache_prefix='+str(Path(temp)/'pycache'),'-m',
                               'model.thermodynamic_ground.packet','--worker',str(path)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if worker.returncode:raise RuntimeError('fresh ground worker failed: '+worker.stderr)
        result=strict_json.loads(path.read_text())
    if source_bindings()!=_IMPORTED_BINDINGS or result['source_bindings']!=_IMPORTED_BINDINGS:
        raise ValueError('source changed across fresh replay')
    validate(result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args()
    value=_produce() if args.worker else report()
    if args.check:
        expected=strict_json.loads(args.check.read_text());validate(expected)
        if canonical(expected)!=canonical(value):raise ValueError('receipt or source bindings changed')
    elif args.write or args.worker:(args.write or args.worker).write_text(canonical(value))
    if not args.worker:
        print(json.dumps({'all_pass':True,'ground_existence':'fixed-spacing original-mu0 interacting3D lattice',
                          'density_interval':['1/4','1'],'source_bindings':len(value['source_bindings'])}))


if __name__=='__main__':main()
