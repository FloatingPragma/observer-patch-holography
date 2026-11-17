"""Cold source-bound witnesses for the all-valley continuum field-net theorem."""
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
from . import algebra,bounds,cells,valleys


ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/continuum_field_net/README.md',
        'source_receipts':['model/thermodynamic_ground/results.json',
                           'model/retarded_spatial_locality/results.json'],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/CONTINUUM_FIELD_NET.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_proof():
    for name in ('thermodynamic_ground','retarded_spatial_locality'):
        value=strict_json.loads((ROOT/'model'/name/'results.json').read_text())
        if value.get('all_pass') is not True:raise ValueError('prior theorem missing')
        pins=value.get('source_bindings')
        if not pins:raise ValueError('prior source custody missing')
        for path,digest in pins.items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior source changed: '+path)
    return {'all_prior_source_hashes_verified':True,
            'state':'a selected translation-invariant full ground at every spacing',
            'infinite_ground_values_computed':False,
            'prior_expensive_proofs_reexecuted':False}


def theorem():
    return {'refinement':'a->0, c(a)=v/a, fixed v>0, 0<=U(a)<=Umax<infinity',
            'all_eight_valleys_on_one_subsequence':True,
            'all_multitime_polynomial_limits_exist_subsequentially':True,
            'equal_time_full_continuum_CAR':True,
            'continuous_space_and_time_translations':True,
            'positive_energy_on_all_polynomial_GNS_vectors':True,
            'operator_graded_locality_from_polynomial_contexts':True,
            'propagation_cone':'sup-norm distance >96v |time difference|',
            'bounded_smeared_fields':True,
            'uniform_low_momentum_inverse_bound':True,
            'full_CAR_free_Duhamel_error':'|t|[U(a)||f||+v a^2 sum_j||partial_j^3 f||/6]',
            'U_to_zero_free_dynamics_and_unique_filled_sea':True,
            'limit_state_computed':False,
            'unique_interacting_limit_proved':False,
            'fixed_positive_U_interaction_survival_proved':False,
            'Lorentz_covariance_proved':False,
            'physical_species_or_masses_identified':False,
            'Einstein_G1_to_G6_completed':False}


def validate(value):
    if value.get('schema')!='oph.continuum_field_net.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid field-net receipt')
    t=value['theorem']
    for key in ('all_eight_valleys_on_one_subsequence','all_multitime_polynomial_limits_exist_subsequentially',
                'equal_time_full_continuum_CAR','continuous_space_and_time_translations',
                'positive_energy_on_all_polynomial_GNS_vectors','operator_graded_locality_from_polynomial_contexts',
                'bounded_smeared_fields','uniform_low_momentum_inverse_bound',
                'U_to_zero_free_dynamics_and_unique_filled_sea'):
        if t[key] is not True:raise ValueError('field-net theorem changed')
    for key in ('limit_state_computed','unique_interacting_limit_proved',
                'fixed_positive_U_interaction_survival_proved','Lorentz_covariance_proved',
                'physical_species_or_masses_identified','Einstein_G1_to_G6_completed'):
        if t[key] is not False:raise ValueError('continuum field net was overpromoted')
    if t['propagation_cone']!='sup-norm distance >96v |time difference|':
        raise ValueError('propagation speed or geometry changed')
    if value['valleys']['valley_count']!=8 or value['valleys']['nearest_even_translation_phase_one'] is not True:
        raise ValueError('a valley or translation phase was discarded')
    if value['bounds']['inverse_weak_limit_interchanged'] is not False:
        raise ValueError('nonlinear inverse limit was assumed')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    value={'schema':'oph.continuum_field_net.v1','all_pass':True,
           'prior_proof':prior_proof(),'theorem':theorem(),
           'contact_norm':algebra.contact_certificate(),'cell_average':cells.certificate(),
           'valleys':valleys.certificate(),'bounds':bounds.certificate(),
           'proof_scope':'Finite exact witnesses audit the local operator, cell and sign identities. The README and fragment provide the uniform estimates, all-word compactness, GNS, positive-energy and graded-locality proofs; the continuum state is not numerically computed.',
           'source_bindings':pins}
    if source_bindings()!=pins:raise ValueError('sources changed during replay')
    validate(value);return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-continuum-net-') as temp:
        path=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        worker=subprocess.run([sys.executable,'-W','error','-X',
                               'pycache_prefix='+str(Path(temp)/'pycache'),'-m',
                               'model.continuum_field_net.packet','--worker',str(path)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if worker.returncode:raise RuntimeError('fresh continuum-net worker failed: '+worker.stderr)
        value=strict_json.loads(path.read_text())
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('sources changed across cold replay')
    validate(value);return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args();value=_produce() if args.worker else report()
    if args.check:
        expected=strict_json.loads(args.check.read_text());validate(expected)
        if canonical(expected)!=canonical(value):raise ValueError('continuum-net receipt changed')
    elif args.write or args.worker:(args.write or args.worker).write_text(canonical(value))
    if not args.worker:
        print(json.dumps({'all_pass':True,'valleys':8,'positive_energy_continuum_net':True,
                          'fixed_U_interaction_survival':False,'Lorentz_or_Einstein_completion':False,
                          'source_bindings':len(value['source_bindings'])}))


if __name__=='__main__':main()
