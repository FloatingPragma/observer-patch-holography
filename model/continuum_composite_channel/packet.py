"""Cold source-bound witnesses for the resolved continuum composite channel."""
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
from . import algebra,bounds


ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/continuum_composite_channel/README.md',
        'source_receipts':['model/continuum_field_net/results.json',
                           'model/continuum_ground_constraints/results.json'],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/CONTINUUM_COMPOSITE_CHANNEL.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def prior_proof():
    for name in ('continuum_field_net','continuum_ground_constraints'):
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
            'composite_d_exact_CAR':True,
            'same_ground_mixed_fluctuations_vanish':True,
            'uniform_degree_controlled_Dyson_approximation':True,
            'mixed_defects_vanish_in_all_dynamical_polynomial_contexts':True,
            'joint_equal_time_original_and_composite_CAR':True,
            'valleys_in_each_family':8,
            'microscopic_modes_added':False,
            'original_positive_energy_cyclic_subrepresentation_preserved':True,
            'positive_energy_time_group_on_entire_enlarged_GNS_proved':False,
            'composite_limiting_time_continuity_proved':False,
            'bounded_frequency_memory_survival_proved':False,
            'Hartree_free_closure_proved':False,
            'statistical_independence_or_Gaussianity_proved':False,
            'physical_species_or_masses_identified':False,
            'Einstein_G1_to_G6_completed':False}


def validate(value):
    if value.get('schema')!='oph.continuum_composite_channel.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid composite-channel receipt')
    t=value['theorem']
    for key in ('composite_d_exact_CAR','same_ground_mixed_fluctuations_vanish',
                'uniform_degree_controlled_Dyson_approximation',
                'mixed_defects_vanish_in_all_dynamical_polynomial_contexts',
                'joint_equal_time_original_and_composite_CAR',
                'original_positive_energy_cyclic_subrepresentation_preserved'):
        if t[key] is not True:raise ValueError('composite theorem changed')
    for key in ('microscopic_modes_added','positive_energy_time_group_on_entire_enlarged_GNS_proved',
                'composite_limiting_time_continuity_proved','bounded_frequency_memory_survival_proved',
                'Hartree_free_closure_proved','statistical_independence_or_Gaussianity_proved',
                'physical_species_or_masses_identified','Einstein_G1_to_G6_completed'):
        if t[key] is not False:raise ValueError('composite channel was overpromoted')
    if t['valleys_in_each_family']!=8:raise ValueError('a retained valley was removed')
    if value['bounds']['composite_time_continuity_proved'] is not False:
        raise ValueError('composite continuity was assumed')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    value={'schema':'oph.continuum_composite_channel.v1','all_pass':True,
           'prior_proof':prior_proof(),'theorem':theorem(),
           'onsite_algebra':algebra.certificate(),
           'full_Fock_covariance':algebra.finite_covariance_control(),'bounds':bounds.certificate(),
           'proof_scope':'Exact local CAR and full-Fock covariance witnesses audit the operator identities. The README and fragment prove averaged same-ground fluctuations, uniform degree-controlled dynamics, all-context vanishing and joint equal-time CAR. Composite time continuity and energy survival are not proved.',
           'source_bindings':pins}
    if source_bindings()!=pins:raise ValueError('sources changed during replay')
    validate(value);return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-composite-channel-') as temp:
        path=Path(temp)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        worker=subprocess.run([sys.executable,'-W','error','-X',
                               'pycache_prefix='+str(Path(temp)/'pycache'),'-m',
                               'model.continuum_composite_channel.packet','--worker',str(path)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if worker.returncode:raise RuntimeError('fresh composite-channel worker failed: '+worker.stderr)
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
        if canonical(expected)!=canonical(value):raise ValueError('composite-channel receipt changed')
    elif args.write or args.worker:(args.write or args.worker).write_text(canonical(value))
    if not args.worker:
        print(json.dumps({'all_pass':True,'valleys':8,'resolved_equal_time_composite_channel':True,
                          'bounded_frequency_memory_survival':False,'Lorentz_or_Einstein_completion':False,
                          'source_bindings':len(value['source_bindings'])}))


if __name__=='__main__':main()
