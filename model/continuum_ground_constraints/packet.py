"""Cold source-bound replay of weak-coupling ground constraints."""
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
from . import bounds

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/continuum_ground_constraints/README.md',
        'source_receipts':['model/thermodynamic_ground/results.json',
                           'model/retarded_spatial_locality/results.json'],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/CONTINUUM_GROUND_CONSTRAINTS.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def validate(value):
    if value.get('schema')!='oph.continuum_ground_constraints.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid ground-constraint receipt')
    result=value['constraints']
    if result['Wick_coefficients']['4']!=8 or result['double_occupancy_interval']!=['0','1/4']:
        raise ValueError('ground constraint changed')
    for key in ('second_moment_equality_passed_to_limit','macroscopic_smeared_field_Gaussianity_inferred',
                'chemical_potential_shift_added','physical_mass_or_species_identified'):
        if result[key] is not False:raise ValueError('unproved limiting conclusion')


def canonical(value):
    return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    for name in ('thermodynamic_ground','retarded_spatial_locality'):
        prior=strict_json.loads((ROOT/f'model/{name}/results.json').read_text())
        if prior.get('all_pass') is not True:raise ValueError('prior proof did not pass')
        for path,digest in prior['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior proof source changed: '+path)
    value={'schema':'oph.continuum_ground_constraints.v1','all_pass':True,
           'constraints':bounds.certificate(),'source_bindings':pins,
           'proof_scope':'Exact rational bounds and countercontrols audit the self-contained variational and CAR proof. No infinite ground values or physical particle assignments are computed.'}
    validate(value)
    if source_bindings()!=pins:raise ValueError('source changed during replay')
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-ground-constraints-') as temp:
        path=Path(temp)/'worker.json';env=dict(os.environ)
        env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']='1'
        subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temp)/'pycache'),
                        '-m','model.continuum_ground_constraints.packet','--worker',str(path)],
                       cwd=ROOT,env=env,check=True,capture_output=True,text=True)
        value=strict_json.loads(path.read_text())
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('source changed across cold replay')
    validate(value)
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',nargs='?',const=OUTPUT,type=Path)
    group.add_argument('--check',nargs='?',const=OUTPUT,type=Path)
    group.add_argument('--worker',type=Path)
    args=parser.parse_args()
    if args.worker:
        args.worker.write_text(canonical(_produce()));return
    value=report()
    if args.write:args.write.write_text(canonical(value))
    else:
        retained=strict_json.loads(args.check.read_text())
        validate(retained)
        if canonical(value)!=canonical(retained):raise ValueError('ground-constraint replay changed')
    print(json.dumps({'all_pass':True,'source_bindings':len(value['source_bindings']),
                      'physical_identification':False}))


if __name__=='__main__':main()
