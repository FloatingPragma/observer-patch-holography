"""Cold source-bound replay of the conditional continuum CPT constraint."""
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

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'
PRIORS=('continuum_field_net','continuum_ground_constraints','spectral_readout')
BOUNDARIES=('Lorentz_covariance_alone_excluded',
            'centered_dynamics_ground_positivity_proved',
            'actual_nonzero_coupling_continuum_measure_computed',
            'physical_mass_or_species_identified',
            'whole_action_or_all_field_dictionaries_excluded')


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/continuum_cpt_constraint/README.md',
        'source_receipts':[f'model/{name}/results.json' for name in PRIORS],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/CONTINUUM_CPT_CONSTRAINT.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def validate(value):
    if value.get('schema')!='oph.continuum_cpt_constraint.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid CPT-constraint receipt')
    result=value['constraints']
    if result['component_count']!=16 or result['same_field_same_clock_CPT_vacuum_excluded_when_U0_positive'] is not True:
        raise ValueError('conditional CPT constraint changed')
    if any(result[key] is not False for key in BOUNDARIES):
        raise ValueError('unproved physical or universal conclusion')
    if value['controls']['centering_failure']['same_reference_is_centered_ground'] is not False:
        raise ValueError('centering positivity counterexample lost')


def canonical(value):
    return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    for name in PRIORS:
        prior=strict_json.loads((ROOT/f'model/{name}/results.json').read_text())
        if prior.get('all_pass') is not True:raise ValueError('prior proof did not pass')
        for path,digest in prior['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior proof source changed: '+path)
    value={'schema':'oph.continuum_cpt_constraint.v1','all_pass':True,
           'constraints':bounds.certificate(),'controls':algebra.certificate(),
           'source_bindings':pins,
           'proof_scope':'Exact rational moment bounds and full-CAR/Clifford controls audit the conditional same-field, same-ground, same-clock CPT obstruction. No Lorentz-alone or whole-action no-go is asserted.'}
    validate(value)
    if source_bindings()!=pins:raise ValueError('source changed during replay')
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-cpt-constraint-') as temp:
        path=Path(temp)/'worker.json';env=dict(os.environ)
        env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']='1'
        subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temp)/'pycache'),
                        '-m','model.continuum_cpt_constraint.packet','--worker',str(path)],
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
        retained=strict_json.loads(args.check.read_text());validate(retained)
        if canonical(value)!=canonical(retained):raise ValueError('CPT-constraint replay changed')
    print(json.dumps({'all_pass':True,'source_bindings':len(value['source_bindings']),
                      'physical_identification':False}))


if __name__=='__main__':main()
