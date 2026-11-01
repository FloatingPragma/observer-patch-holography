"""Fresh-process exact regional-ground certificate and source custody."""
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
from . import regional

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/ground_regional_modular/README.md',
        'source_receipts':['model/spin_resolved_spectrum/results.json'],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/GROUND_REGIONAL_MODULAR.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def validate(value):
    if value.get('schema')!='oph.ground_regional_modular.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid ground regional modular receipt')
    result=value['regional_ground']
    if len(result['all_site_pairs'])!=28 or len(result['single_sites'])!=8:
        raise ValueError('incomplete regional source census')
    for row in result['all_site_pairs']:
        if row['full_dimension']!=16 or row['support_rank']!=8 or row['kernel_rank']!=8:
            raise ValueError('the regional zero eigenspace must be retained')
    for key in ('faithful_on_full_regional_CAR','full_CAR_first_law_inferred',
                'geometric_cap_or_Einstein_identified','physical_particles_or_masses_identified'):
        if result[key] is not False:raise ValueError('unsupported ground modular promotion')


def canonical(value):
    return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    prior=strict_json.loads((ROOT/'model/spin_resolved_spectrum/results.json').read_text())
    if prior.get('all_pass') is not True:raise ValueError('source ground proof did not pass')
    for path,digest in prior['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('source ground proof changed: '+path)
    value={'schema':'oph.ground_regional_modular.v1','all_pass':True,
           'regional_ground':regional.certificate(),'source_bindings':pins,
           'proof_scope':'Exact signed full-output regional traces and a self-contained Slater/spin-twirl modular proof on the actual finite particle ground. No geometric or physical-spectrum identification.'}
    validate(value)
    if source_bindings()!=pins:raise ValueError('source changed during regional replay')
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-ground-regional-') as temp:
        path=Path(temp)/'worker.json';env=dict(os.environ)
        env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']='1'
        subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temp)/'pycache'),
                        '-m','model.ground_regional_modular.packet','--worker',str(path)],
                       cwd=ROOT,env=env,check=True,capture_output=True,text=True)
        value=strict_json.loads(path.read_text())
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('source changed across cold regional replay')
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
        if canonical(value)!=canonical(retained):raise ValueError('regional ground replay changed')
    print(json.dumps({'all_pass':True,'source_bindings':len(value['source_bindings']),
                      'full_site_pairs':28,'physical_identification':False}))


if __name__=='__main__':main()
