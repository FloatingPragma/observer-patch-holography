"""Cold source-bound replay for the changed clock and its actual ground."""
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
from . import bounds,cube,dimer,faithfulness

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'
PRIORS=('interacting_vacuum','thermodynamic_ground','continuum_field_net',
        'continuum_ground_constraints','spectral_readout','continuum_cpt_constraint')


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/centered_source_ground/README.md',
        'source_receipts':[f'model/{name}/results.json' for name in PRIORS],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/CENTERED_SOURCE_GROUND.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def validate(value):
    if value.get('schema')!='oph.centered_source_ground.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid centered-source receipt')
    row=value['cube']['bounds']
    if row['centered_global_ground_number']!=8 or row['old_N4_ground_remains_centered_ground'] is not False:
        raise ValueError('centered ground sector changed')
    if row['unique_N8_ground_proved_by_imported_Lieb'] is not True:
        raise ValueError('imported finite uniqueness result lost')
    if row['quantified_neutral_gap_proved'] is not False or row['thermodynamic_ground_uniqueness_proved'] is not False:
        raise ValueError('finite uniqueness promoted beyond its domain')
    lieb=value['cube']['imported_Lieb_theorem']
    if lieb['proof_imported_not_reproved'] is not True:
        raise ValueError('external theorem attribution lost')
    if (lieb['real_symmetric_hopping_checked'] is not True
            or lieb['connected_vertices']!=list(range(8))
            or lieb['bipartition_sizes']!=[4,4] or lieb['strict_positive_U']!='1'
            or lieb['particle_number']!=8 or lieb['even_half_filling'] is not True
            or lieb['ground_hidden_spin']!=0 or lieb['ground_multiplicity']!=1):
        raise ValueError('imported spin theorem hypotheses or conclusion changed')
    local=value['cube_local_state']
    if (local['full_four_state_local_CAR_faithful'] is not True
            or local['double_occupation_strict_lower']!='1/15000'
            or local['double_occupation_strict_upper']!='1/8'
            or local['infinite_ground_local_faithfulness_inferred'] is not False
            or local['geometric_modular_charge_identified'] is not False):
        raise ValueError('finite centered local-state theorem changed or promoted')
    for key in ('old_ground_reused_without_check','interacting_Lorentz_continuum_proved','physical_mass_or_species_identified'):
        if value['bounds'][key] is not False:raise ValueError('unproved centered physical conclusion')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    for name in PRIORS:
        prior=strict_json.loads((ROOT/f'model/{name}/results.json').read_text())
        if prior.get('all_pass') is not True:raise ValueError('prior proof did not pass')
        for path,digest in prior['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior proof source changed: '+path)
    value={'schema':'oph.centered_source_ground.v1','all_pass':True,
           'cube':cube.certificate(),'cube_local_state':faithfulness.certificate(),
           'dimer':dimer.certificate(),'bounds':bounds.certificate(),
           'source_bindings':pins,
           'proof_scope':'New charge-clock and new ground. Exact full-source symmetry, imported all-N lower bound, variational centered charge-sector theorem, explicitly imported Lieb uniqueness/spin theorem with exact finite hypotheses, full-CAR centered-cube one-site faithfulness and entropy control, complete dimer readout and centered thermodynamic/continuum construction. Ritz diagnostics do not determine admission.'}
    validate(value)
    if source_bindings()!=pins:raise ValueError('source changed during replay')
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-centered-ground-') as temp:
        path=Path(temp)/'worker.json';env=dict(os.environ)
        env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']='1'
        subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temp)/'pycache'),
                        '-m','model.centered_source_ground.packet','--worker',str(path)],
                       cwd=ROOT,env=env,check=True,capture_output=True,text=True)
        value=strict_json.loads(path.read_text())
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('source changed across cold replay')
    validate(value);return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',nargs='?',const=OUTPUT,type=Path)
    group.add_argument('--check',nargs='?',const=OUTPUT,type=Path)
    group.add_argument('--worker',type=Path)
    group.add_argument('--diagnostic',action='store_true')
    args=parser.parse_args()
    if args.diagnostic:print(canonical(cube.numerical_diagnostic()));return
    if args.worker:args.worker.write_text(canonical(_produce()));return
    value=report()
    if args.write:args.write.write_text(canonical(value))
    else:
        retained=strict_json.loads(args.check.read_text());validate(retained)
        if canonical(value)!=canonical(retained):raise ValueError('centered-source replay changed')
    print(json.dumps({'all_pass':True,'source_bindings':len(value['source_bindings']),
                      'physical_identification':False}))


if __name__=='__main__':main()
