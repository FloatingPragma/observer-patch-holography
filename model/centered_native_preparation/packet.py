"""Cold complete native replay with lossless full-row outcome archive."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import target,tree,bounds,experiment,execute

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'
PRIORS=('centered_particle_spectrum','native_spectral_detector')
REPLAY_ATOL=2e-11


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths=integration_source_paths(ROOT,{'domain_and_proof':'model/centered_native_preparation/README.md',
        'source_receipts':[f'model/{n}/results.json' for n in PRIORS],'independent_checks':[]})
    paths.update({'model/native_spectral_detector/contract.cpp',
                  'model/centered_particle_spectrum/witness.npz',
                  'paper/tex_fragments/CENTERED_NATIVE_PREPARATION.tex'})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def custody():
    result={}
    for name in PRIORS:
        row=strict_json.load(ROOT/f'model/{name}/results.json')
        mismatch=[p for p,d in row['source_bindings'].items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=d]
        result[name]={'historical_all_pass':row.get('all_pass') is True,
                      'current_source_mismatches':mismatch,'current_source_verified':not mismatch,
                      'historical_expensive_experiment_reexecuted':False}
    return result


def scope():
    return {'finite_attachment':'FA89','reference':'unique centered N8 hidden-spin0 eta0 ground',
            'operation':'complete native signed-state preparation and M4 charge-tagged energy instrument',
            'clock':'K=H-N/2, referenced to its certified ground energy',
            'actual_native_preparation_executed':True,'actual_native_M4_readout_executed':True,
            'all_first_export_failures_retained':True,'all_field_and_clock_rows_retained':True,
            'ground_state_selected_by_autonomous_repair':False,
            'continuing_bare_field_repair_included':False,
            'all228_lines_resolved_by_native_bins':False,'physical_species_or_masses_identified':False,
            'interacting_continuum_or_Einstein_requirements_completed':False,
            'source_law_changed':False,'integrated_gate_accepted':False}


def compare(expected,actual,path='receipt'):
    if type(expected) is not type(actual):raise ValueError(path+': JSON type changed')
    if type(actual) is dict:
        if expected.keys()!=actual.keys():raise ValueError(path+': keys changed')
        for k in actual:compare(expected[k],actual[k],path+'.'+k)
    elif type(actual) is list:
        if len(expected)!=len(actual):raise ValueError(path+': length changed')
        for j,(x,y) in enumerate(zip(expected,actual,strict=True)):compare(x,y,f'{path}[{j}]')
    elif type(actual) is float:
        if not math.isfinite(actual) or not math.isfinite(expected):raise ValueError(path+': nonfinite')
        tolerance=min(REPLAY_ATOL,max(1e-42,abs(actual)*1e-8))
        if abs(actual-expected)>tolerance:raise ValueError(path+': numerical evidence changed')
    elif type(actual) not in (str,int,bool,type(None)) or actual!=expected:raise ValueError(path+': value changed')


def validate(value,arrays):
    if value.get('schema')!='oph.centered_native_preparation.v1' or value.get('new_checks_pass') is not True:
        raise ValueError('invalid centered native receipt')
    for key,expected in (('scope',scope()),('complete_error',bounds.certificate()),('prior_custody',custody()),
                         ('exact_spectral_probability_intervals',bounds.spectral_probabilities()),
                         ('finite_shot_bound',bounds.old.finite_shot_certificate()),
                         ('full_source_programs',experiment.program_metadata())):
        if canonical(value.get(key))!=canonical(expected):raise ValueError('changed exact scope, bounds or custody: '+key)
    prior=value['prior_custody'];accepted=all(r['historical_all_pass'] and r['current_source_verified'] for r in prior.values())
    if value.get('dependency_custody_pass') is not accepted or value.get('all_pass') is not accepted:
        raise ValueError('inherited pending custody hidden')
    if value.get('source_bindings')!=source_bindings():raise ValueError('current source binding changed')
    shapes={'hidden_success':(65536,1),'prepared_ground':(65536,1),
            'tree_failure_weights':(262137,4),'frame_failure_weights':(184,4),
            'tree_block_offsets':(32765,2),'charge_success':(131072,1),
            'charge_failure_weights':(13,4),'clock_success':(262144,2),
            'clock_failure_weights':(607071,4),'pointer_probabilities':(4,2),
            'literal_Gray_coefficients':(65535,),'checked_conditional_angles':(65535,)}
    if set(arrays)!=set(shapes):raise ValueError('complete archive key coverage changed')
    for key,shape in shapes.items():
        a=arrays[key]
        expected=np.int64 if key=='tree_block_offsets' else np.complex128 if key.endswith('success') or key=='prepared_ground' else np.float64
        if a.shape!=shape or a.dtype!=expected or not np.isfinite(a).all():raise ValueError('full native array type or shape changed: '+key)
        if key.endswith('failure_weights') and np.min(a,initial=0)<0:raise ValueError('negative failure probability')
    for key,actual in [('literal_Gray_coefficients',np.concatenate([c for _,_,c in tree.ground_design()])),
                       ('checked_conditional_angles',np.concatenate([a for _,a,_ in tree.ground_design()]))]:
        if not np.array_equal(arrays[key],actual):raise ValueError('literal angle archive changed')
    for name,events in experiment.programs().items():
        pairs=np.fromiter((e['kind']=='J' for e in events),dtype=bool,count=shapes[name][0])
        if np.any(arrays[name][pairs,2:]!=0):raise ValueError('nonphysical padding is not an outcome')
    out=value['execution'];prep=out['preparation'];clock=out['clock']
    if any(type(actual) is not int or actual!=expected for actual,expected in (
        (prep['native_tree_pulses'],262137),(prep['native_frame_pulses'],184),
        (clock['clock_pulses'],607071),(clock['clock_full_rows'],262144),(clock['classical_charge_columns'],2))):
        raise ValueError('native resource or carrier count changed')
    for key,expected in (('full_centered_ground_preparation_executed',True),('complete_M4_centered_phase_readout_executed',True),
        ('source_clock_centering_explicit',True),('coherent_tag_reused_only_after_classical_readout',True),
        ('all228_lines_separately_resolved',False),('spectral_or_number_projection_used_in_native_run',False),
        ('native_success_postselected_or_renormalized',False),('autonomous_ground_selection_or_bare_repair_protection_proved',False),
        ('physical_particle_species_or_masses_identified',False),('finite_shot_outcome_samples_generated',False)):
        if out.get(key) is not expected:raise ValueError('executed scope changed: '+key)
    for key,reported in [('tree_failure_weights',prep['tree_first_failure_probability']),('frame_failure_weights',prep['frame_first_failure_probability']),
                         ('charge_failure_weights',clock['charge_failure_probability']),('clock_failure_weights',clock['clock_failure_probability'])]:
        if type(reported) is not float or abs(float(arrays[key].sum())-reported)>max(1e-40,abs(reported)*1e-8):raise ValueError('failure mass changed')
    compare(out['all_charge_clock_outcomes'],experiment.pointer_records(arrays['clock_success']))
    actual=np.sum(abs(arrays['clock_success'].reshape(4,65536,2))**2,axis=1)
    if not np.allclose(actual,arrays['pointer_probabilities'],atol=1e-15,rtol=0):raise ValueError('pointer probabilities do not match full state')
    stop=sum(float(arrays[k].sum()) for k in arrays if k.endswith('failure_weights'))
    total=float(np.vdot(arrays['clock_success'],arrays['clock_success']).real)+stop
    if abs(total-1)>execute.ADMISSION or abs(total-out['complete_probability'])>REPLAY_ATOL:raise ValueError('complete native probability changed')
    if abs(stop-out['complete_first_failure_probability'])>max(1e-40,abs(stop)*1e-8):raise ValueError('a stop branch disappeared')
    if len(out['all_charge_clock_outcomes'])!=216:raise ValueError('raw pointer outcome omitted')


def produce():
    pins=source_bindings();execution,arrays=experiment.run()
    prior=custody();accepted=all(r['historical_all_pass'] and r['current_source_verified'] for r in prior.values())
    value={'schema':'oph.centered_native_preparation.v1','new_checks_pass':True,'dependency_custody_pass':accepted,'all_pass':accepted,
           'scope':scope(),'complete_error':bounds.certificate(),'exact_spectral_probability_intervals':bounds.spectral_probabilities(),
           'finite_shot_bound':bounds.old.finite_shot_certificate(),'full_source_programs':experiment.program_metadata(),
           'execution':execution,'prior_custody':prior,
           'source_bindings':pins,'replay_float_absolute_tolerance':REPLAY_ATOL}
    if source_bindings()!=pins:raise ValueError('source changed during complete native run')
    validate(value,arrays);return value,arrays


def write_pair(path,value,arrays):
    archive=path.with_suffix('.npz')
    with archive.open('wb') as stream:np.savez_compressed(stream,**arrays)
    data=dict(value,complete_record_archive={'filename':archive.name,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
              'format':'lossless full-row NumPy archive, no pickle, no number or spectral truncation'})
    path.write_text(canonical(data))


def read_pair(path):
    value=strict_json.load(path);entry=value.pop('complete_record_archive')
    if Path(entry['filename']).name!=entry['filename']:raise ValueError('local archive filename required')
    archive=path.parent/entry['filename']
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('complete archive hash changed')
    with np.load(archive,allow_pickle=False) as data:arrays={k:data[k] for k in data.files}
    validate(value,arrays);return value,arrays


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-centered-native-') as tmp:
        path=Path(tmp)/'worker.json';env=dict(os.environ)
        for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[key]='1'
        subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(tmp)/'cache'),
            '-m','model.centered_native_preparation.packet','--worker',str(path)],cwd=ROOT,env=env,check=True,capture_output=True,text=True)
        value,arrays=read_pair(path)
    if source_bindings()!=pins:raise ValueError('source changed across cold replay')
    return value,arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT);group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--worker',type=Path);args=parser.parse_args()
    value,arrays=produce() if args.worker else report()
    if args.check:
        old,retained=read_pair(args.check);compare(old,value)
        for key,current in arrays.items():
            atol,rtol=(1e-42,1e-8) if key.endswith('failure_weights') else (REPLAY_ATOL,0)
            if not np.allclose(retained[key],current,atol=atol,rtol=rtol):raise ValueError('native array replay changed: '+key)
    else:write_pair(args.worker or args.write,value,arrays)
    if not args.worker:print(canonical({k:value[k] for k in ('new_checks_pass','dependency_custody_pass','all_pass')}),end='')


if __name__=='__main__':main()
