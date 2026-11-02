"""Cold source-bound replay of complete native local observer instruments."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import numpy as np
from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import algebra,bounds,circuits,execute,experiment,precision

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
PRIORS=('protected_spin_limit','centered_native_preparation')
FRAGMENT='paper/tex_fragments/NATIVE_SPIN_OBSERVER.tex'
REPLAY_ATOL=2e-11


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def source_bindings():
    paths=integration_source_paths(ROOT,{'domain_and_proof':'model/native_spin_observer/README.md',
        'source_receipts':[f'model/{p}/results.json' for p in PRIORS],'independent_checks':[]})
    paths.update({FRAGMENT,'model/native_spectral_detector/contract.cpp'})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def custody():
    result={}
    for name in PRIORS:
        row=strict_json.load(ROOT/f'model/{name}/results.json')
        changed=[p for p,d in row['source_bindings'].items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=d]
        result[name]={'historical_all_pass':row['all_pass'] is True,'current_source_mismatches':changed,
                      'current_source_verified':not changed,'old_expensive_execution_repeated':False}
    return result


def theorem():
    return {'attachment':'FA93','encoding':'one fermion per two-mode site, logical 0=down and 1=up',
            'direct_exchange_native_pulses':62,'complete_C_E_field_repair_executed':True,
            'nonadjacent_CAR_sign_and_reference_test_executed':True,
            'complete_quantum_reference_and_all_outcomes_retained':True,
            'controlled_local_history_trace_bound':'input_error+(2m+1)B_R(T)+2AJ^2T*delta+58N*2^-p',
            'trace_norm_convention':'unhalved Schatten one norm; probability total variation is at most half',
            'arbitrary_accuracy_interface':'exact dyadic integrated port commands from certified source-read coefficient algorithm',
            'fixed_binary64_control_error_sent_to_zero':False,
            'resource_plan_itself_executed':False,'specified_finite_Neel_state_prepared':True,
            'infinite_ground_state_prepared':False,'arbitrary_ground_marginal_oracle_used':False,
            'bare_repair_during_each_constituent_pulse_protected':False,
            'new_physical_noise_or_memory_law_derived':False,'physical_particle_species_identified':False,
            'relativistic_interacting_continuum_or_Einstein_dynamics_derived':False,
            'source_generator_law_changed':False,'integrated_gate_accepted':False}


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
        # Every floating field here is a norm, mass or probability. An absolute
        # replay tolerance must not admit negative values near zero.
        if min(actual,expected)<0:raise ValueError(path+': negative native diagnostic')
        if path.rsplit('.',1)[-1] in ('complete_first_failure_probability','rounded_negative_mass'):
            if abs(math.sqrt(actual)-math.sqrt(expected))>REPLAY_ATOL:
                raise ValueError(path+': native failure amplitude changed')
        elif abs(actual-expected)>REPLAY_ATOL:raise ValueError(path+': native diagnostic changed')
    elif type(actual) not in (str,int,bool,type(None)) or actual!=expected:raise ValueError(path+': value changed')


def compare_arrays(retained,current,key):
    """Compare every output; failure weights use a single total-amplitude budget.

    This is host replay agreement, not a replacement for validate's independently
    reconstructed target, exact source certificate, or probability checks.
    """
    if retained.shape!=current.shape or retained.dtype!=current.dtype:
        raise ValueError('native replay array carrier changed: '+key)
    if not np.isfinite(retained).all() or not np.isfinite(current).all():
        raise ValueError('nonfinite native replay array: '+key)
    if key.endswith('failure_weights'):
        if np.min(retained,initial=0)<0 or np.min(current,initial=0)<0:
            raise ValueError('negative native replay failure weight')
        if np.linalg.norm(np.sqrt(retained)-np.sqrt(current))>REPLAY_ATOL:
            raise ValueError('native replay total failure amplitude changed: '+key)
    elif not np.allclose(retained,current,atol=REPLAY_ATOL,rtol=0):
        raise ValueError('native replay array changed: '+key)


def validate(value,arrays):
    keys={'schema','new_checks_pass','dependency_custody_pass','all_pass','theorem','algebra','precision',
          'bounds','native_error_certificates','execution','prior_custody','source_bindings'}
    if set(value)!=keys or value.get('schema')!='oph.native_spin_observer.v1' or value.get('new_checks_pass') is not True:
        raise ValueError('invalid native spin receipt')
    for k,expected in [('theorem',theorem()),('algebra',algebra.certificate()),('precision',precision.certificate()),
                       ('bounds',bounds.certificate()),('native_error_certificates',experiment.analytic_errors()),
                       ('prior_custody',custody()),('source_bindings',source_bindings())]:
        if canonical(value[k])!=canonical(expected):raise ValueError('exact theorem, type or source changed: '+k)
    accepted=all(r['historical_all_pass'] and r['current_source_verified'] for r in value['prior_custody'].values())
    if value['dependency_custody_pass'] is not accepted or value['all_pass'] is not accepted:
        raise ValueError('inherited pending custody hidden')
    plans=experiment.plans()
    if set(value['execution'])!=set(plans) or set(arrays)!={n+s for n in plans for s in ('_success','_failure_weights')}:
        raise ValueError('complete experiment coverage changed')
    for name,template in plans.items():
        bits,initial,target=experiment.initial_and_target(name)
        success=arrays[name+'_success'];weights=arrays[name+'_failure_weights'];row=value['execution'][name]
        if success.shape!=initial.shape or success.dtype!=np.complex128 or weights.shape!=(len(template),4) or weights.dtype!=np.float64:
            raise ValueError('complete native array carrier changed')
        if not np.isfinite(success).all() or not np.isfinite(weights).all() or np.min(weights,initial=0)<0:
            raise ValueError('invalid native output')
        pair=np.fromiter((k=='J' for k,_,_ in template),bool,len(template))
        if np.any(weights[pair,2:]!=0):raise ValueError('padding became a physical failure outcome')
        count,mass=row.get('rounded_negative_count'),row.get('rounded_negative_mass')
        if type(count) is not int or not 0<=count<=4*len(template) or type(mass) is not float or not 0<=mass<1e-35:
            raise ValueError('invalid host-roundoff diagnostic')
        residual=abs(float(np.vdot(success,success).real)+float(weights.sum())-float(np.vdot(initial,initial).real))
        if residual>execute.ADMISSION:raise ValueError('complete source probability lost')
        run=SimpleNamespace(bits=bits,success=success,failure_weights=weights,
            failure_probability=float(weights.sum()),word_sha256=experiment.analytic_errors()[name]['word_sha256'],
            native_failure_slots=4*len(template)-2*int(pair.sum()),normalization_residual=residual,
            rounded_negative_count=count,rounded_negative_mass=mass)
        compare(experiment.report_case(name,run,target),row,'execution.'+name)


def produce():
    pins=source_bindings();records,arrays=experiment.run();prior=custody()
    accepted=all(r['historical_all_pass'] and r['current_source_verified'] for r in prior.values())
    value={'schema':'oph.native_spin_observer.v1','new_checks_pass':True,'dependency_custody_pass':accepted,'all_pass':accepted,
           'theorem':theorem(),'algebra':algebra.certificate(),'precision':precision.certificate(),
           'bounds':bounds.certificate(),'native_error_certificates':experiment.analytic_errors(),
           'execution':records,'prior_custody':prior,'source_bindings':pins}
    if source_bindings()!=pins:raise ValueError('source changed during native execution')
    validate(value,arrays);return value,arrays


def write_pair(path,value,arrays):
    path=Path(path);archive=path.with_suffix('.npz')
    with archive.open('wb') as stream:np.savez_compressed(stream,**arrays)
    entry={'filename':archive.name,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
           'format':'lossless complete native success/reference and first-failure arrays, no pickle'}
    path.write_text(canonical(dict(value,complete_record_archive=entry)))


def read_pair(path):
    path=Path(path);value=strict_json.load(path);entry=value.pop('complete_record_archive')
    if Path(entry['filename']).name!=entry['filename']:raise ValueError('local archive filename required')
    archive=path.parent/entry['filename']
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('archive hash changed')
    with np.load(archive,allow_pickle=False) as data:arrays={k:data[k] for k in data.files}
    validate(value,arrays);return value,arrays


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-native-spin-') as tmp:
        path=Path(tmp)/'worker.json';env=dict(os.environ)
        for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[key]='1'
        completed=subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(tmp)/'cache'),
                       '-m','model.native_spin_observer.packet','--worker',str(path)],cwd=ROOT,env=env,capture_output=True,text=True)
        if completed.returncode:raise RuntimeError('cold native worker failed:\n'+completed.stderr)
        value,arrays=read_pair(path)
    if source_bindings()!=pins:raise ValueError('source changed across cold replay')
    return value,arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT);group.add_argument('--worker',type=Path)
    args=parser.parse_args();value,arrays=produce() if args.worker else report()
    if args.check:
        old,retained=read_pair(args.check);compare(old,value)
        for key,current in arrays.items():
            compare_arrays(retained[key],current,key)
    else:write_pair(args.worker or args.write,value,arrays)
    if not args.worker:print(canonical({k:value[k] for k in ('new_checks_pass','dependency_custody_pass','all_pass')}),end='')


if __name__=='__main__':main()
