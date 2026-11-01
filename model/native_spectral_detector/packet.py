"""Cold native finite-clock replay with complete states and output records."""
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
from . import preparation, circuits, bounds, experiment, execution


ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'
REPLAY_ATOL=2e-11


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/native_spectral_detector/README.md',
        'source_receipts':['model/spin_resolved_spectrum/results.json',
                           'model/spectral_readout/results.json',
                           'model/ground_spectral_source/results.json'],
        'independent_checks':[]})
    paths.update(['model/native_spectral_detector/contract.cpp',
                  'paper/tex_fragments/NATIVE_SPECTRAL_DETECTOR.tex'])
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


_IMPORTED_BINDINGS=source_bindings()


def authenticate_parents():
    for package in ('spin_resolved_spectrum','spectral_readout','ground_spectral_source'):
        value=strict_json.load(ROOT/f'model/{package}/results.json')
        if value.get('all_pass') is not True:raise ValueError('parent receipt not accepted')
        for name,digest in value['source_bindings'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('frozen parent source changed: '+name)


def compare(expected,actual,path='receipt'):
    if type(expected) is not type(actual):raise ValueError(path+': type changed')
    if type(actual) is dict:
        if expected.keys()!=actual.keys():raise ValueError(path+': keys changed')
        for key in actual:compare(expected[key],actual[key],path+'.'+key)
    elif type(actual) is list:
        if len(expected)!=len(actual):raise ValueError(path+': length changed')
        for j,(a,b) in enumerate(zip(expected,actual,strict=True)):compare(a,b,f'{path}[{j}]')
    elif type(actual) is float:
        if not math.isfinite(expected) or not math.isfinite(actual):raise ValueError(path+': nonfinite')
        exact=any(word in path.rsplit('.',1)[-1] for word in ('bound','admission','tolerance'))
        tolerance=0 if exact else min(REPLAY_ATOL,max(1e-42,abs(actual)*1e-8))
        if abs(expected-actual)>tolerance:raise ValueError(path+': numerical evidence changed')
    elif type(actual) not in (str,int,bool,type(None)) or expected!=actual:
        raise ValueError(path+': value changed')


def validate_accounting(value,arrays):
    if value.get('schema')!='oph.native_spectral_detector.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid native spectral detector receipt')
    result=value['executed_native_detector']
    if result['invariant_ground_preparation_executed'] is not True or result['complete_binary_phase_clock_source_executed'] is not True:
        raise ValueError('both preparation and clock must execute')
    records=result['source_records']
    names=['ground_tree','selector']+[f'conditional_ground_{j}' for j in range(15)]+['charge_instrument','phase_clock']
    if [r['stage'] for r in records]!=names:raise ValueError('source stage coverage changed')
    if [r['pulses'] for r in records]!=[5379,147]+[376]*15+[13,606751]:
        raise ValueError('native source pulse coverage changed')
    keys={r['failure_weight_array'] for r in records}|{'prepared_ground_ensemble','phase_clock_complete_success_state'}
    if set(arrays)!=keys:raise ValueError('archive key coverage changed')
    for record in records:
        table=arrays[record['failure_weight_array']]
        if not math.isfinite(record['first_failure_probability']) or record['first_failure_probability']<0:
            raise ValueError('invalid native failure probability')
        if (table.dtype!=np.float64 or table.shape!=(record['pulses'],4) or
                not np.isfinite(table).all() or np.any(table<0)):
            raise ValueError('native failure table invalid')
        if record['native_abort_slots']+record['padding_slots']!=table.size:
            raise ValueError('native failure/padding slot coverage changed')
        if abs(float(table.sum())-record['first_failure_probability'])>max(1e-40,abs(record['first_failure_probability'])*1e-8):
            raise ValueError('native failure mass changed')
    prepared=arrays['prepared_ground_ensemble'];state=arrays['phase_clock_complete_success_state']
    for array,shape in ((prepared,(65536,15)),(state,(262144,30))):
        if array.dtype!=np.complex128 or array.shape!=shape or not np.isfinite(array).all():
            raise ValueError('complete native quantum array invalid')
    if (len(result['all_selector_pointer_outcomes'])!=1296 or
            len(result['all_joint_selector_charge_clock_outcomes'])!=3240 or
            len(result['all_charge_clock_marginal_outcomes'])!=216):
        raise ValueError('complete measurement outcome coverage changed')
    compare(result['all_selector_pointer_outcomes'],experiment.selector_pointer_probabilities())
    compare(result['all_joint_selector_charge_clock_outcomes'],experiment.joint_pointer_probabilities(state))
    compare(result['all_charge_clock_marginal_outcomes'],experiment.pointer_probabilities(state))
    if not math.isfinite(result['invalid_selector_stop_probability']) or result['invalid_selector_stop_probability']<0:
        raise ValueError('invalid selector stop probability')
    stopped=math.fsum(r['first_failure_probability'] for r in records)+result['invalid_selector_stop_probability']
    if abs(stopped-result['first_failure_and_invalid_selector_stop_mass'])>1e-40:
        raise ValueError('complete first-stop accounting changed')
    total=stopped+math.fsum(r['probability'] for r in result['all_charge_clock_marginal_outcomes'])
    if not math.isfinite(result['complete_probability']) or abs(total-1)>execution.ADMISSION or abs(total-result['complete_probability'])>REPLAY_ATOL:
        raise ValueError('complete detector probability changed')
    prep_stopped=math.fsum(r['first_failure_probability'] for r in records[:-2])+result['invalid_selector_stop_probability']
    if abs(float(np.vdot(prepared,prepared).real)+prep_stopped-1)>execution.ADMISSION:
        raise ValueError('complete prepared mixture probability changed')
    for key in ('host_diagnostics_prove_literal_control_accuracy','particle_number_or_spectral_projection_used',
                'success_charge_or_bin_postselection_used','all78_lines_separately_resolved_by_native_bins',
                'finite_shot_measurement_records_generated','physical_species_or_relativistic_mass_identified'):
        if result[key] is not False:raise ValueError('unproved detector scope promotion')
    if value['source_law_changed'] is not False or value['physical_identification'] is not False or value['integrated_gate_accepted'] is not False:
        raise ValueError('unproved integrated claim')


def _produce():
    pins=source_bindings()
    if pins!=_IMPORTED_BINDINGS:raise ValueError('source changed since import; use a fresh interpreter')
    authenticate_parents()
    value={'schema':'oph.native_spectral_detector.v1','all_pass':True,
           'native_mixed_ground_preparation':preparation.certificate(),
           'complete_native_error_certificate':bounds.certificate(),
           'exact_spectral_probability_intervals':bounds.ideal_probabilities(),
           'finite_shot_bound':bounds.finite_shot_certificate(),
           'executed_native_detector':experiment.report(),
           'source_bindings':pins,'replay_float_absolute_tolerance':REPLAY_ATOL,
           'source_law_changed':False,'physical_identification':False,'integrated_gate_accepted':False}
    arrays=experiment.arrays()
    if source_bindings()!=pins:raise ValueError('source changed during native detector replay')
    validate_accounting(value,arrays);json.dumps(value,allow_nan=False)
    return value,arrays


def write_pair(target,value,arrays):
    archive=target.with_suffix('.npz')
    shards=[target.with_name(target.stem+f'.state-{j}.npz') for j in range(5)]
    if any(p.exists() for p in [target,archive,*shards]):raise FileExistsError('retain existing receipt/archive')
    state_key='phase_clock_complete_success_state'
    records={key:array for key,array in arrays.items() if key!=state_key}
    with archive.open('xb') as stream:np.savez_compressed(stream,**records)
    files=[{'filename':archive.name,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
            'role':'records'}]
    for j,path in enumerate(shards):
        start,stop=6*j,6*j+6
        with path.open('xb') as stream:
            np.savez_compressed(stream,state_columns=arrays[state_key][:,start:stop])
        files.append({'filename':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                      'role':'state_columns','column_start':start,'column_stop':stop})
    data=dict(value,complete_record_archive={'files':files,
             'format':'lossless NumPy NPZ records and five ordered full-row state shards; no pickle or truncation'})
    with target.open('x',encoding='utf-8') as stream:stream.write(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')


def read_pair(target):
    value=strict_json.load(target);metadata=value.pop('complete_record_archive')
    files=metadata['files']
    if len(files)!=6 or files[0].get('role')!='records':raise ValueError('complete six-file archive required')
    arrays={};columns=[]
    for index,entry in enumerate(files):
        if Path(entry['filename']).name!=entry['filename']:raise ValueError('local archive filename required')
        archive=target.parent/entry['filename']
        if hashlib.sha256(archive.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('archive hash changed')
        with np.load(archive,allow_pickle=False) as stored:
            if index==0:arrays={name:stored[name] for name in stored.files}
            else:
                start,stop=6*(index-1),6*index
                if (entry.get('role')!='state_columns' or entry.get('column_start')!=start or
                        entry.get('column_stop')!=stop or stored.files!=['state_columns']):
                    raise ValueError('state shard column coverage changed')
                part=stored['state_columns']
                if part.shape!=(262144,6) or part.dtype!=np.complex128:raise ValueError('full-row state shard invalid')
                columns.append(part)
    arrays['phase_clock_complete_success_state']=np.column_stack(columns)
    validate_accounting(value,arrays)
    return value,arrays


def report():
    pins=source_bindings()
    if pins!=_IMPORTED_BINDINGS:raise ValueError('source changed since import')
    with tempfile.TemporaryDirectory(prefix='oph-native-spectrum-') as temporary:
        target=Path(temporary)/'worker.json';environment=dict(os.environ)
        for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):environment[name]='1'
        result=subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temporary)/'pycache'),
                               '-m','model.native_spectral_detector.packet','--worker',str(target)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if result.returncode:raise RuntimeError('cold native detector replay failed: '+result.stderr[-5000:])
        value,arrays=read_pair(target)
    if source_bindings()!=pins or value['source_bindings']!=pins:raise ValueError('source changed across cold replay')
    return value,arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    action.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    action.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args();value,arrays=_produce() if args.worker else report()
    if args.check:
        expected,stored=read_pair(args.check);compare(expected,value)
        if stored.keys()!=arrays.keys():raise ValueError('archive keys changed')
        for name,actual in arrays.items():
            atol,rtol=(1e-42,1e-8) if name.endswith('failure_weights') else (REPLAY_ATOL,0)
            if stored[name].shape!=actual.shape or not np.allclose(stored[name],actual,atol=atol,rtol=rtol):
                raise ValueError('archive replay changed: '+name)
    else:write_pair(args.worker or args.write,value,arrays)
    if not args.worker:print('Verified complete native invariant-ground preparation and four-bin energy clock; no physical masses assigned.')


if __name__=='__main__':main()
