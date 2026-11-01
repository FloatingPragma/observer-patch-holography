"""Cold native-ground replay with complete source records and lossless states."""
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
from model.cube_eta_protection import repair
from . import preparation, circuits, bounds, experiment


ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'
REPLAY_ATOL=1e-12


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/ground_spectral_source/README.md',
        'source_receipts':['model/retarded_matter_response/results.json',
                           'model/interacting_vacuum/results.json'],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/GROUND_SPECTRAL_SOURCE.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


_IMPORTED_BINDINGS=source_bindings()


def compare(expected,actual,path='receipt'):
    if type(expected) is not type(actual):raise ValueError(path+': receipt type changed')
    if type(actual) is dict:
        if expected.keys()!=actual.keys():raise ValueError(path+': receipt keys changed')
        for key in actual:compare(expected[key],actual[key],path+'.'+key)
    elif type(actual) is list:
        if len(expected)!=len(actual):raise ValueError(path+': receipt length changed')
        for j,(a,b) in enumerate(zip(expected,actual,strict=True)):compare(a,b,f'{path}[{j}]')
    elif type(actual) is float:
        if not math.isfinite(expected) or not math.isfinite(actual):raise ValueError(path+': nonfinite evidence')
        exact=any(s in path.rsplit('.',1)[-1] for s in ('bound','admission','tolerance'))
        tolerance=0 if exact else min(REPLAY_ATOL,max(1e-42,abs(actual)*1e-8))
        if abs(expected-actual)>tolerance:raise ValueError(path+': numerical evidence changed')
    elif type(actual) not in (str,int,bool,type(None)) or expected!=actual:
        raise ValueError(path+': receipt value changed')


def validate_accounting(value,arrays):
    if value.get('schema')!='oph.ground_spectral_source.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid ground source receipt')
    result=value['executed_native_readout'];records=result['source_prefix_records']
    if result['ground_preparation_executed'] is not True or result['full17bit_charged_Ramsey_executed'] is not True:
        raise ValueError('ground preparation and charged readout must both execute')
    if [(r['stage'],r['pulses']) for r in records]!=[
        ('ground_tree',5379),('ground_frame',184),('charged_pointer',27),
        ('field_first',337920),('field_second',337920)]:
        raise ValueError('native stage or pulse coverage changed')
    settings=result['settings']
    if len(settings)!=2 or [r['phase_over_pi'] for r in settings]!=['0','1/2']:
        raise ValueError('both complete analyzer settings required')
    all_records=records+[s['analyzer_source_records'] for s in settings]
    keys={r['failure_weight_array'] for r in all_records}
    keys.update(['prepared_ground_success_state','analyzer_0_full_success_state',
                 'analyzer_1_full_success_state'])
    if set(arrays)!=keys:raise ValueError('complete archive key set changed')
    for row in all_records:
        table=arrays[row['failure_weight_array']]
        if (table.shape!=(row['pulses'],4) or table.dtype!=np.dtype('float64') or
                not np.isfinite(table).all() or np.any(table<0)):
            raise ValueError('native failure array shape/type/value changed')
        if row['native_abort_slots']+row['padding_slots']!=4*row['pulses']:
            raise ValueError('native failure/padding coverage changed')
        if abs(float(table.sum())-row['first_failure_probability'])>row['failure_sum_admission']+1e-40:
            raise ValueError('native failure mass changed')
    prepared=arrays['prepared_ground_success_state']
    if (prepared.shape!=(65536,1) or prepared.dtype!=np.dtype('complex128') or
            not np.isfinite(prepared).all()):raise ValueError('complete prepared state changed')
    prep_total=float(np.vdot(prepared,prepared).real)+sum(r['first_failure_probability'] for r in records[:2])
    if abs(prep_total-1)>result['host_numerical_admission']:
        raise ValueError('prepared ground loses unconditioned probability')
    native=result['all130_repair_outcome_slots']
    expected=[(law,[repair.PORTS[j] for j in seam],['native',coin,raw])
              for law in ('C','E') for seam in repair.SEAMS
              for coin in (0,1) for raw in range(4 if law=='C' else 8)]
    actual=[(r['law'],r['seam'],r['record']) for r in native if r['record'][0]=='native']
    outside=[(r['law'],r['seam']) for r in native if r['record']==['outside']]
    if (len(native)!=130 or actual!=expected or outside!=
            [(law,[repair.PORTS[j] for j in seam]) for law in ('C','E') for seam in repair.SEAMS]):
        raise ValueError('complete repair records changed')
    if any(type(r['probability']) is not float or not math.isfinite(r['probability']) or
           not 0<=r['probability']<=1 for r in native):raise ValueError('invalid repair probability')
    native_mass=math.fsum(r['probability'] for r in native if r['record'][0]=='native')
    if abs(native_mass-result['repair_native_first_stop_mass'])>1e-40:
        raise ValueError('native first-stop probability changed')
    prefix=math.fsum(r['first_failure_probability'] for r in records)+native_mass
    for j,setting in enumerate(settings):
        state=arrays[f'analyzer_{j}_full_success_state']
        if (state.shape!=(131072,1) or state.dtype!=np.dtype('complex128') or
                not np.isfinite(state).all()):raise ValueError('complete final state changed')
        probabilities=setting['all_six_pointer_probabilities']
        if len(probabilities)!=6 or any(not math.isfinite(x) or not 0<=x<=1 for x in probabilities):
            raise ValueError('complete pointer probabilities changed')
        calculated=[float(np.vdot(out,out).real) for out in experiment.pointer_outputs(state)]
        if max(abs(a-b) for a,b in zip(calculated,probabilities,strict=True))>REPLAY_ATOL:
            raise ValueError('pointer marginal differs from archived complete state')
        total=prefix+setting['analyzer_source_records']['first_failure_probability']+math.fsum(probabilities)
        if abs(total-setting['complete_probability'])>REPLAY_ATOL or abs(total-1)>result['host_numerical_admission']:
            raise ValueError('complete unconditioned probability changed')
    for key in ('success_or_charge_postselection_used','all78_lines_native_resolved',
                'ground_selected_by_passive_scheduler','physical_species_or_relativistic_mass_identified',
                'host_diagnostics_admit_literal_control_theorem'):
        if result[key] is not False:raise ValueError('native source scope promoted')
    if (value['physical_identification'] is not False or value['source_law_changed'] is not False or
            value['integrated_gate_accepted'] is not False):raise ValueError('unproved integrated claim')


def _produce():
    pins=source_bindings()
    if pins!=_IMPORTED_BINDINGS:raise ValueError('source changed since import; use a fresh interpreter')
    # This authenticates every bound input of the original 78-line source.
    bounds.spectral_rows()
    value={'schema':'oph.ground_spectral_source.v1','all_pass':True,
           'source_preparation':preparation.certificate(),
           'charged_probe':circuits.exact_probe_certificate(),
           'localized_defect':circuits.localized_defect_certificate(),
           'complete_native_error_certificate':bounds.certificate(),
           'spectral_probability_intervals':bounds.probability_intervals(),
           'executed_native_readout':experiment.report(),
           'source_bindings':pins,'replay_float_absolute_tolerance':REPLAY_ATOL,
           'source_law_changed':False,'physical_identification':False,
           'integrated_gate_accepted':False}
    arrays=experiment.arrays()
    if source_bindings()!=pins:raise ValueError('source changed during native ground replay')
    validate_accounting(value,arrays);json.dumps(value,allow_nan=False)
    return value,arrays


def write_pair(target,value,arrays):
    archive=target.with_suffix('.npz')
    if target.exists() or archive.exists():raise FileExistsError('retain existing receipt and archive')
    with archive.open('xb') as stream:np.savez_compressed(stream,**arrays)
    value=dict(value,complete_record_archive={'filename':archive.name,
               'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
               'format':'lossless NumPy NPZ; no object arrays or pickle'})
    with target.open('x',encoding='utf-8',newline='\n') as stream:
        stream.write(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')


def read_pair(target):
    value=strict_json.load(target);metadata=value.pop('complete_record_archive')
    if Path(metadata['filename']).name!=metadata['filename']:raise ValueError('archive must be a local filename')
    archive=target.parent/metadata['filename']
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=metadata['sha256']:
        raise ValueError('archive hash changed')
    with np.load(archive,allow_pickle=False) as stored:arrays={name:stored[name] for name in stored.files}
    validate_accounting(value,arrays)
    return value,arrays


def report():
    pins=source_bindings()
    if pins!=_IMPORTED_BINDINGS:raise ValueError('source changed since import; use a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-ground-source-') as temporary:
        target=Path(temporary)/'worker.json';environment=dict(os.environ)
        environment['OPENBLAS_NUM_THREADS']=environment['OMP_NUM_THREADS']='1'
        result=subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temporary)/'pycache'),
                               '-m','model.ground_spectral_source.packet','--worker',str(target)],
                              cwd=ROOT,env=environment,capture_output=True,text=True)
        if result.returncode:raise RuntimeError('cold native ground replay failed: '+result.stderr[-5000:])
        value,arrays=read_pair(target)
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('source changed across cold native replay')
    return value,arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    actions.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    actions.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args();value,arrays=_produce() if args.worker else report()
    if args.check:
        retained,stored=read_pair(args.check);compare(retained,value)
        if stored.keys()!=arrays.keys():raise ValueError('archive keys changed')
        for name,actual in arrays.items():
            atol,rtol=(1e-42,1e-8) if name.endswith('failure_weights') else (REPLAY_ATOL,0)
            if stored[name].shape!=actual.shape or not np.allclose(stored[name],actual,atol=atol,rtol=rtol):
                raise ValueError('archive replay changed: '+name)
    else:write_pair(args.worker or args.write,value,arrays)
    if not args.worker:print('Verified native original-ground preparation and complete charged Ramsey output; no physical mass assigned.')


if __name__=='__main__':main()
