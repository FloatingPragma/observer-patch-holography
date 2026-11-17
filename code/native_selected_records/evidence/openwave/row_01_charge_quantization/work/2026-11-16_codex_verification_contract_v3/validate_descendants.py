"""Bounded subprocess tests of genuine predicates, saved-data comparisons and modes."""
import argparse
import ast
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

HERE=Path(__file__).resolve().parent
INTEGRATION=HERE.parents[2]/'row_01_charge_quantization/reviews/2026-11-16_codex_resubmission_integration'
SOURCE=INTEGRATION/'source_export_v2/source'
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
ENV.pop('PYTHONOPTIMIZE',None)


def snapshot(folder):
    return {str(p.relative_to(folder)):(hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mtime_ns)
            for p in folder.rglob('*') if p.is_file()}


def replace(text, old, new):
    if text.count(old)!=1:raise RuntimeError('Mutation not unique: '+old)
    return text.replace(old,new)


def invoke(command,cwd,label,output,expect,env=ENV):
    begin=time.monotonic()
    try:
        run=subprocess.run(command,cwd=cwd,env=env,capture_output=True,text=True,timeout=60)
        result={'label':label,'command':command,'exit_code':run.returncode,'stdout':run.stdout,
                'stderr':run.stderr,'expected_success':expect,'passed':(run.returncode==0)==expect}
    except subprocess.TimeoutExpired as error:
        result={'label':label,'command':command,'timeout':True,'passed':False,'error':str(error)}
    result['seconds']=time.monotonic()-begin
    (output/(label+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('stdout','stderr')}),flush=True)
    return result


def payloads(output):
    results=[]
    sites={
        'charge':('entries == {"native": 120, "swap_hold": 0}','entries == {"native": 121, "swap_hold": 0}'),
        'energy':('np.max(abs(R @ R - np.eye(33))) == 0','np.max(abs(R @ R - np.eye(33))) > 1'),
        'spin':('visibility < 0 if law == "native" and turns == 1 else visibility > 0','visibility > 0 if law == "native" and turns == 1 else visibility < 0'),
        'common':('impurity < .002 and pair_impurity < .01','impurity < -1 and pair_impurity < -1'),
        'pair':('energy[1]<OLD_PAIR','energy[1]<1000'),
    }
    with tempfile.TemporaryDirectory(prefix='contract-payloads-') as temp:
        work=Path(temp);shutil.copytree(HERE/'descendants/scripts',work/'scripts');shutil.copytree(HERE/'descendants/data',work/'data')
        for name,(old,new) in sites.items():
            path=work/'scripts'/f'm10_4_{name}.py';original=path.read_text();data=work/'data'/f'm10_4_{name}_summary.json'
            before=snapshot(work/'data')
            results.append(invoke([sys.executable,'-B',str(path)],work,name+'_verify',output,True))
            if snapshot(work/'data')!=before:raise RuntimeError('Verification wrote saved evidence')
            # Mutating an existing scientific condition must reject both interpreter modes.
            path.write_text(replace(original,old,new))
            for optimized in (False,True):
                command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(path)]
                results.append(invoke(command,work,name+('_mutant_optimized' if optimized else '_mutant_normal'),output,False))
            path.write_text(original)
            # Changed scientific evidence, not a changed expected verdict label.
            saved=data.read_bytes();value=json.loads(saved)
            def alter(obj):
                if isinstance(obj,dict):
                    for key,item in obj.items():
                        if isinstance(item,float):obj[key]=item+.125;return True
                        if isinstance(item,(dict,list)) and alter(item):return True
                    for key,item in obj.items():
                        if type(item) is int:obj[key]=item+1;return True
                elif isinstance(obj,list):
                    for i,item in enumerate(obj):
                        if isinstance(item,float):obj[i]=item+.125;return True
                        if isinstance(item,(dict,list)) and alter(item):return True
                return False
            if not alter(value):raise RuntimeError('No scientific saved number')
            data.write_text(json.dumps(value));before=snapshot(work/'data')
            results.append(invoke([sys.executable,'-B',str(path)],work,name+'_changed_saved',output,False))
            if snapshot(work/'data')!=before:raise RuntimeError('Rejected replay wrote saved evidence')
            data.write_bytes(saved)
        # The identical first-record helper is embedded in all five scripts.
        charge=work/'scripts/m10_4_charge.py';newrecord=work/'first.json'
        results.append(invoke([sys.executable,'-B',str(charge),'--record-to',str(newrecord)],work,'first_record',output,True))
        before=(newrecord.read_bytes(),newrecord.stat().st_mtime_ns)
        results.append(invoke([sys.executable,'-B',str(charge),'--record-to',str(newrecord)],work,'first_record_overwrite',output,False))
        if before!=(newrecord.read_bytes(),newrecord.stat().st_mtime_ns):raise RuntimeError('First record overwritten')
        results.append(invoke([sys.executable,'-B',str(charge),'--verify',str(work/'missing.json')],work,'missing_saved_record',output,False))
    return results


def producers(output):
    results=[]
    with tempfile.TemporaryDirectory(prefix='contract-producers-') as temp:
        work=Path(temp);shutil.copytree(HERE/'descendants/producers',work/'producers')
        for name,script in [('reservoir','experiment.py'),('pair','pair_source.py')]:
            folder=work/'producers'/name;path=folder/script;before=snapshot(folder)
            results.append(invoke([sys.executable,'-B',str(path)],folder,name+'_verify',output,True))
            if snapshot(folder)!=before:raise RuntimeError('Producer verification wrote source/evidence')
        folder=work/'producers/reservoir';path=folder/'experiment.py'
        path.write_text(replace(path.read_text(),
            "'conditional_transport_prediction': abs(observed-float(expected)) <= threshold}",
            "'conditional_transport_prediction': abs(observed-float(expected)) < -1}"))
        target=work/'failed-scientific-record'
        results.append(invoke([sys.executable,'-B','-O',str(path),'--record-to',str(target)],folder,'failed_gate_record_nonzero',output,False))
        recorded=json.loads((target/'RECEIPT.json').read_text())
        if recorded['passed'] or recorded['checks']['conditional_transport_prediction']:
            raise RuntimeError('Did not exercise actual failed scientific gate')
        before=snapshot(target)
        results.append(invoke([sys.executable,'-B',str(path),'--record-to',str(target)],folder,'producer_overwrite_rejected',output,False))
        if snapshot(target)!=before:raise RuntimeError('Producer overwrote first record')
        folder=work/'producers/pair';path=folder/'pair_source.py'
        path.write_text(replace(path.read_text(),'out==reversed_out and (pp,hh,qq,ll)==(p2,h2,q2,l2)',
                                'out!=reversed_out and (pp,hh,qq,ll)==(p2,h2,q2,l2)'))
        results.append(invoke([sys.executable,'-B','-O',str(path)],folder,'pair_joint_composition_gate',output,False))
    return results


def harness(output):
    results=[]
    # Full 28-job replay, comparison with saved scientific output, and no-write audit.
    receipt=output/'SOURCE_REPLAY.json'
    command=[sys.executable,'-B','-O',str(HERE/'descendants/verify_source.py'),'--source',str(SOURCE),'--record-to',str(receipt)]
    # Orchestration is longer than one bounded job. Each child has its own 60s limit.
    start=time.monotonic();p=subprocess.run(command,env=dict(ENV,PYTHONOPTIMIZE='2'),capture_output=True,text=True,timeout=300)
    result={'label':'full_harness_28','command':command,'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr,
            'passed':p.returncode==0,'seconds':time.monotonic()-start}
    (output/'full_harness_28.json').write_text(json.dumps(result,indent=2)+'\n');results.append(result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('stdout','stderr')}),flush=True)
    # Real parser/check function, tested in separate processes. These are contract
    # probes, not a competing source law. Preserve legitimate historical failures.
    with tempfile.TemporaryDirectory(prefix='contract-parser-') as temp:
        folder=Path(temp);module=str(HERE/'descendants/verify_source.py')
        cases=[('earlier_failure',[{'passed':False},{'passed':True}]),
               ('nested_gate_failure',[{'passed':True,'gates':{'selection':False}}]),
               ('changed_observable',[{'passed':True,'observation':2}]),
               ('historical_failure_allowed',[{'passed':True,'historical_primary_passed':False,'observation':1}])]
        for label,actual in cases:
            expected=[{'passed':True,'observation':1}]
            if label=='historical_failure_allowed':expected=actual
            contract={'expected_output':expected,'required_verdict_paths':[[0,'passed']],'runtime_paths':[]}
            script=folder/'probe.py';script.write_text('import runpy\nf=runpy.run_path('+repr(module)+')["check_output"]\nf('+repr(actual)+','+repr(contract)+')\n')
            results.append(invoke([sys.executable,'-B','-O',str(script)],folder,label,output,label=='historical_failure_allowed'))
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('group',choices=['payloads','producers','harness']);args=parser.parse_args()
    output=HERE/('validation_'+args.group);output.mkdir(exist_ok=False)
    result=globals()[args.group](output)
    summary={'group':args.group,'passed':all(r['passed'] for r in result),
             'results':[{k:v for k,v in r.items() if k not in ('stdout','stderr')} for r in result]}
    (output/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    raise SystemExit(0 if summary['passed'] else 1)
