"""Focused real-process checks on temporary copies; source evidence stays read-only."""
import argparse
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
SOURCE=HERE.parents[1]/'work/2026-11-16_codex_amplified_source_interval_v1'


def require(ok,message):
    if not ok:raise ValueError(message)


def snapshot():
    names=list(json.loads((SOURCE/'PREREG.json').read_text())['files'])+['PREREG.json','HEADLINE.json','CERTIFICATE.json']
    return {name:[hashlib.sha256((SOURCE/name).read_bytes()).hexdigest(),(SOURCE/name).stat().st_mtime_ns,(SOURCE/name).stat().st_ino] for name in names}


def run():
    before=snapshot();results=[]
    for name in ['clean_optimized','failed_impurity_gate_optimized','changed_source_column','changed_saved_interval','refused_overwrite']:
        with tempfile.TemporaryDirectory(prefix='amplified-interval-audit-') as d:
            temp=Path(d)
            for p in SOURCE.iterdir():
                if p.is_file():shutil.copy2(p,temp/p.name)
            command=[sys.executable,'-B','-O','certificate.py'];status=1
            if name=='clean_optimized':status=0;needle='"passed": true'
            elif name=='failed_impurity_gate_optimized':
                p=temp/'evaluator.py';s=p.read_text();needle='native single off-band gate';old="require(single_off < real('.002'), 'native single off-band gate')"
                require(s.count(old)==1,'mutation target');p.write_text(s.replace(old,"require(single_off < real('-1'), 'native single off-band gate')"))
            elif name=='changed_source_column':
                p=temp/'SOURCE_INPUTS.json';x=json.loads(p.read_text());x['sources']['-1']['native_columns'][0]['current']+=1;p.write_text(json.dumps(x));needle='calibrated column decoder'
            elif name=='changed_saved_interval':
                p=temp/'HEADLINE.json';x=json.loads(p.read_text());r=x['cases'][0]['signed_mode_weights'][0]['lo_dyadic'];r[1]=str(int(r[1])+1);p.write_text(json.dumps(x));needle='saved interval evidence differs'
            else:command+=['--record-to','HEADLINE.json'];needle='refusing to overwrite evidence'
            start=time.monotonic();process=subprocess.run(command,cwd=temp,text=True,capture_output=True,timeout=60,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'});seconds=time.monotonic()-start
            record=dict(name=name,command=command,expected_status=status,actual_status=process.returncode,intended_diagnostic=needle,seconds=seconds,stdout=process.stdout,stderr=process.stderr)
            require(process.returncode==status and needle in process.stdout+process.stderr,'unexpected process result '+json.dumps(record))
            results.append(record)
            print(json.dumps({'case':name,'passed':True,'seconds':seconds}),flush=True)
    require(before==snapshot(),'producer bytes/mtime/inodes changed')
    return {'passed':True,'case_count':len(results),'processes':results,'producer_files_bytes_mtime_inodes_unchanged':True}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--record-to',type=Path);args=p.parse_args()
    if args.record_to is not None and args.record_to.exists():raise FileExistsError(args.record_to)
    result=run()
    if args.record_to is not None:
        with args.record_to.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    else:
        saved=json.loads((HERE/'PROCESSES.json').read_text());fields=lambda r:[(x['name'],x['actual_status'],x['intended_diagnostic']) for x in r['processes']]
        require(fields(saved)==fields(result),'saved process verdicts differ')
    print(json.dumps({k:v for k,v in result.items() if k!='processes'},indent=2))
