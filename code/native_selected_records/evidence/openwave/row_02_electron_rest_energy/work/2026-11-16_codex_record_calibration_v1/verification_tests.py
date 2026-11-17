"""Real process/status tests, mutating only temporary copied evidence/code."""
import argparse
from concurrent.futures import ThreadPoolExecutor
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


def require(ok,message):
    if not ok:raise ValueError(message)


def fingerprint():
    paths=list(json.loads((HERE/'PREREG.json').read_text())['files'])+['PREREG.json','experiment_v2.py','VERIFICATION_REPAIR_PREREG.json','records/EVIDENCE.json.gz','records/RECEIPT.json','INDEPENDENT_CHECK.json']
    return {p:(hashlib.sha256((HERE/p).read_bytes()).hexdigest(),(HERE/p).stat().st_mtime_ns,(HERE/p).stat().st_ino) for p in paths}


def case(name):
    with tempfile.TemporaryDirectory(prefix='continuing-calibration-') as directory:
        root=Path(directory)
        for p in HERE.iterdir():
            if p.is_file():shutil.copy2(p,root/p.name)
        shutil.copytree(HERE/'records',root/'records')
        command=[sys.executable,'-B'];expected=0;needle='"passed": true'
        if name.startswith('scientific_predicate_'):
            p=root/'experiment_v2.py';source=p.read_text();old="require(native!=control,'same observer failed to distinguish the probe law')"
            require(source.count(old)==1,'mutant source site')
            p.write_text(source.replace(old,"require(native==control,'same observer failed to distinguish the probe law')"))
            q=root/'VERIFICATION_REPAIR_PREREG.json';data=json.loads(q.read_text());data['experiment_v2_sha256']=hashlib.sha256(p.read_bytes()).hexdigest();q.write_text(json.dumps(data))
            if name.endswith('optimized'):command.append('-O')
            command+=['experiment_v2.py'];expected=1;needle='same observer failed to distinguish the probe law'
        elif name=='changed_receipt':
            p=root/'records/RECEIPT.json';d=json.loads(p.read_text());d['literal_calibration_calls']+=1;p.write_text(json.dumps(d))
            command+=['experiment_v2.py'];expected=1;needle='saved receipt differs'
        elif name=='changed_saved_evidence':
            p=root/'records/EVIDENCE.json.gz';d=json.loads(gzip.decompress(p.read_bytes()));d['trials'][0]['audit']['events'][0][1][3]+=1;p.write_bytes(gzip.compress(json.dumps(d).encode(),mtime=0))
            command+=['experiment_v2.py'];expected=1;needle='saved evidence differs'
        elif name=='exclusive_existing_record':
            command+=['experiment_v2.py','--record-to','records'];expected=1;needle='FileExistsError'
        elif name=='explicit_first_record':command+=['experiment_v2.py','--record-to','new_record']
        elif name=='clean_optimized':command+=['-O','experiment_v2.py']
        elif name=='independent_clean':command+=['-O','independent_check.py']
        elif name.startswith('independent_'):
            p=root/'mutant.json.gz';d=json.loads(gzip.decompress((root/'records/EVIDENCE.json.gz').read_bytes()))
            if name=='independent_drop_trial':d['trials'].pop();needle='exact 5256 trial coverage'
            elif name=='independent_duplicate_trial':d['trials'][-1]=d['trials'][0];needle='exact 5256 trial coverage'
            elif name=='independent_missing_column':d['trials'][0]['observer']['calibration']['columns'].pop();needle='complete observer calibration differs'
            elif name=='independent_first_calibration_record':
                t=d['trials'][0];i=next(i for i,(stage,r) in enumerate(t['audit']['events']) if stage.startswith('calibration_'));t['audit']['events'][i][1][3]+=1;needle='literal native current'
            elif name=='independent_compiler_reversal':
                r=next(iter(d['compiled_catalogs'].values()))['cooling5']['R_folded_permutation'];r[0],r[1]=r[1],r[0];needle='measured reversal/primitive'
            else:raise ValueError(name)
            p.write_bytes(gzip.compress(json.dumps(d).encode(),mtime=0));command+=['-O','independent_check.py','--evidence',str(p)];expected=1
        elif name=='invalid_coin_optimized':
            command+=['-O','-c',"from controller import calibrate;calibrate(lambda e,s:[*e,2,0],[0,1,2],1,1)"]
            expected=1;needle='address/coin record'
        elif name=='duplicate_column_optimized':
            command+=['-O','-c',"import json;from compiler import table;r=json.load(open('records/RECEIPT.json'));c=r['representative_measured_tables']['native:1']['columns'];c[-1]=c[0];table(c)"]
            expected=1;needle='duplicate column'
        else:raise ValueError(name)
        started=time.monotonic();run=subprocess.run(command,cwd=root,text=True,capture_output=True,timeout=60,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'});elapsed=time.monotonic()-started
        actual=run.stdout+run.stderr
        require(run.returncode==expected and needle in actual,'wrong process outcome '+name+' '+actual[-1000:])
        if name=='explicit_first_record':
            data=json.loads(gzip.decompress((root/'new_record/EVIDENCE.json.gz').read_bytes()));old=json.loads(gzip.decompress((root/'records/EVIDENCE.json.gz').read_bytes()));data.pop('seconds');old.pop('seconds');require(data==old,'first-record full result')
        return dict(name=name,command=[str(Path(c).name) if c==sys.executable else c for c in command],expected_exit=expected,actual_exit=run.returncode,intended_diagnostic=needle,seconds=elapsed,stdout=run.stdout,stderr=run.stderr,passed=True)


def run():
    before=fingerprint();names=['clean_optimized','independent_clean','scientific_predicate_normal','scientific_predicate_optimized','changed_receipt','changed_saved_evidence','exclusive_existing_record','explicit_first_record','independent_drop_trial','independent_duplicate_trial','independent_missing_column','independent_first_calibration_record','independent_compiler_reversal','invalid_coin_optimized','duplicate_column_optimized']
    started=time.monotonic()
    with ThreadPoolExecutor(max_workers=2) as pool:
        records=list(pool.map(case,names))
    require(fingerprint()==before,'sealed inputs/evidence changed during tests')
    return dict(passed=True,tests=records,case_count=len(records),inputs_evidence_bytes_mtime_inode_unchanged=True,seconds=time.monotonic()-started)


def main():
    p=argparse.ArgumentParser();p.add_argument('--record-to',type=Path);args=p.parse_args()
    if args.record_to is not None and args.record_to.exists():raise FileExistsError(args.record_to)
    result=run()
    if args.record_to is not None:
        with args.record_to.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    else:
        saved=json.loads((HERE/'VERIFICATION_TESTS.json').read_text())
        require([(r['name'],r['actual_exit'],r['passed']) for r in saved['tests']]==[(r['name'],r['actual_exit'],r['passed']) for r in result['tests']],'saved verification verdicts')
    print(json.dumps({k:v for k,v in result.items() if k!='tests'},indent=2))


if __name__=='__main__':main()
