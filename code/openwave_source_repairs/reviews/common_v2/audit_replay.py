"""Temporary-copy read-only and deliberately false named-gate verification."""
from pathlib import Path
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,tempfile,time
HERE=Path(__file__).resolve().parent
def inventory(p):return {str(x.relative_to(p)):hashlib.sha256(x.read_bytes()).hexdigest()for x in p.rglob('*')if x.is_file()}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);a=ap.parse_args()
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1');rows=[]
    with tempfile.TemporaryDirectory(prefix='openwave-common-review-')as d:
        d=Path(d);p=d/'producer';v=d/'review';v.mkdir()
        shutil.copytree(a.source,p,ignore=shutil.ignore_patterns('__pycache__'))
        for n in ('check.py','RESULT.json','PREREG.json'):shutil.copy2(HERE/n,v/n)
        before=inventory(d)
        command=[sys.executable,'-B',str(v/'check.py'),str(p)]
        start=time.perf_counter();x=subprocess.run(command,env=env,capture_output=True,text=True)
        rows.append(dict(test='clean_read_only',exit=x.returncode,seconds=time.perf_counter()-start,stdout=x.stdout,stderr=x.stderr))
        assert x.returncode==0 and before==inventory(d)
        target=p/'RECEIPT.json';original=target.read_bytes();r=json.loads(original)
        for gate in r['gates']:
            changed=json.loads(original);changed['gates'][gate]=False;assert changed['all_gates_pass'] is True
            target.write_text(json.dumps(changed));start=time.perf_counter()
            x=subprocess.run(command,env=env,capture_output=True,text=True)
            rows.append(dict(test='false_named_gate',gate=gate,exit=x.returncode,seconds=time.perf_counter()-start,stderr=x.stderr))
            assert x.returncode!=0 and 'AssertionError' in x.stderr,gate
        target.write_bytes(original);assert before==inventory(d)
    result=dict(all_passed=True,read_only_bytes_unchanged=True,false_named_gate_cases=len(r['gates']),python=sys.version,platform=platform.platform(),rows=rows)
    with(HERE/'REPLAY_AUDIT.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items()if k!='rows'},sort_keys=True))
if __name__=='__main__':main()
