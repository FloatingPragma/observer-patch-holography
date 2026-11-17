"""Clean, read-only replay and deliberately false scientific-gate checks."""
from pathlib import Path
import argparse, hashlib, json, os, platform, shutil, subprocess, sys, tempfile, time
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inventory(p):return {str(x.relative_to(p)):sha(x) for x in p.rglob('*') if x.is_file()}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);args=ap.parse_args()
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONINTMAXSTRDIGITS='0',PYTHONDONTWRITEBYTECODE='1')
    rows=[]
    with tempfile.TemporaryDirectory(prefix='openwave-energy-review-')as d:
        d=Path(d);p=d/'producer';v=d/'review';v.mkdir()
        shutil.copytree(args.source,p,ignore=shutil.ignore_patterns('__pycache__','.git'))
        for n in ('check_math.py','RESULT.json','check_clock.py','CLOCK_RESULT.json'):shutil.copy2(HERE/n,v/n)
        before=inventory(d)
        for checker in ('check_math.py','check_clock.py'):
            start=time.monotonic();x=subprocess.run([sys.executable,'-B',str(v/checker),str(p)],env=env,capture_output=True,text=True)
            rows.append(dict(test='clean_'+checker,exit=x.returncode,seconds=time.monotonic()-start,stdout=x.stdout,stderr=x.stderr))
            assert x.returncode==0,x.stderr
        assert before==inventory(d),'verification changed saved bytes'
        for receipt,checker in [('attempts/v2/RECEIPT.json','check_math.py'),('observer/RECEIPT.json','check_math.py'),('clock/RECEIPT.json','check_clock.py')]:
            target=p/receipt;original=target.read_bytes();record=json.loads(original)
            # Mutate one actual named verdict, retaining a falsely optimistic aggregate.
            # The reviewer must inspect the constituent checks, not only that aggregate.
            for gate in record['gates']:
                changed=json.loads(original);changed['gates'][gate]=False
                assert changed['all_gates_pass'] is True
                target.write_text(json.dumps(changed))
                start=time.monotonic();x=subprocess.run([sys.executable,'-B',str(v/checker),str(p)],env=env,capture_output=True,text=True)
                rows.append(dict(test='false_named_gate',receipt=receipt,gate=gate,exit=x.returncode,seconds=time.monotonic()-start,stderr=x.stderr[-1200:]))
                assert x.returncode!=0 and 'AssertionError' in x.stderr,(receipt,gate)
            target.write_bytes(original)
        assert before==inventory(d),'temporary mutation restoration failed'
    result=dict(python=sys.version,platform=platform.platform(),threads=1,read_only_clean_replays=2,false_named_gates=sum(x['test']=='false_named_gate'for x in rows),all_passed=True,rows=rows)
    with(HERE/'REPLAY_AUDIT.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},sort_keys=True))
if __name__=='__main__':main()
