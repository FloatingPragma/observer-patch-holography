"""Targeted verification-only rejection tests with Python assertions disabled."""
from pathlib import Path
import copy,gzip,hashlib,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
def run():
    original={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'records').iterdir()}
    d=json.loads(gzip.decompress((HERE/'records/EVIDENCE.json.gz').read_bytes()));results=[]
    with tempfile.TemporaryDirectory(prefix='admission-optimized-') as tmp:
        for name in ('wrong_current','substituted_density'):
            x=copy.deepcopy(d)
            if name=='wrong_current':x['rows'][0]['audit']['events'][0][1][3]=999
            else:next(iter(x['complete_prefixes'].values()))['clock_density']['real'][0][0]=.5
            p=Path(tmp)/(name+'.json.gz');p.write_bytes(gzip.compress(json.dumps(x).encode(),mtime=0));before=hashlib.sha256(p.read_bytes()).hexdigest()
            command=[sys.executable,'-B','-O',str(HERE/'check.py'),'--evidence',str(p)]
            r=subprocess.run(command,capture_output=True,text=True);out=json.loads(r.stdout)
            if r.returncode==0 or out.get('passed') is not False:raise ValueError('optimized mutant accepted')
            if hashlib.sha256(p.read_bytes()).hexdigest()!=before:raise ValueError('mutant rewritten')
            results.append({'name':name,'exit_code':r.returncode,'output':out})
    after={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'records').iterdir()}
    if after!=original:raise ValueError('original records changed')
    return {'passed':True,'verification_only':True,'scientific_checks_without_assertions':results,'saved_records_unchanged':True}
if __name__=='__main__':
    try:print(json.dumps(run(),indent=2))
    except Exception as e:print(json.dumps({'passed':False,'error':repr(e)}));raise SystemExit(1)
