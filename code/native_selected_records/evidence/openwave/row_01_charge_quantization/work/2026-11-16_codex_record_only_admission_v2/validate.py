"""Verification-only clean copy and scientific mutants; stdout report, temporary writes."""
import copy,gzip,hashlib,json,random,shutil,subprocess,sys,tempfile,time
from pathlib import Path

HERE=Path(__file__).resolve().parent
def require(x,message):
    if not x:raise ValueError(message)
def inventory(p):return {str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in p.rglob('*') if f.is_file() and '__pycache__' not in f.parts}
def invoke(p,script,extra=(),optimized=False):
    cmd=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(p/script),*map(str,extra)]
    start=time.perf_counter();r=subprocess.run(cmd,cwd=p,capture_output=True,text=True)
    try:out=json.loads(r.stdout)
    except ValueError:out={'stdout':r.stdout,'stderr':r.stderr}
    return {'command':cmd,'exit_code':r.returncode,'seconds':time.perf_counter()-start,'result':out,'stderr':r.stderr}

def run():
    started=time.perf_counter();original=inventory(HERE)
    with tempfile.TemporaryDirectory(prefix='record-admission-') as tmp:
        root=Path(tmp);p=root/'packet';shutil.copytree(HERE,p,ignore=shutil.ignore_patterns('__pycache__'))
        before=inventory(p);clean=[]
        for script in ('experiment.py','check.py'):
            r=invoke(p,script);require(r['exit_code']==0 and r['result']['passed'],'clean replay failed');clean.append(r)
        require(inventory(p)==before,'clean verification wrote saved files')
        d=json.loads(gzip.decompress((p/'records/EVIDENCE.json.gz').read_bytes()));mutants=[]
        def test(name,change,expected):
            x=copy.deepcopy(d);change(x);f=root/(name+'.json.gz')
            f.write_bytes(gzip.compress(json.dumps(x).encode(),mtime=0));h=hashlib.sha256(f.read_bytes()).hexdigest()
            r=invoke(p,'check.py',['--evidence',f]);require((r['exit_code']==0)==expected,name)
            require(hashlib.sha256(f.read_bytes()).hexdigest()==h,'mutant was rewritten')
            require(inventory(p)==before,'verification changed original copy')
            mutants.append({'name':name,'expected_pass':expected,**r})
        test('false_scientific_verdict',lambda x:x.update(passed=False),False)
        test('wrong_current',lambda x:x['rows'][0]['audit']['events'][0][1].__setitem__(3,999),False)
        test('invented_reservoir_load',lambda x:x['rows'][0]['audit']['final'].__setitem__(2,999),False)
        test('wrong_extracted_sign',lambda x:x['rows'][0]['controller'].__setitem__('sign',1),False)
        test('missing_counter_outcome',lambda x:next(iter(x['complete_prefixes'].values()))['records'][-1]['all_probabilities'].pop(),False)
        test('substituted_clock_density',lambda x:next(iter(x['complete_prefixes'].values()))['clock_density']['real'][0].__setitem__(0,.5),False)
        def bless_control(x):
            row=next(r for r in x['rows'] if r['control'] and r['controller']['outcome']=='repair_failed');row['controller']['outcome']='signed_target_unit'
        test('false_control_success',bless_control,False)
        def blind(x):
            random.Random(2026).shuffle(x['rows'])
            for i,r in enumerate(x['rows']):r.update(seed_audit_only=-999-i,total_audit_only=10000+i,filename=f'wrong-total-{i}.json',run_order=-i)
        test('blind_seed_total_filename_order',blind,True)
        receipt_path=p/'records/RECEIPT.json';saved=receipt_path.read_bytes();receipt=json.loads(saved);receipt['passed']=False
        receipt_path.write_text(json.dumps(receipt));changed=inventory(p)
        failures=[]
        for script in ('experiment.py','check.py'):
            r=invoke(p,script);require(r['exit_code']!=0,'false receipt accepted');require(inventory(p)==changed,'false receipt repaired silently');failures.append(r)
        receipt_path.write_bytes(saved);require(inventory(p)==before,'restoration mismatch')
        r=invoke(p,'experiment.py',optimized=True);require(r['exit_code']!=0,'optimized legacy checks bypassed');failures.append(r)
        r=invoke(p,'experiment.py',['--record-to',p/'records']);require(r['exit_code']!=0,'existing first record overwritten');failures.append(r)
        require(inventory(p)==before,'failed mode changed records')
    require(inventory(HERE)==original,'original packet changed')
    return {'passed':True,'clean_copy':clean,'scientific_mutants':mutants,'strict_failure_exits':failures,
            'source_and_records_unchanged':True,'seconds':time.perf_counter()-started,'verification_only':True,'external_G6':False}

if __name__=='__main__':
    try:print(json.dumps(run(),indent=2))
    except Exception as e:print(json.dumps({'passed':False,'error':repr(e)}));raise SystemExit(1)
