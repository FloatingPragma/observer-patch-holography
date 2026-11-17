#!/usr/bin/env python3
"""Post-outcome integrity audit; mutations are not scientific alternatives."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,tempfile,time
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def snapshot(p):return {str(x.relative_to(p)):sha(x)for x in p.rglob('*')if x.is_file()}
def main():
    before=snapshot(HERE);env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    outputs={};results={}
    with tempfile.TemporaryDirectory(prefix='openwave-common-v2-audit-')as tmp:
        clean=Path(tmp)/'clean';shutil.copytree(HERE,clean)
        start=time.monotonic();r=subprocess.run([sys.executable,str(clean/'headline.py')],cwd=clean,env=env,capture_output=True,text=True)
        outputs['clean_replay.log']=r.stdout+r.stderr
        results['clean_read_only']={'exit_code':r.returncode,'seconds':time.monotonic()-start,'all_bytes_unchanged':snapshot(clean)==before}
        assert r.returncode==0 and snapshot(clean)==before
        start=time.monotonic();r=subprocess.run([sys.executable,str(clean/'integrate.py'),'--record'],cwd=clean,env=env,capture_output=True,text=True)
        outputs['exclusive_record_rejection.log']=r.stdout+r.stderr
        results['exclusive_record']={'exit_code':r.returncode,'seconds':time.monotonic()-start,'all_bytes_unchanged':snapshot(clean)==before}
        assert r.returncode!=0 and 'FileExistsError'in r.stderr and snapshot(clean)==before
        mutant=Path(tmp)/'mutant';shutil.copytree(HERE,mutant)
        source=(mutant/'integrate.py').read_text();old="'finite_joint_return_bound':2*spatial_error<F(11,10000)";new="'finite_joint_return_bound':2*spatial_error<F(1,10**12)"
        assert source.count(old)==1;source=source.replace(old,new);(mutant/'integrate.py').write_text(source)
        reg=json.loads((mutant/'PREREG.json').read_text());reg['files']['integrate.py']=sha(mutant/'integrate.py')
        (mutant/'PREREG.json').write_text(json.dumps(reg,indent=2,sort_keys=True)+'\n')
        for name in ['RECEIPT.json','COMMON_MAPS.npz']:(mutant/name).unlink()
        start=time.monotonic();r=subprocess.run([sys.executable,str(mutant/'integrate.py'),'--record'],cwd=mutant,env=env,capture_output=True,text=True)
        outputs['impossible_threshold.log']=r.stdout+r.stderr
        data=json.loads((mutant/'RECEIPT.json').read_text())
        assert r.returncode!=0 and not data['all_gates_pass']and data['gates']['finite_joint_return_bound']is False
        assert sum(not x for x in data['gates'].values())==1
        outputs['impossible_threshold_receipt.json']=(mutant/'RECEIPT.json').read_text()
        results['real_scientific_threshold_mutation']={'exit_code':r.returncode,'seconds':time.monotonic()-start,
          'old_expression':old,'mutated_expression':new,'mutated_code_sha256':sha(mutant/'integrate.py'),
          'mutated_prereg_sha256':sha(mutant/'PREREG.json'),'exactly_one_failed_gate':'finite_joint_return_bound',
          'mutation_is_integrity_test_not_alternative_law':True}
    assert snapshot(HERE)==before
    result={'post_outcome_verification':True,'original_all_bytes_unchanged':True,'primary_prereg_sha256':sha(HERE/'PREREG.json'),
      'scientific_receipt_sha256':sha(HERE/'RECEIPT.json'),'checks':results}
    out=HERE/'workflow_audit';out.mkdir(exist_ok=False)
    for name,content in outputs.items():(out/name).write_text(content)
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
