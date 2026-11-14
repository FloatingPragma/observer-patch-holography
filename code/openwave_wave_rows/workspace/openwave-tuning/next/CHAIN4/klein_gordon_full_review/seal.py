"""Replay the independent review with original-workspace reads forbidden."""
from pathlib import Path
import hashlib,json,os,shutil,site,subprocess,sys,sysconfig,tempfile

HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[3]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    inputs=json.loads((HERE/'INPUTS.json').read_text())['files_sha256']
    for relative,h in inputs.items():assert sha(HERE/'inputs'/relative)==h
    with tempfile.TemporaryDirectory(prefix='oph-kg-full-review-')as temporary:
        target=Path(temporary)/'review'
        shutil.copytree(HERE,target,ignore=shutil.ignore_patterns('__pycache__'))
        (target/'MANIFEST.json').unlink(missing_ok=True)
        regenerated=('EXACT_ALGEBRA.json','REVIEW.json','IMPORTS.json')
        for name in regenerated:(target/name).unlink(missing_ok=True)
        launcher='''from pathlib import Path
import sys,runpy,json
package=Path(sys.argv[1]).resolve();forbidden=Path(sys.argv[2]).resolve()
sys.path.extend(json.loads(sys.argv[3]));reads=set()
def guard(event,args):
    if event=='open'and args and isinstance(args[0],(str,bytes)):
        p=Path(args[0]).resolve()
        if p.is_relative_to(forbidden):raise RuntimeError('Mutable workspace read denied: '+str(p))
        if p.is_relative_to(package)and p.is_file():reads.add(str(p.relative_to(package)))
sys.addaudithook(guard)
sys.path.insert(0,str(package))
runpy.run_path(str(package/'review.py'),run_name='__main__')
(package/'ISOLATED_READS.json').write_text(json.dumps(sorted(reads),indent=2)+'\\n')
'''
        env=dict(os.environ);env.pop('PYTHONPATH',None)
        env.update(PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
        sites=sorted({str(p)for p in (sysconfig.get_path('purelib'),sysconfig.get_path('platlib'),site.getusersitepackages(),*site.getsitepackages())if p and Path(p).is_dir()})
        result=subprocess.run([sys.executable,'-I','-B','-S','-c',launcher,str(target),str(WORKSPACE),json.dumps(sites)],cwd=target,env=env,text=True,capture_output=True)
        (HERE/'ISOLATED_RUN.log').write_text(result.stdout+result.stderr)
        if result.returncode:raise RuntimeError('Independent review replay failed; see ISOLATED_RUN.log')
        for relative,h in inputs.items():assert sha(target/'inputs'/relative)==h
        for name in (*regenerated,'ISOLATED_READS.json'):shutil.copyfile(target/name,HERE/name)
    review=json.loads((HERE/'REVIEW.json').read_text())
    assert review['decision']=='PASS_FULL_LITERAL_CRITERION_ON_DECLARED_NATIVE_ARCHITECTURE'
    for relative,h in inputs.items():assert sha(HERE/'inputs'/relative)==h
    receipt={'status':'PASS','original_mutable_workspace_reads_denied':True,'python':sys.version.split()[0],
      'all_candidate_and_provider_bytes_unchanged':True,'frozen_inputs':len(inputs),
      'deleted_and_regenerated':{name:sha(HERE/name)for name in regenerated},
      'native_large_streams_executed':False,'official_rows_changed':[]}
    (HERE/'PORTABILITY.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    files={str(p.relative_to(HERE)):sha(p)for p in sorted(HERE.rglob('*'))if p.is_file()and not p.is_relative_to(HERE/'inputs')and'__pycache__'not in p.parts and p.name!='MANIFEST.json'}
    manifest={'schema':'oph.row30.full_review.manifest.v1','files_sha256':files,'producer_files':inputs,
      'candidate_manifest_sha256':review['candidate_manifest_sha256'],'decision':review['decision'],
      'review_type':'full_criterion','row':30,'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'manifest_sha256':sha(HERE/'MANIFEST.json'),'review_sha256':sha(HERE/'REVIEW.json'),'files':len(files),'inputs':len(inputs),'decision':review['decision']}),flush=True)

if __name__=='__main__':main()
