#!/usr/bin/env python3
"""One command for the frozen historical source and checked chain facade."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys
HERE=Path(__file__).resolve().parent
SNAPSHOT=HERE/'frozen_workspace'
REL=Path('openwave-tuning/next/CHAIN/integration')


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def verify():
    manifest=json.loads((HERE/'SNAPSHOT_MANIFEST.json').read_text())
    for relative,digest in manifest['files'].items():
        path=SNAPSHOT/relative
        if path.is_symlink() or not path.is_file() or sha(path)!=digest:
            raise ValueError('frozen workspace changed: '+relative)
    if sha(SNAPSHOT/REL/'SNAPSHOT_MANIFEST.json')!=sha(HERE/'SNAPSHOT_MANIFEST.json'):
        raise ValueError('snapshot manifest copies disagree')
    return manifest


def main():
    manifest=verify()
    env=dict(os.environ);env.pop('PYTHONPATH',None);env['PYTHONDONTWRITEBYTECODE']='1'
    env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']=env['VECLIB_MAXIMUM_THREADS']='1'
    bootstrap="""
import sys,sysconfig,runpy,importlib.util,site
from pathlib import Path
here=Path(sys.argv[1]).resolve()
for key in ('purelib','platlib'):
    value=sysconfig.get_path(key)
    if value not in sys.path:sys.path.append(value)
user=site.getusersitepackages()
if user not in sys.path:sys.path.append(user)
spec=importlib.util.spec_from_file_location('snapshot_guard',here/'snapshot_guard.py')
guard=importlib.util.module_from_spec(spec);sys.modules['snapshot_guard']=guard;spec.loader.exec_module(guard)
guard.install(here)
runpy.run_path(str(here/'integrate.py'),run_name='__main__')
"""
    print('Executing immutable historical workspace; mutable project reads are denied.',flush=True)
    outcome=subprocess.run([sys.executable,'-I','-S','-c',bootstrap,str(SNAPSHOT/REL)],cwd=SNAPSHOT,env=env)
    if outcome.returncode:raise SystemExit(outcome.returncode)
    verify()
    for name in ('RESULT.json','ALL31_GATE.json','CONTRACT.json','MANIFEST.json'):
        shutil.copyfile(SNAPSHOT/REL/name,HERE/name)
    result=json.loads((HERE/'RESULT.json').read_text())
    execution={'schema':'oph.chain.snapshot_execution.v1','exit_code':outcome.returncode,
        'snapshot_manifest_sha256':sha(HERE/'SNAPSHOT_MANIFEST.json'),'snapshot_file_count':manifest['file_count'],
        'all_snapshot_bytes_unchanged':True,'isolation':result['snapshot_execution']['runtime_isolation'],
        'child_manifest_sha256':sha(SNAPSHOT/REL/'MANIFEST.json'),
        'outputs':{name:sha(HERE/name) for name in ('RESULT.json','ALL31_GATE.json','CONTRACT.json')}}
    (HERE/'SNAPSHOT_EXECUTION.json').write_text(json.dumps(execution,indent=2)+'\n')
    package=json.loads((HERE/'MANIFEST.json').read_text())
    package['files']['SNAPSHOT_EXECUTION.json']=sha(HERE/'SNAPSHOT_EXECUTION.json')
    package['snapshot_child_manifest_sha256']=execution['child_manifest_sha256']
    (HERE/'MANIFEST.json').write_text(json.dumps(package,indent=2)+'\n')
    print('Historical snapshot custody and runtime isolation passed.',flush=True)


if __name__=='__main__':main()
