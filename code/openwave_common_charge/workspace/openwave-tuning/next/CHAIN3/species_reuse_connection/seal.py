"""Verify inputs, replay in an isolated copy, and seal this component scope."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import site
import subprocess
import sys
import sysconfig
import tempfile

HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[3]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inputs=json.loads((HERE/'INPUTS.json').read_text())['files_sha256']
    for relative,h in inputs.items():assert sha(HERE/'inputs'/relative)==h,relative
    with tempfile.TemporaryDirectory(prefix='oph-species-replay-')as temporary:
        target=Path(temporary)/'workspace/openwave-tuning/next/CHAIN3/species_reuse_connection'
        shutil.copytree(HERE,target,ignore=shutil.ignore_patterns('__pycache__'))
        (target/'MANIFEST.json').unlink(missing_ok=True)
        for name in ('REGISTER.json','RESULT.json','NATIVE_PROGRAMS.json','NATIVE_RECORDS.json','FULL90_SPECTRUM.json','RUNTIME_INPUTS.json'):
            (target/name).unlink(missing_ok=True)
        launcher='''from pathlib import Path
import sys, runpy, json
package=Path(sys.argv[1]).resolve(); forbidden=Path(sys.argv[2]).resolve()
sys.path.extend(json.loads(sys.argv[3]))
reads=set()
def guard(event,args):
    if event=='open' and args and isinstance(args[0],(str,bytes)):
        p=Path(args[0]).resolve()
        if p.is_relative_to(forbidden): raise RuntimeError('Mutable workspace read denied: '+str(p))
        if p.is_relative_to(package) and p.is_file(): reads.add(str(p.relative_to(package)))
sys.addaudithook(guard)
sys.path.insert(0,str(package))
runpy.run_path(str(package/'run.py'),run_name='__main__')
import json
(package/'ISOLATED_READS.json').write_text(json.dumps(sorted(reads),indent=2)+'\\n')
'''
        environment=dict(os.environ);environment.pop('PYTHONPATH',None)
        environment.update(PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
        runtime_sites=sorted({str(p)for p in (sysconfig.get_path('purelib'),sysconfig.get_path('platlib'),site.getusersitepackages(),*site.getsitepackages())if p and Path(p).is_dir()})
        p=subprocess.run([sys.executable,'-I','-B','-S','-c',launcher,str(target),str(WORKSPACE),json.dumps(runtime_sites)],
             cwd=target,env=environment,text=True,capture_output=True)
        (HERE/'ISOLATED_RUN.log').write_text(p.stdout+p.stderr)
        if p.returncode:raise RuntimeError('isolated replay failed; see ISOLATED_RUN.log')
        for relative,h in inputs.items():assert sha(target/'inputs'/relative)==h,relative
        regenerated={}
        for name in ('REGISTER.json','RESULT.json','NATIVE_PROGRAMS.json','NATIVE_RECORDS.json','FULL90_SPECTRUM.json','RUNTIME_INPUTS.json'):
            shutil.copyfile(target/name,HERE/name)
            regenerated[name]=sha(HERE/name)
        shutil.copyfile(target/'ISOLATED_READS.json',HERE/'ISOLATED_READS.json')
    result=json.loads((HERE/'RESULT.json').read_text())
    assert result['status']=='PASS_SCOPED_FINITE_REGISTER_AND_NATIVE_MASS_SLOT'
    assert result['native_mass_and_clock']['actual_provider_suffix_first48operations_equal']
    assert result['native_mass_and_clock']['mass_and_clock_both_have_native_C_E_fine_readout']
    assert len(result['flavored_field_complete_bounds'])==2
    for relative,h in inputs.items():assert sha(HERE/'inputs'/relative)==h
    receipt={'status':'PASS','isolated_python':sys.version.split()[0],
        'original_mutable_workspace_reads_denied':True,'all_frozen_inputs_unchanged':True,
        'deleted_and_regenerated':regenerated,'frozen_input_files':len(inputs),
        'full_field_word_executed':False,'physical_rows_validated':[]}
    (HERE/'PORTABILITY.json').write_text(json.dumps(receipt,indent=2)+'\n')
    local={str(p.relative_to(HERE)):sha(p)for p in sorted(HERE.rglob('*'))if p.is_file()and p.name!='MANIFEST.json'and'__pycache__'not in p.parts and not p.is_relative_to(HERE/'inputs')}
    manifest={'schema':'oph.species_connection.freeze.v1','files_sha256':local,
        'producer_files':inputs,'physical_rows_validated':[],'official_rows_changed':[],
        'scope':result['status']}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'manifest_sha256':sha(HERE/'MANIFEST.json'),'result_sha256':sha(HERE/'RESULT.json'),
                     'frozen_files':len(inputs),'local_files':len(local)}),flush=True)


if __name__=='__main__':main()
