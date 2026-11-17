#!/usr/bin/env python3
"""Offline custody, bounded reconstruction and isolated native-reader replay."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time
HERE=Path(__file__).resolve().parent
REL=Path('openwave/row_02_electron_rest_energy/work/2026-11-15_codex_source_repair_v2')
SOURCE=HERE/'workspace'/REL
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def verify():
    bound=read(HERE/'EXPORT_BINDINGS.json')
    assert sha(SOURCE/'MANIFEST.sha256')==bound['original_producer_manifest_sha256']
    for line in (SOURCE/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);assert sha(SOURCE/n)==h,n
    deps=read(SOURCE/'DEPENDENCIES.json')
    for n,m in deps['final_native_coupled_and_routed_replay_files'].items():
        assert sha(HERE/'workspace/openwave'/n)==m['sha256'],n
    for m in deps['historical_fullsource_compilation_dependency_maps'].values():
        assert sha(HERE/'workspace/openwave'/m['path'])==m['sha256']
        for n,h in m['files'].items():assert sha(HERE/'workspace'/n)==h,n
    for n,m in bound['historical_archive_files'].items():
        assert not (HERE/n).is_symlink();assert sha(HERE/n)==m['sha256']
    for n,m in bound['complete_suppliers_and_reviews'].items():
        directory=HERE/n;assert sha(directory/'MANIFEST.sha256')==m['original_manifest_sha256']
        for line in (directory/'MANIFEST.sha256').read_text().splitlines():
            h,f=line.split('  ',1);assert sha(directory/f)==h,(n,f)
    return bound,deps
def run(script,cwd):
    env=os.environ.copy()
    for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):env[k]='1'
    result=subprocess.run([sys.executable,script],cwd=cwd,env=env,
                          capture_output=True,text=True,timeout=55)
    assert result.returncode==0,(script,result.returncode,result.stdout,result.stderr)
    return {'exit_code':result.returncode,'stdout':result.stdout,'stderr':result.stderr}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['verify','headline','independent','native'])
    args=parser.parse_args();start=time.monotonic();binding,deps=verify()
    result={'mode':args.mode,'custody_verified':True,
        'original_producer_manifest_sha256':binding['original_producer_manifest_sha256']}
    if args.mode=='headline':
        completed=run('headline.py',SOURCE);result['execution']=completed
        result['headline']=json.loads(completed['stdout'])
    elif args.mode=='independent':
        review=HERE/'workspace/openwave/row_02_electron_rest_energy/reviews/2026-11-16_codex_energy_join'
        result['independent_execution']=run('check.py',review)
    elif args.mode=='native':
        # Only the explicitly bounded 98-file native-reader dependency tree is
        # copied. The archived millions-of-pulse execution is not replayed.
        with tempfile.TemporaryDirectory(prefix='oph-energy-portable-') as temporary:
            root=Path(temporary);target=root/REL
            retained=str(REL.relative_to('openwave')/'coupled_energy/RECEIPT.json')
            for name,info in deps['final_native_coupled_and_routed_replay_files'].items():
                if name==retained:continue
                src=HERE/'workspace/openwave'/name;dst=root/'openwave'/name
                dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
            one=run('run_coupled_energy.py',target)
            expected=read(SOURCE/'coupled_energy/RECEIPT.json')
            actual=read(target/'coupled_energy/RECEIPT.json')
            assert actual['scientific']==expected['scientific']
            assert actual['artifacts']==expected['artifacts']
            # The frozen routed stage binds elapsed metadata too. Restore that
            # receipt only after matching its scientific content and NPZ hash.
            shutil.copy2(SOURCE/'coupled_energy/RECEIPT.json',target/'coupled_energy/RECEIPT.json')
            two=run('join_routed_source.py',target)
            expected=read(SOURCE/'routed_join/RECEIPT.json');actual=read(target/'routed_join/RECEIPT.json')
            assert actual['scientific']==expected['scientific']
            assert actual['maps_sha256']==expected['maps_sha256']
            result.update({'native_execution':one,'routed_execution':two,
                'scientific_json_identical':True,'npz_bytes_identical':True,
                'fresh_full_native_flow_claimed':False})
    result['seconds']=time.monotonic()-start
    result['environment']={'python':sys.version,'executable':sys.executable,
        'numpy':__import__('numpy').__version__,'scipy':__import__('scipy').__version__,
        'BLAS_threads':1}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
