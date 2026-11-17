#!/usr/bin/env python3
"""Read-only portable custody and bounded scientific checks."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time
HERE=Path(__file__).resolve().parent
def verify(root):
    entries=(root/'MANIFEST.sha256').read_text().splitlines()
    for line in entries:
        h,name=line.split('  ',1);p=root/name
        if p.is_symlink()or root not in p.resolve().parents:raise ValueError('nonlocal evidence path')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise ValueError('hash mismatch: '+name)
    return len(entries)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['verify','headline','independent'],default='headline',nargs='?')
    args=parser.parse_args();start=time.monotonic();binding=json.loads((HERE/'EXPORT_BINDINGS.json').read_text())
    source=HERE/binding['producer_path'];review=HERE/binding['independent_review_path']
    for root,key in [(source,'producer_manifest_sha256'),(review,'independent_review_manifest_sha256')]:
        if hashlib.sha256((root/'MANIFEST.sha256').read_bytes()).hexdigest()!=binding[key]:raise ValueError('manifest binding changed')
    files={'producer':verify(source),'review':verify(review)}
    if(HERE/'MANIFEST.sha256').exists():files['outer']=verify(HERE)
    if args.mode!='verify':
        checkdir=source/'independent_headline'if args.mode=='headline'else review
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONINTMAXSTRDIGITS='0')
        for script in ['check_math.py','check_clock.py']:
            subprocess.run([sys.executable,str(checkdir/script),str(source)],env=env,check=True)
    print(json.dumps({'mode':args.mode,'verified_files':files,'seconds':time.monotonic()-start,
                      'read_only':True,'full_pair_propagation_reexecuted':False,'literal_new_native_calls':0}))
if __name__=='__main__':main()
