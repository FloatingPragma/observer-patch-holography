"""Retained real-check coverage mutations; only disposable evidence copies change."""
import argparse
import copy
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import time
from probe_gaps import OLD, mutations, need

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(folder):
    return {str(p.relative_to(folder)):(sha(p),p.stat().st_size,p.stat().st_mtime_ns,p.stat().st_ino)
            for p in folder.rglob('*') if p.is_file()}


def calculate():
    before=snapshot(OLD)
    evidence=OLD/'records/EVIDENCE.json.gz'
    data=json.loads(gzip.decompress(evidence.read_bytes()))
    results=[]
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    env.pop('PYTHONOPTIMIZE',None)

    def execute(label,path,optimized,expected,cause=None):
        command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(HERE/'check.py')]
        if path is not None:command+=['--evidence',str(path)]
        saved=snapshot(path.parent) if path is not None else snapshot(evidence.parent)
        process=subprocess.run(command,env=env,capture_output=True,text=True,timeout=60)
        need(process.returncode==expected,label+' exit status')
        answer=json.loads(process.stdout);answer.pop('seconds',None)
        need(answer['passed']==(expected==0),label+' scientific verdict')
        if cause:need(cause in answer.get('error',''),label+' wrong rejection cause')
        need(not process.stderr,label+' unexpected stderr')
        need(saved==(snapshot(path.parent) if path is not None else snapshot(evidence.parent)),label+' wrote evidence')
        results.append({'label':label,'optimized':optimized,'exit_code':process.returncode,
                        'result':answer,'evidence_bytes_and_metadata_unchanged':True})

    execute('original_clean',None,False,0)
    execute('original_clean_optimized',None,True,0)
    causes={'drop_column':'incomplete global consumer catalog',
            'duplicate_column':'missing or duplicated global consumer key',
            'drop_input_without_replacement':'missing input cases',
            'duplicate_stream_input':'decoded input/public-stream multiplicity differs',
            'duplicate_other_total_input':'decoded input/public-stream multiplicity differs'}
    with tempfile.TemporaryDirectory(prefix='admission-coverage-validation-') as tmp:
        folder=Path(tmp)
        for label,mutant,_ in mutations(data):
            path=folder/(label+'.json.gz')
            path.write_bytes(gzip.compress(json.dumps(mutant,separators=(',',':')).encode(),compresslevel=1,mtime=0))
            for optimized in (False,True):execute(label,path,optimized,1,causes[label])
            path.unlink()
        # Registered preparation and scientific records stay fixed; irrelevant
        # audit labels and row order may change arbitrarily without affecting it.
        blinded=copy.deepcopy(data)
        for i,row in enumerate(blinded['rows']):
            row.update(seed_audit_only=-i-999,total_audit_only='wrong',filename='fake_Q999.json',run_order=-i)
        random.Random(3617).shuffle(blinded['rows'])
        path=folder/'metadata_blind.json.gz'
        path.write_bytes(gzip.compress(json.dumps(blinded,separators=(',',':')).encode(),compresslevel=1,mtime=0))
        execute('false_labels_shuffled_inputs',path,True,0)
    need(before==snapshot(OLD),'Sealed ancestor packet changed')
    return {'passed':True,'same_group':True,'external_G6':False,'ancestor_unchanged':True,
            'ancestor_manifest':sha(OLD/'MANIFEST.sha256'),'source_evidence_sha256':sha(evidence),
            'checks':results,'fresh_source_runs':0,'native_law_or_controller_changed':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-to',type=Path)
    args=parser.parse_args()
    if args.record_to:need(not args.record_to.exists(),'First-record destination exists')
    start=time.perf_counter();result=calculate()
    if args.record_to:
        with args.record_to.open('x') as stream:stream.write(json.dumps(result,indent=2)+'\n')
    else:need(result==json.loads((HERE/'VALIDATION.json').read_text()),'Saved coverage validation differs')
    print(json.dumps({'passed':True,'checks':len(result['checks']),'seconds':time.perf_counter()-start}))
