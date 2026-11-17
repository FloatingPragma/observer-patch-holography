"""Read-only, schema-bound replay of the 28 historical b6df jobs.

Original packets remain historical. This runner enforces their assertion-enabled
execution and compares every output with the previously saved clean replay.
It does not repair the original direct producer interfaces or award science.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def pointer(value, path):
    for key in path:value=value[key]
    return value


def records(text):
    decoder=json.JSONDecoder();output=[];text=text.strip()
    while text:
        value,end=decoder.raw_decode(text);output.append(value);text=text[end:].lstrip()
    require(bool(output),'No JSON result')
    return output


def check_output(actual, contract):
    expected=contract['expected_output']
    require(len(actual)==len(expected),'Missing or additional result record')
    for record in actual:
        require(isinstance(record,dict),'Non-object result')
        require(record.get('passed') is not False,'A current result declares failure')
        if isinstance(record.get('gates'),dict):
            require(all(v is True for v in record['gates'].values()),'A current scientific gate failed')
        if 'verdict' in record:require(record['verdict']=='PASS','Current verdict is not PASS')
    for path in contract['required_verdict_paths']:
        require(pointer(actual,path)==pointer(expected,path),'Required verdict differs: '+str(path))
    # Runtime-only fields have explicit prebound paths; no scientific field is dropped.
    normalized=json.loads(json.dumps(actual));reference=json.loads(json.dumps(expected))
    for path in contract['runtime_paths']:
        require(isinstance(pointer(actual,path),(int,float)) and pointer(actual,path)>=0,'Invalid runtime')
        for value in (normalized,reference):
            parent=pointer(value,path[:-1]);parent[path[-1]]=0
    require(json.dumps(normalized,sort_keys=True,allow_nan=False)==
            json.dumps(reference,sort_keys=True,allow_nan=False),'Saved scientific output differs')


def snapshot(folder):
    result={}
    for path in folder.rglob('*'):
        require(not path.is_symlink(),'Source symlink rejected: '+str(path))
        if path.is_file():
            stat=path.stat()
            result[str(path.relative_to(folder))]=(hashlib.sha256(path.read_bytes()).hexdigest(),
                                                    stat.st_size,stat.st_mtime_ns,stat.st_ino)
    return result


def custody(source):
    inventory=json.loads((source/'EVIDENCE.json').read_text());leaves=set()
    for provider in inventory['providers']:
        base=(source/provider['snapshot']).resolve()
        require(base.is_relative_to(source),'Provider escapes source')
        manifest=base/'MANIFEST.sha256'
        require(hashlib.sha256(manifest.read_bytes()).hexdigest()==provider['manifest_sha256'],'Manifest mismatch')
        entries=manifest.read_text().splitlines()
        require(len(entries)==provider['files'],'Wrong leaf count')
        for line in entries:
            digest,name=line.split('  ',1);path=(base/name).resolve()
            require(path.is_relative_to(base),'Leaf escapes provider')
            require(hashlib.sha256(path.read_bytes()).hexdigest()==digest,'Changed source evidence: '+str(path))
            leaves.add(str(path))
    return {'providers':len(inventory['providers']),'distinct_leaves':len(leaves)}


def run_one(source, contract):
    folder=source/'evidence/openwave'/contract['directory']
    command=[sys.executable,'-B',contract['script'],*contract['arguments']]
    # PYTHONOPTIMIZE=0 is inherited by all child Python processes. A parent -O
    # cannot disable checks in these original assertion-based historical jobs.
    env=dict(os.environ,PYTHONOPTIMIZE='0',PYTHONDONTWRITEBYTECODE='1',
             OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    begin=time.monotonic();result={'index':contract['index'],'command':command}
    try:
        process=subprocess.run(command,cwd=folder,env=env,text=True,capture_output=True,timeout=60)
        result.update(exit_code=process.returncode,stdout=process.stdout,stderr=process.stderr)
        require(process.returncode==0,'Subprocess failed')
        require(not process.stderr,'Unexpected stderr')
        output=records(process.stdout);check_output(output,contract)
        result.update(output=output,passed=True)
    except (AssertionError,ValueError,KeyError,TypeError,subprocess.TimeoutExpired) as error:
        result.update(passed=False,error=str(error))
    result['seconds']=time.monotonic()-begin
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--contracts',type=Path,default=HERE/'SOURCE_JOB_CONTRACTS.json')
    parser.add_argument('--job',type=int,action='append')
    parser.add_argument('--record-to',type=Path,help='Explicit exclusive creation of a new replay receipt')
    args=parser.parse_args();source=args.source.resolve()
    if args.record_to:
        require(not args.record_to.resolve().is_relative_to(source),'Receipt must be outside saved source')
        require(not args.record_to.exists(),'First-record destination already exists')
    contracts=json.loads(args.contracts.read_text())['jobs']
    if args.job is not None:
        require(len(set(args.job))==len(args.job) and all(0<=i<len(contracts) for i in args.job),'Invalid job selection')
        contracts=[contracts[i] for i in args.job]
    before=snapshot(source);checked=custody(source);begin=time.monotonic()
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs=list(pool.map(lambda item:run_one(source,item),contracts))
    unchanged=before==snapshot(source)
    output={'jobs':jobs,'custody':checked,'saved_evidence_unchanged':unchanged,
            'passed':unchanged and all(j['passed'] for j in jobs),'seconds':time.monotonic()-begin,
            'same_group':True,'external_G6':False,'mode':'read-only saved-evidence comparison'}
    payload=json.dumps(output,indent=2,allow_nan=False)+'\n'
    if args.record_to:
        with args.record_to.open('x') as stream:stream.write(payload)
    print(payload,end='')
    raise SystemExit(0 if output['passed'] else 1)


if __name__=='__main__':main()
