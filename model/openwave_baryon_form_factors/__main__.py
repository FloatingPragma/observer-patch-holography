import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .source import certificate

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def pins():
    paths={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')or p.name=='inputs.json'}
    for name in ('openwave_common_hyperfine_spectrum',):
        p=ROOT/'model'/name/'results.json';paths.add(str(p.relative_to(ROOT)))
        parent=json.loads(p.read_text())
        for key,digest in parent['source_sha256'].items():
            assert hashlib.sha256((ROOT/key).read_bytes()).hexdigest()==digest,(name,key)
            paths.add(key)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(paths))}

def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();before=pins();value=certificate();assert before==pins()
    value['source_sha256']=before;p=HERE/'results.json'
    if args.write:p.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(p.read_text()),value)
    print(value['verdict'])

if __name__=='__main__':main()
