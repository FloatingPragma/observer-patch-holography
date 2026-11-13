import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .certificate import calculate

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PARENTS=('openwave_common_hyperfine_spectrum','openwave_strong_step_scaling')


def pins():
    files={str(p.relative_to(ROOT))for p in HERE.iterdir()
           if p.suffix in('.py','.md','.tex')or p.name in('inputs.json','witnesses.json')}
    for parent in PARENTS:
        receipt=ROOT/'model'/parent/'results.json';files.add(str(receipt.relative_to(ROOT)))
        bindings=json.loads(receipt.read_text())['source_sha256']
        for path,digest in bindings.items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
        files.update(bindings)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(files))}


def main():
    p=argparse.ArgumentParser();mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');a=p.parse_args()
    before=pins();result=calculate();assert before==pins();result['source_sha256']=before
    path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),result)
    print(result['verdict'])


if __name__=='__main__':main()
