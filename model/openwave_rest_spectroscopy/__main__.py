import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .certificate import calculate as certificate

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PARENTS=('openwave_spin_gradient_nuclei','openwave_dressed_baryon_lines')

def pins():
    files={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')or p.name in('inputs.json','spectral_witnesses.json')}
    for package in PARENTS:
        receipt=ROOT/'model'/package/'results.json';files.add(str(receipt.relative_to(ROOT)))
        files.update(json.loads(receipt.read_text())['source_sha256'])
    files=import_closure(files)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(files)}

def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    before=pins();result=certificate();assert before==pins();result['source_sha256']=before;path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),result)
    print(result['verdict'])

if __name__=='__main__':main()
