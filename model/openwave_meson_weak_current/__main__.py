import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_larmor_dynamics.certificate import equivalent
from .decay import certificate
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
    seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')}
    seeds.update(('model/openwave_native_composites/PROOF.md','model/native_flavour_transport/PROOF.md',
                  'model/openwave_weak_decay/PROOF.md'))
    seeds.update(SOURCE_FILES)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}
def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    before=pins();result=certificate();assert before==pins();result['source_sha256']=before;path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),result)
    print(result['verdict'])
if __name__=='__main__':main()
