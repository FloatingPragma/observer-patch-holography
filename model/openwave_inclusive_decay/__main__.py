import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .decay import certificate
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
    seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')}
    seeds.update(('model/openwave_weak_decay/PROOF.md','model/lepton_mass_compiler/PROOF.md',
                  'model/lepton_cycle/experimental_inputs.json'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}
def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    before=pins();value=certificate();assert before==pins();value['source_sha256']=before;path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),value),'inclusive decay replay differs'
    print(value['verdict'])
if __name__=='__main__':main()
