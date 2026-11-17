import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .experiment import run
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
    seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')}
    seeds.update(('model/openwave_angular_torque/PROOF.md','model/openwave_rotational_g_limit/PROOF.md'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}
def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    before=pins();value=run();assert pins()==before;value['source_sha256']=before;path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),value),'quantum source replay differs'
    print(value['verdict'])
if __name__=='__main__':main()
