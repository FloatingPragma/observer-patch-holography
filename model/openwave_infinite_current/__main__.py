import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .current import certificate
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
 seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')}
 seeds.update(('model/openwave_quantum_rotational/PROOF.md','model/openwave_quantum_rotational/results.json','model/openwave_infinite_photon/PROOF.md','model/openwave_native_quantum_charge/PROOF.md'))
 return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();before=pins();v=certificate();assert pins()==before;v['source_sha256']=before;path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
 else:assert equivalent(json.loads(path.read_text()),v)
 print(v['verdict'])
if __name__=='__main__':main()
