"""Run the complete transverse construction and independently replay it."""
import argparse,hashlib,json,math
from pathlib import Path
from .records import run
from model.particle_common_contract.provenance import import_closure
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def pins():
    paths=set()
    for parent in ('model/openwave_common_gravity_records/results.json','model/openwave_universal_energy_geometry/results.json','model/openwave_string_positive_cone/results.json','model/openwave_common_hyperfine_spectrum/results.json'):
        paths.add(parent)
        for p,h in json.loads((ROOT/parent).read_text())['source_sha256'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale parent',p)
            paths.add(p)
    paths|={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex','.json')and p.name!='results.json'}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(paths))}


def equivalent(a,b):
    if isinstance(b,dict):return isinstance(a,dict)and set(a)==set(b)and all(equivalent(a[k],v)for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list)and len(a)==len(b)and all(equivalent(x,y)for x,y in zip(a,b,strict=True))
    if isinstance(b,float):return type(a)is float and math.isfinite(a)and math.isfinite(b)and abs(a-b)<=2e-9+1e-9*abs(b)
    return type(a)is type(b)and a==b


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    before=pins();r=run();assert before==pins();r['source_sha256']=before
    out=HERE/'results.json'
    if a.write:out.write_text(json.dumps(r,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:assert equivalent(json.loads(out.read_text()),r),'full transverse producer replay mismatch'
    print(r['verdict'])
if __name__=='__main__':main()
