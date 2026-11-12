"""Source custody and numerical comparisons for complete precession records."""
from pathlib import Path
import hashlib
from model.particle_common_contract.provenance import import_closure
from .experiment import run
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def inputs():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(('model/openwave_intrinsic_magnetic/PROOF.md','model/openwave_intrinsic_magnetic/results.json','model/requirements.txt'))
    for name in seeds:
        if not (ROOT/name).is_file():raise FileNotFoundError(name)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def certificate():
    value=run();value['input_sha256']=inputs();return value


def equivalent(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equivalent(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if isinstance(b,float):return isinstance(a,(float,int)) and abs(a-b)<=max(2e-9,abs(b)*2e-6)
    return type(a) is type(b) and a==b
