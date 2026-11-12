from pathlib import Path
import hashlib
from model.particle_common_contract.provenance import import_closure
from .experiment import run
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def inputs():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(str(p.relative_to(ROOT)) for p in (HERE/'witnesses').glob('*.npz'))
    seeds.update(('model/openwave_intrinsic_magnetic/PROOF.md','model/openwave_larmor_dynamics/PROOF.md','model/requirements.txt'))
    for p in seeds:
        if not (ROOT/p).is_file():raise FileNotFoundError(p)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def certificate():
    value=run();value['input_sha256']=inputs();return value


def equivalent(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equivalent(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if isinstance(b,float):return isinstance(a,(float,int)) and abs(a-b)<=max(1e-10,abs(b)*1e-7)
    return type(a) is type(b) and a==b
