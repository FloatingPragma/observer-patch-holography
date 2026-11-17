from pathlib import Path
import hashlib
from model.particle_common_contract.provenance import import_closure
from .experiment import run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def inputs():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(str(p.relative_to(ROOT)) for p in (HERE/'witnesses').glob('*.npz'))
    seeds.update(('model/openwave_intrinsic_magnetic/PROOF.md','model/openwave_larmor_dynamics/PROOF.md','model/requirements.txt'))
    for p in seeds:
        if not (ROOT/p).is_file():raise FileNotFoundError(p)
    return bound_pins({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))})


def certificate():
    value=run();value['input_sha256']=bound_pins(inputs());return value


def equivalent(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equivalent(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if isinstance(b,float):return isinstance(a,(float,int)) and abs(a-b)<=max(1e-10,abs(b)*1e-7)
    return type(a) is type(b) and a==b
