"""Source custody and numerical comparisons for complete precession records."""
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
    seeds.update(('model/openwave_intrinsic_magnetic/PROOF.md','model/openwave_intrinsic_magnetic/results.json','model/requirements.txt'))
    for name in seeds:
        if not (ROOT/name).is_file():raise FileNotFoundError(name)
    return bound_pins({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))})


def certificate():
    value=run();value['input_sha256']=bound_pins(inputs());return value


def equivalent(a,b):
    if isinstance(a,dict) and isinstance(b,dict):
        # Recorded pin maps are compared without prose and test entries.
        a={k:bound_pins(v) if k in('input_sha256','source_sha256') else v for k,v in a.items()}
        b={k:bound_pins(v) if k in('input_sha256','source_sha256') else v for k,v in b.items()}
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equivalent(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if isinstance(b,float):return isinstance(a,(float,int)) and abs(a-b)<=max(2e-9,abs(b)*2e-6)
    return type(a) is type(b) and a==b
