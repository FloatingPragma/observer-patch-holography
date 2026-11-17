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
def certificate():
    value=run();seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in ('.py','.md')}
    seeds.update(('model/openwave_pinned_coulomb/PROOF.md','model/openwave_pinned_coulomb/results.json'))
    value['source_sha256']=bound_pins({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))})
    return value
