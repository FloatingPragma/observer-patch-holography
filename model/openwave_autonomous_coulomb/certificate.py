from pathlib import Path
import hashlib
from model.particle_common_contract.provenance import import_closure
from .experiment import run
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def certificate():
    value=run();seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in ('.py','.md')}
    seeds.update(('model/openwave_pinned_coulomb/PROOF.md','model/openwave_pinned_coulomb/results.json'))
    value['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}
    return value
