"""Immutable source binding for the finite magnetostatic envelope result."""
from pathlib import Path
import hashlib
from model.particle_common_contract.provenance import import_closure
from .experiment import run
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DEPENDENCIES=('model/native_mass_channels/README.md','model/openwave_field_observables/PROOF.md',
              'model/openwave_magnetic_response/PROOF.md','model/requirements.txt')


def inputs():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(DEPENDENCIES)
    for name in seeds:
        if not (ROOT/name).is_file():raise FileNotFoundError('missing stationary magnetic input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def certificate():
    result=run();result['input_sha256']=inputs();return result


def validate_pins(value):
    for name,expected in value['input_sha256'].items():
        path=Path(name)
        if path.is_absolute() or '..' in path.parts or hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
            raise ValueError('stationary magnetic input changed: '+name)
