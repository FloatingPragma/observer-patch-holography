"""Source custody for native Gauss projection and the finite energy choice."""
from pathlib import Path
import hashlib
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from .result import certificate as result
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DEPENDENCIES=tuple('model/'+name+'/'+file for name in
    ('openwave_reference_binding','openwave_native_charge_flux','openwave_field_observables')
    for file in ('PROOF.md','results.json'))


def input_hashes():
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(DEPENDENCIES);seeds.update(SOURCE_FILES);seeds.add('model/requirements.txt')
    for name in seeds:
        if not (ROOT/name).is_file():raise FileNotFoundError('missing native current/Coulomb input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def certificate():
    value=result();value['input_sha256']=input_hashes();return value


def validate_pins(value):
    for name,expected in value['input_sha256'].items():
        path=Path(name)
        if path.is_absolute() or '..' in path.parts or hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
            raise ValueError('native current/Coulomb source changed: '+name)
