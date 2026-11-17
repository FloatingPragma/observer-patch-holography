"""Recomputed finite common-source and complete reference receipt."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from .model import certificate_data
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins

ROOT=Path(__file__).resolve().parents[3]


def input_hashes():
    here=Path(__file__).parent
    seeds={str(p.relative_to(ROOT)) for p in here.iterdir()
           if p.suffix in ('.py','.md') and not p.name.startswith('test_')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt','model/coherent_defect_wave/chiral_source/PROOF.md',
                  'model/coherent_defect_wave/chiral_regulator/PROOF.md',
                  'model/coherent_defect_wave/chiral_higgs/PROOF.md',
                  'model/coherent_defect_wave/field_identity_audit/PROOF.md'))
    return bound_pins({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(seeds)) if (ROOT/name).is_file()})


def certificate():
    return {'schema':'oph.coherent_defect_wave.chiral_joint.v1',**certificate_data(),
            'input_sha256':bound_pins(input_hashes())}
