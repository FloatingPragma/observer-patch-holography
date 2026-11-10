"""Recomputed finite common-source and complete reference receipt."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from .model import certificate_data

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
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(seeds)) if (ROOT/name).is_file()}


def certificate():
    return {'schema':'oph.coherent_defect_wave.chiral_joint.v1',**certificate_data(),
            'input_sha256':input_hashes()}
