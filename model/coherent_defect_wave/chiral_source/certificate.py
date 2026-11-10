"""Exact receipt for the local chiral carrier and gauge-current attachment."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import certificate_data
from model.common_source.chiral_candidate import action

ROOT=Path(__file__).resolve().parents[3]


def input_hashes():
    here=Path(__file__).parent
    seeds={str(p.relative_to(ROOT)) for p in here.iterdir()
           if p.suffix in ('.py','.md') and not p.name.startswith('test_')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
                  'model/coherent_defect_wave/field_identity_audit/PROOF.md',
                  'model/coherent_defect_wave/physical_attachment/PROOF.md',
                  'model/standard_model_source_map/FERMION_ATTACHMENT.md',
                  'model/standard_model_source_map/GAUGE_ATTACHMENT.md',
                  'model/standard_model_source_map/KINETIC_ATTACHMENT.md',
                  'Lean/Screen/ExteriorSelection.lean','Lean/Screen/WeylYukawaConventions.lean'))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(seeds)) if (ROOT/name).is_file()}


def certificate():
    return {'schema':'oph.coherent_defect_wave.chiral_source.v1',**certificate_data(),
            'input_sha256':input_hashes()}
