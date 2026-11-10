"""Recomputed actual-current source binding, with pinned theorem interfaces."""
import hashlib
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from .model import dark_module,face_charge_obstruction,bound_transport,antisymmetric_clouds,mass_freedom

ROOT=Path(__file__).resolve().parents[3]


def input_hashes():
    here=Path(__file__).parent
    seeds={str(p.relative_to(ROOT)) for p in here.iterdir()
           if p.suffix in ('.py','.md') and not p.name.startswith('test_')}
    seeds.update(('model/requirements.txt','model/coherent_defect_wave/species/PROOF.md',
                  'model/coherent_defect_wave/native_routing_join/PROOF.md',
                  'model/coherent_defect_wave/physical_attachment/PROOF.md',
                  'model/coherent_defect_wave/mobile_flux/ROUTING_COMPLETION.md',
                  'model/standard_model_source_map/FERMION_ATTACHMENT.md',
                  'model/standard_model_source_map/GAUGE_ATTACHMENT.md',
                  'model/standard_model_source_map/candidates/PROOF.md',
                  'Lean/Screen/ExteriorSelection.lean','Lean/Screen/FermionSectorAssembly.lean',
                  'Lean/Screen/WeylYukawaConventions.lean','Lean/Screen/ChargedCurrentDictionary.lean',
                  'Lean/Geometry/MassShellKinematics.lean','Lean/QFT/InheritanceMatrix.lean'))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(seeds)) if (ROOT/name).is_file()}


def certificate():
    return {'schema':'oph.coherent_defect_wave.field_identity_audit.v1',
            'native_dark_current_module':dark_module(),
            'face_charge_obstruction':face_charge_obstruction(),
            'bound_composite_transport':bound_transport(),
            'whole_object_exterior_extension':antisymmetric_clouds(),
            'mass_normalization_and_chirality':mass_freedom(),
            'scope':{'actual_native_R15_current_binding':True,
                     'finite_complete_charge_preserving_bound_transport':True,
                     'native_idle_law_selects_binding':False,
                     'chiral_SM_particle_identification':False,
                     'interacting_relativistic_limit':False,
                     'source_GeV_masses':False},
            'input_sha256':input_hashes()}
