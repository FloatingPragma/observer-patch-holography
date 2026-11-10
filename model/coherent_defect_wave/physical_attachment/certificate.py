"""Source custody for the spin, current and chiral-mass attachments."""
import hashlib
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from . import spin,chiral,compensator

ROOT=Path(__file__).resolve().parents[3]


def input_hashes():
    local=Path(__file__).parent
    seeds={str(p.relative_to(ROOT)) for p in local.glob('*.py') if not p.name.startswith('test_')}
    seeds.update(str((local/name).relative_to(ROOT)) for name in ('PROOF.md','README.md'))
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt','model/coherent_defect_wave/species/PROOF.md',
                  'model/coherent_defect_wave/spatial/ROUTE.md',
                  'model/standard_model_source_map/FERMION_ATTACHMENT.md',
                  'model/standard_model_source_map/kinetic_attachment.py',
                  'Lean/Geometry/CanonicalLorentzModule.lean',
                  'Lean/Geometry/MassShellKinematics.lean',
                  'evidence/native_geometric_source_20260925/spec.json',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(seeds)) if (ROOT/name).is_file()}


def certificate():
    return {'schema':'oph.coherent_defect_wave.physical_attachment.v1',
            'input_sha256':input_hashes(),
            'spin_and_current_attachment':spin.certificate(),
            'chiral_internal_attachment':chiral.certificate(),
            'charge_compensating_memory':compensator.certificate(),
            'scope':{'declared_free_Dirac_route_spin_half':True,
                     'declared_route_exact_finite_charge_continuity':True,
                     'source_internal_representation_mass_constraints':True,
                     'native_U_relativistic_spin_identification':False,
                     'interacting_Lorentz_limit':False,
                     'dynamical_chiral_gauge_theory':False,
                     'SM_particle_identification':False,'GeV_masses':False}}
