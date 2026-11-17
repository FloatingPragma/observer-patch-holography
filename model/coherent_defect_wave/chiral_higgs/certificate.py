"""Primary source custody for the exact chiral Higgs response."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import certificate as calculate
from .apparatus import certificate as apparatus_certificate
from model.common_source.chiral_candidate.action import SOURCE_FILES

ROOT=Path(__file__).resolve().parents[3]


def certificate():
    folder=Path(__file__).parent
    paths={str(p.relative_to(ROOT)) for p in folder.glob('*.py') if not p.name.startswith('test_')}
    paths.update(('model/coherent_defect_wave/chiral_higgs/PROOF.md',
                  'model/coherent_defect_wave/chiral_higgs/GAUSS_APPARATUS.md',
                  'model/coherent_defect_wave/physical_attachment/PROOF.md',
                  'model/coherent_defect_wave/field_identity_audit/PROOF.md',
                  'model/coherent_defect_wave/chiral_source/PROOF.md',
                  'model/coherent_defect_wave/chiral_source/model.py',
                  'model/coherent_defect_wave/chiral_regulator/PROOF.md',
                  'model/coherent_defect_wave/chiral_regulator/model.py',
                  'model/requirements.txt'))
    paths.update(SOURCE_FILES)
    value=calculate()
    value['compensated_gauge_apparatus']=apparatus_certificate(value['quantum_spectrum']['complete_driven_reservoir_energy_records'])
    value['input_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                           for name in sorted(import_closure(paths)) if (ROOT/name).is_file()}
    return value
