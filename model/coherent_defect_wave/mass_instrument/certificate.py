"""Source and apparatus custody for exact mass-readout records."""
import hashlib
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from .model import certificate as calculate

ROOT=Path(__file__).resolve().parents[3]


def certificate():
    folder=Path(__file__).parent
    paths={str(p.relative_to(ROOT)) for p in folder.glob('*.py') if not p.name.startswith('test_')}
    paths.update(('model/coherent_defect_wave/mass_instrument/PROOF.md',
                  'model/coherent_defect_wave/many_body_route/PROOF.md',
                  'model/coherent_defect_wave/mobile_flux/ROUTING_COMPLETION.md',
                  'model/requirements.txt'))
    names=import_closure(paths)
    value=calculate()
    value['input_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                           for name in sorted(names) if (ROOT/name).is_file()}
    return value
