"""Source and apparatus custody for exact mass-readout records."""
import hashlib
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from .model import certificate as calculate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins

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
    value['input_sha256']=bound_pins({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                           for name in sorted(names) if (ROOT/name).is_file()})
    return value
