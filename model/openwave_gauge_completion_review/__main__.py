from pathlib import Path
import json
from model.openwave_pure_su2_source.custody import execute, pin_package
from .review import run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE = Path(__file__).resolve().parent
def pins():
    extra = []
    for name in ('openwave_transverse_completion', 'openwave_quark_flux_records'):
        extra.extend(('model/'+name+'/results.json', 'model/'+name+'/PROOF.md'))
        extra.extend(str(p.relative_to(HERE.parents[1])) for p in (HERE.parent/name).glob('*.py'))
        extra.extend(bound_pins(json.loads((HERE.parent/name/'results.json').read_text())['source_sha256']))
    return pin_package(HERE, extra)
if __name__ == '__main__':
    execute(HERE, pins, run)
