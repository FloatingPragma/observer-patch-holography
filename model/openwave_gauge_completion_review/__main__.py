from pathlib import Path
import json
from model.openwave_pure_su2_source.custody import execute, pin_package
from .review import run
HERE = Path(__file__).resolve().parent
def pins():
    extra = []
    for name in ('openwave_transverse_completion', 'openwave_quark_flux_records'):
        extra.extend(('model/'+name+'/results.json', 'model/'+name+'/PROOF.md'))
        extra.extend(str(p.relative_to(HERE.parents[1])) for p in (HERE.parent/name).glob('*.py'))
        extra.extend(json.loads((HERE.parent/name/'results.json').read_text())['source_sha256'])
    return pin_package(HERE, extra)
if __name__ == '__main__':
    execute(HERE, pins, run)
