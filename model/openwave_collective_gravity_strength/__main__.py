from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT, pin_package, execute
from .source import run
HERE = Path(__file__).resolve().parent
PARENTS = ['model/openwave_replicated_geometry/results.json', 'model/openwave_gravity_unit_constraint/results.json']
def pins():
    paths = set(PARENTS)
    for name in PARENTS:
        for p, h in json.loads((ROOT/name).read_text())['source_sha256'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, ('stale parent', p)
            paths.add(p)
    paths.add('model/openwave_string_positive_cone/results.json')
    return pin_package(HERE, paths)
if __name__ == '__main__':
    execute(HERE, pins, run)
