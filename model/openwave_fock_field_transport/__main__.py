from hashlib import sha256
from pathlib import Path
import json
from model.openwave_pure_su2_source.custody import ROOT, execute, pin_package
from .transport import run

HERE = Path(__file__).resolve().parent
PARENTS = ('model/openwave_common_electromagnetism/results.json',
           'model/openwave_dirac_field_connection/results.json')


def pins():
    paths = set(PARENTS)
    paths.add('model/particle_common_contract/PROOF.md')
    for rel in PARENTS:
        for p, h in json.loads((ROOT/rel).read_text())['source_sha256'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, ('stale parent', p)
            paths.add(p)
    return pin_package(HERE, paths)


if __name__ == '__main__':
    execute(HERE, pins, run)
