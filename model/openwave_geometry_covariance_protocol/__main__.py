"""Original source-bound producer and cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT, pin_package, execute
from .source import run
HERE=Path(__file__).resolve().parent
PARENTS=('model/openwave_gravity_noise_normalization/results.json',
         'model/openwave_fock_field_transport/results.json')
def pins():
    paths=set(PARENTS)
    for parent in PARENTS:
        for path,h in json.loads((ROOT/parent).read_text())['source_sha256'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==h,('stale parent',path)
            paths.add(path)
    return pin_package(HERE,paths)
if __name__=='__main__':execute(HERE,pins,run)
