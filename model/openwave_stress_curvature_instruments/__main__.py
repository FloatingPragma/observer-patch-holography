from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,execute,pin_package
from .source import run

HERE=Path(__file__).resolve().parent
PARENTS=('model/openwave_tensor_record_transport/results.json',)

def pins():
    paths=set(PARENTS)|{'model/particle_geometric_response/PROOF.md','model/particle_field_state_response/PROOF.md'}
    for rel in PARENTS:
        for p,h in json.loads((ROOT/rel).read_text())['source_sha256'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale parent',p)
            paths.add(p)
    return pin_package(HERE,paths)

if __name__=='__main__':execute(HERE,pins,run)
