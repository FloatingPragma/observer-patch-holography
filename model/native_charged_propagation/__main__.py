"""Original writer and complete-source cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,execute,pin_package
from .source import run
HERE=Path(__file__).resolve().parent
PARENTS=('model/openwave_active_native_charge/results.json','model/openwave_native_charge_readout/results.json','model/openwave_coherent_native_architecture/results.json','model/openwave_coherent_native_propagation/results.json')
def pins():
    paths=set(PARENTS)
    for rel in PARENTS:
        for p,h in json.loads((ROOT/rel).read_text())['source_sha256'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale inherited source',p)
            paths.add(p)
    return pin_package(HERE,paths)
if __name__=='__main__':execute(HERE,pins,run)
