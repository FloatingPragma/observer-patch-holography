"""Original writer and complete-source cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,execute,pin_package
from .source import run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent
PARENTS=('model/openwave_active_native_charge/results.json','model/openwave_native_charge_readout/results.json','model/openwave_coherent_native_architecture/results.json','model/openwave_coherent_native_propagation/results.json')
def pins():
    paths=set(PARENTS)
    for rel in PARENTS:
        for p,h in bound_pins(json.loads((ROOT/rel).read_text())['source_sha256']).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale inherited source',p)
            paths.add(p)
    return pin_package(HERE,paths)
if __name__=='__main__':execute(HERE,pins,run)
