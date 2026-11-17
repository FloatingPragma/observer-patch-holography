"""Authenticated current spectral-measure writer and cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent
PARENTS=['model/openwave_full_strong_current/results.json','model/openwave_lambda_wavepacket/results.json','model/openwave_cluster_weyl_sequence/results.json']
def pins():
    paths=set(PARENTS)|{str((HERE/'witnesses.json').relative_to(ROOT))}
    for name in PARENTS:
        parent=json.loads((ROOT/name).read_text())
        for path,expected in bound_pins(parent['source_sha256']).items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==expected,('stale parent',path)
            paths.add(path)
    return pin_package(HERE,paths)
if __name__=='__main__':execute(HERE,pins,run)
