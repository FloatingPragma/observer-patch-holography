"""Original writer and exact cold replay for the active charge construction."""
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

def pins():
    paths={'model/exchange_sector/statistics.py','model/defect_protection/protection.py'}
    for parent in ('model/openwave_native_process/results.json','model/openwave_native_defect_records/results.json'):
        paths.add(parent)
        for name,h in bound_pins(json.loads((ROOT/parent).read_text())['source_sha256']).items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,('stale native parent',name)
            paths.add(name)
    return pin_package(HERE,paths)

if __name__=='__main__':
    execute(HERE,pins,lambda:json.loads(json.dumps(run())))
