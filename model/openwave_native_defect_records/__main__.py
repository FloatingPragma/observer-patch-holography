"""Original exact current-record writer and cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT, pin_package, execute
from .source import run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent


def pins():
    parent = 'model/openwave_native_process/results.json'
    paths = {parent, 'code/a5_closure/native_defect_stability_gate.py',
             'model/defect_protection/protection.py'}
    for name, expected in bound_pins(json.loads((ROOT/parent).read_text())['source_sha256']).items():
        assert sha256((ROOT/name).read_bytes()).hexdigest() == expected, ('stale native process', name)
        paths.add(name)
    return pin_package(HERE, paths)


if __name__ == '__main__':
    execute(HERE, pins, lambda: json.loads(json.dumps(run())))
