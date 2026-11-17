from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT, execute, pin_package
from .source import run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
PARENTS = ('model/openwave_native_response/results.json',
           'model/openwave_native_event_clock/results.json')


def pins():
    paths = set(PARENTS)
    for name in PARENTS:
        for path, digest in bound_pins(json.loads((ROOT/name).read_text())['source_sha256']).items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
            paths.add(path)
    return pin_package(HERE, paths)


if __name__ == '__main__':
    execute(HERE, pins, run)
