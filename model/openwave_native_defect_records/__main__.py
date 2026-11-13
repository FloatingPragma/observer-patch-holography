"""Original exact current-record writer and cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT, pin_package, execute
from .source import run

HERE = Path(__file__).resolve().parent


def pins():
    parent = 'model/openwave_native_process/results.json'
    paths = {parent, 'code/a5_closure/native_defect_stability_gate.py',
             'model/defect_protection/protection.py'}
    for name, expected in json.loads((ROOT/parent).read_text())['source_sha256'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest() == expected, ('stale native process', name)
        paths.add(name)
    return pin_package(HERE, paths)


if __name__ == '__main__':
    execute(HERE, pins, lambda: json.loads(json.dumps(run())))
