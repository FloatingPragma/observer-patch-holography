"""Actual source producer with independent complete cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import execute,pin_package
from .source import run
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PARENTS=['model/openwave_light_nuclear_curve/results.json', 'model/openwave_string_positive_cone/results.json']

def pins():
    paths=set(PARENTS)
    for name in PARENTS:
        parent=json.loads((ROOT/name).read_text())
        for path,expected in parent['source_sha256'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==expected,('stale parent',path)
            paths.add(path)
    paths.update(['model/openwave_transverse_gravity/inputs.json','model/openwave_lepton_widths/inputs.json'])
    return pin_package(HERE,paths)

if __name__=='__main__':
    execute(HERE,pins,run)
