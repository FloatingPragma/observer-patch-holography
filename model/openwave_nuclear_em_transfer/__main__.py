"""Authenticated nuclear Gauss-extension writer and cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run
HERE=Path(__file__).resolve().parent
PARENTS=['model/openwave_light_nuclear_curve/results.json','model/openwave_nuclear_isospin/results.json','model/openwave_transverse_completion/results.json']
def pins():
    paths=set(PARENTS)
    for name in PARENTS:
        parent=json.loads((ROOT/name).read_text())
        for path,expected in parent['source_sha256'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==expected,('stale parent',path)
            paths.add(path)
    return pin_package(HERE,paths)
if __name__=='__main__':execute(HERE,pins,run)
