"""Source-authenticated producer and separate cold replay."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run
HERE=Path(__file__).resolve().parent
PARENTS=['model/openwave_dirac_maxwell_refinement/results.json']
EXTRA=[]

def pins():
    paths=set(PARENTS+EXTRA)
    for rel in PARENTS:
        value=json.loads((ROOT/rel).read_text())
        for p,h in value['source_sha256'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale parent',p)
            paths.add(p)
    return pin_package(HERE,paths)

if __name__=='__main__':
    execute(HERE,pins,run)
