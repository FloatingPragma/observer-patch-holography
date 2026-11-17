"""Source-authenticated producer and separate cold replay."""
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
PARENTS=['model/openwave_su2_volume_gap/results.json']
EXTRA=['Lean/ObserverPatchHolography/SourceScaleCompatibility.lean', 'Lean/ObserverPatchHolography/CapacityStress.lean', 'Lean/ObserverPatchHolography/CapacityNonidentifiability.lean', 'Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean']

def pins():
    paths=set(PARENTS+EXTRA)
    for rel in PARENTS:
        value=json.loads((ROOT/rel).read_text())
        for p,h in bound_pins(value['source_sha256']).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale parent',p)
            paths.add(p)
    return pin_package(HERE,paths)

if __name__=='__main__':
    execute(HERE,pins,run)
