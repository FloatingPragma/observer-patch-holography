"""Actual source producer with complete parent custody and cold replay."""
from pathlib import Path
import hashlib,json
from model.openwave_pure_su2_source.custody import pin_package,execute,ROOT
from .source import certificate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent
assert HERE==ROOT/'model/openwave_clock_phase_resource', 'canonical source path required'
PARENTS=['model/openwave_autonomous_clock_channel/results.json', 'model/openwave_passive_clock_records/results.json']

def pins():
    extra=['Lean/ObserverPatchHolography/ClockPhaseCovariance.lean','model/openwave_clock_phase_resource/Transfer.lean']
    for rel in PARENTS:
        extra.append(rel)
        r=json.loads((ROOT/rel).read_text())
        for name,expected in bound_pins(r['source_sha256']).items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,('stale parent',name)
            extra.append(name)
    extra.extend(str(p.relative_to(ROOT))for p in HERE.glob('*.json')if p.name!='results.json')
    return pin_package(HERE,extra)

if __name__=='__main__':execute(HERE,pins,certificate)
