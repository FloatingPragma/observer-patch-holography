"""Strict source-bound replay of same-source interaction and calibrated gaps."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run as source_run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent
PARENTS=('model/openwave_native_coherent_clock/results.json','model/openwave_single_valley_leptons/results.json')
def pins():
 paths=set(PARENTS)
 for p in PARENTS:
  for n,h in bound_pins(json.loads((ROOT/p).read_text())['source_sha256']).items():
   assert sha256((ROOT/n).read_bytes()).hexdigest()==h,('stale parent',n)
   paths.add(n)
 return pin_package(HERE,paths|{'model/native_transfer_control/PROOF.md'})
def run():return json.loads(json.dumps(source_run()))
if __name__=='__main__':execute(HERE,pins,run)
