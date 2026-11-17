"""Original writer and cold replay bound to the common native source."""
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
PARENT='model/openwave_native_process/results.json'
def pins():
 paths={PARENT}
 for p,h in bound_pins(json.loads((ROOT/PARENT).read_text())['source_sha256']).items():
  assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale native parent',p)
  paths.add(p)
 return pin_package(HERE,paths)
if __name__=='__main__':execute(HERE,pins,run)
