"""Original writer and cold replay bound to the common native source."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run
HERE=Path(__file__).resolve().parent
PARENT='model/openwave_native_process/results.json'
def pins():
 paths={PARENT}
 for p,h in json.loads((ROOT/PARENT).read_text())['source_sha256'].items():
  assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale native parent',p)
  paths.add(p)
 return pin_package(HERE,paths)
if __name__=='__main__':execute(HERE,pins,run)
