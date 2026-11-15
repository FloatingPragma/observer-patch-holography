"""Strict source-bound replay of same-source interaction and calibrated gaps."""
from pathlib import Path
from hashlib import sha256
import json
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from .source import run as source_run
HERE=Path(__file__).resolve().parent
PARENTS=('model/openwave_native_coherent_clock/results.json','model/openwave_single_valley_leptons/results.json')
def pins():
 paths=set(PARENTS)
 for p in PARENTS:
  for n,h in json.loads((ROOT/p).read_text())['source_sha256'].items():
   assert sha256((ROOT/n).read_bytes()).hexdigest()==h,('stale parent',n)
   paths.add(n)
 return pin_package(HERE,paths|{'model/native_transfer_control/PROOF.md'})
def run():return json.loads(json.dumps(source_run()))
if __name__=='__main__':execute(HERE,pins,run)
