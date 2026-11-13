"""Strict source-bound writer and cold replay for repeated clock records."""
from pathlib import Path
import json
from hashlib import sha256
from model.openwave_pure_su2_source.custody import ROOT,pin_package,execute
from model.common_source import repair_continuation as rc
from .source import run as source_run
HERE=Path(__file__).resolve().parent
def pins():
 paths=set(rc.source_bindings())|{'model/native_transfer_control/PROOF.md','model/native_effect_synthesis/PROOF.md'}
 for parent in ('model/openwave_native_charge_readout/results.json','model/openwave_coherent_native_architecture/results.json'):
  paths.add(parent)
  for p,h in json.loads((ROOT/parent).read_text())['source_sha256'].items():
   assert sha256((ROOT/p).read_bytes()).hexdigest()==h,('stale parent',p)
   paths.add(p)
 return pin_package(HERE,paths)
def run():return json.loads(json.dumps(source_run()))
if __name__=='__main__':execute(HERE,pins,run)
