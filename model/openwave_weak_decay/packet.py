"""Replay the complete finite source certificate with input custody."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import definition_inputs,definition_hashes
from .decay import certificate
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ANCESTORS=['model/openwave_lepton_mass_shell/PROOF.md', 'model/lepton_mass_compiler/PROOF.md']

def source_bindings():
 seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')}
 return definition_hashes(definition_inputs(seeds|set(ANCESTORS)))

def report():
 pins=source_bindings();result=certificate();result['input_sha256']=pins
 if pins!=source_bindings():raise ValueError('source changed during complete replay')
 return result

def main():
 parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
 group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
 args=parser.parse_args();result=report();path=HERE/'results.json'
 if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
 elif json.loads(path.read_text())!=result:raise ValueError('complete source receipt differs from replay')
 print(result['verdict'])
