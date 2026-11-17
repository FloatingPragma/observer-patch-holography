"""Produce and verify the complete Higgs-Wilson source receipt."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent
from .spectrum import certificate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def pins():
    seeds=[str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')]
    seeds+=['model/openwave_full_higgs_floor/PROOF.md','model/openwave_lepton_mass_shell/PROOF.md']
    return bound_pins({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))})
def main():
    ap=argparse.ArgumentParser();g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=ap.parse_args()
    before=pins();v=certificate();assert before==pins();v['source_sha256']=bound_pins(before);p=HERE/'results.json'
    if a.write:p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
    elif not equivalent(json.loads(p.read_text()),v):raise ValueError('single-valley receipt replay differs')
    print(v['verdict'])
if __name__=='__main__':main()
