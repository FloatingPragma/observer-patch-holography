"""Reproduce complete cubic spin and field/frame records."""
import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .records import certificate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
    seeds=[str(path.relative_to(ROOT))for path in HERE.iterdir()if path.suffix in ('.py','.md','.tex')]
    seeds+=['model/openwave_lepton_mass_shell/PROOF.md','model/openwave_lepton_mass_shell/results.json','model/openwave_kinetic_phase_cycle/PROOF.md','model/openwave_kinetic_phase_cycle/results.json',
            'model/lepton_mass_compiler/PROOF.md','model/lepton_mass_compiler/results.json']
    return bound_pins({path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()for path in sorted(import_closure(seeds))})
def main():
    parser=argparse.ArgumentParser(description=__doc__);choice=parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--write',action='store_true');choice.add_argument('--check',action='store_true');args=parser.parse_args()
    before=pins();result=certificate();result['source_sha256']=bound_pins(before)
    if pins()!=before:raise ValueError('source changed during spin-frame replay')
    path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    elif json.loads(path.read_text())!=result:raise ValueError('spin-frame receipt differs')
    print(result['verdict'])
if __name__=='__main__':main()
