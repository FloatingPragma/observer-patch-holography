"""Replay finite source bridges using only the official local RER package."""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('Bridge verification requires enabled assertions')
import argparse
import hashlib
import json
from pathlib import Path
import time
from ..verify import compare, read_result, scientific_inputs as runtime_inputs

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent
ROOT=PACKAGE.parents[1]


def scientific_inputs():
    paths=[p for p in HERE.rglob('*') if p.suffix in ('.py','.lean','.md') and not p.name.startswith('test_')]
    paths+=[HERE/'LEAN_RESULTS.json', PACKAGE/'RESULTS.json', ROOT/'Lean/QFT/CoarseEnergyWard.lean',
            ROOT/'Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean']
    result={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    for name,digest in runtime_inputs().items():
        result[str((PACKAGE/name).resolve().relative_to(ROOT))]=digest
    from .verify_lean import inputs as lean_inputs
    result.update(lean_inputs())
    return result


def run():
    from .verify_lean import validate_saved as validate_lean
    validate_lean()
    before=scientific_inputs()
    runtime=read_result(PACKAGE/'RESULTS.json')
    if runtime.get('scientific_inputs_sha256') != runtime_inputs() or not runtime.get('all_passed'):
        raise ValueError('Current full native runtime receipt required')
    if runtime.get('specification') != 'A-ccf2-local-capacity-v2 finite matter, field and response construction':
        raise ValueError('Runtime receipt must identify the local-capacity law')
    from .kinematics import bridge as kinematics, phase
    from .matter import bridge as matter
    from .reference import check as reference
    from .gravity import coframe, refinement, ward
    started=time.monotonic()
    result={'model':'CCF-LT-2','candidate':'A-ccf2-local-capacity-v2',
            'kinematics':kinematics.run(),'phase_instrument':phase.run()}
    result['matter'],result['matter_columns']=matter.execute()
    result['reference']={'actual_source':reference.actual_reference(),'complete_fixture':reference.fixture(),
                         'stationarity_boundary':reference.stationarity_boundary(result['phase_instrument'])}
    result['gravity']={'coframe':coframe.audit(),'refinement':refinement.audit(),'energy_ward':ward.audit()}
    result['physical_rows_awarded']=[]
    result['all_passed']=True
    if scientific_inputs()!=before: raise ValueError('Scientific inputs changed during bridge execution')
    result['scientific_inputs_sha256']=before
    result['elapsed_seconds']=time.monotonic()-started
    return json.loads(json.dumps(result))


def validate_saved(path=None):
    result=read_result(path or HERE/'RESULTS.json')
    if result.get('scientific_inputs_sha256')!=scientific_inputs():
        raise ValueError('Stale bridge scientific inputs')
    if result.get('model')!='CCF-LT-2' or result.get('candidate')!='A-ccf2-local-capacity-v2' or result.get('all_passed') is not True:
        raise ValueError('Invalid bridge receipt')
    if result.get('physical_rows_awarded')!=[]:
        raise ValueError('Finite bridge cannot award physical rows')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--output',type=Path)
    group.add_argument('--verify',type=Path)
    args=parser.parse_args()
    result=run()
    if args.output:args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:compare(read_result(args.verify or HERE/'RESULTS.json'),result)
    print(json.dumps({'all_passed':True,'candidate':result['candidate'],'elapsed_seconds':result['elapsed_seconds']},indent=2))


if __name__=='__main__':main()
