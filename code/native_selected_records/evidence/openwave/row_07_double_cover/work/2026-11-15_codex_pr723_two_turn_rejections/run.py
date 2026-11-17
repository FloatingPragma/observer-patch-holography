#!/usr/bin/env python3
"""Bound supplement; preserve the original one-slow-turn failure ledger."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
BOUND=HERE.parent/'2026-11-15_codex_pr723_rotation/berry'
MANIFEST='8add06ee1425a84e00997a7311f37ba222e7f77199487040506838cc21f00b7f'

def calculate():
    assert hashlib.sha256((BOUND/'MANIFEST.sha256').read_bytes()).hexdigest()==MANIFEST
    for line in(BOUND/'MANIFEST.sha256').read_text().splitlines():
        h,name=line.split('  ',1);assert hashlib.sha256((BOUND/name).read_bytes()).hexdigest()==h,name
    sys.path.insert(0,str(BOUND))
    # Import by name after removing this launcher name from consideration.
    spec=importlib.util.spec_from_file_location('run',BOUND/'run.py')
    module=importlib.util.module_from_spec(spec);sys.modules['run']=module;spec.loader.exec_module(module)
    import outcome_audit as audit
    receipt=json.loads((BOUND/'ROTOR_RECEIPT.json').read_text())
    source=json.loads((BOUND/'ROTOR_RECORDS.json').read_text())
    rows=[r for r in receipt['geometry_records']if r['coupling']=='polar']
    weights=[tuple(r['native_ticket_weights'])for r in rows]
    maps={w:audit.maps_from_records(table)for w,table in zip(weights,source['polar'])}
    sequence=weights+weights;reference=[(32,16,16)]*192
    outcomes=audit.outcomes(sequence,reference,maps)
    exact=json.loads((BOUND/'ROTOR_AUDIT.json').read_text())['two_turn_exact_check']
    assert abs(outcomes['log10_final_success_probability']-exact['log10_success'])<2e-8
    assert len(outcomes['stages'])==192
    assert all(len(r['conditional_probabilities_by_eraser_flag0_through383'])==384 for r in outcomes['stages'])
    result={'bound_manifest_sha256':MANIFEST,'additional_native_calls':0,
            'stages':192,'eraser_outcomes_per_stage':384,
            'log10_success':outcomes['log10_final_success_probability'],
            'probability_conservation_error':outcomes['probability_conservation_error'],
            'scope':'Original same96-step polar program twice; previously reported slow one-turn polar192 ledger retained unchanged',
            'G6_awarded':False,'G9_awarded':False}
    return result,outcomes

if __name__=='__main__':
    result,outcomes=calculate()
    if '--record'in sys.argv:
        for name,obj in [('RECEIPT.json',result),('TWO_TURN_OUTCOMES.json',outcomes)]:
            with(HERE/name).open('x')as f:f.write(json.dumps(obj,sort_keys=True,indent=2)+'\n')
    else:
        assert result==json.loads((HERE/'RECEIPT.json').read_text())
        assert outcomes==json.loads((HERE/'TWO_TURN_OUTCOMES.json').read_text())
    print(json.dumps(result,sort_keys=True))
