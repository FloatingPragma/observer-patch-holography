#!/usr/bin/env python3
"""Offline complete component replay; frozen integer rotor inputs are portable."""
import hashlib
import json
import sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
for line in(HERE/'MANIFEST.sha256').read_text().splitlines():
    digest,name=line.split('  ',1)
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
for filename in('PREREG.json','BALANCED_PREREG.json','ROTOR_PREREG.json','ROTOR_AUDIT_PREREG.json'):
    obj=json.loads((HERE/filename).read_text())
    for name,h in obj['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,(filename,name)
import run
import independent_check
import outcome_audit
import rotor_audit

def close(a,b):
    if isinstance(a,float):assert abs(a-b)<=2e-10,(a,b)
    elif isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:close(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):close(x,y)
    else:assert a==b,(a,b)

receipt,records=run.execute()
assert run.canonical(records)==(HERE/'RECORDS.json').read_text()
close(receipt,json.loads((HERE/'RECEIPT.json').read_text()))
check=independent_check.check()
close(check,json.loads((HERE/'INDEPENDENT.json').read_text()))
balanced,outcomes,corrected=outcome_audit.execute()
close(balanced,json.loads((HERE/'BALANCED_RECEIPT.json').read_text()))
close(outcomes,json.loads((HERE/'ALL_OUTCOMES.json').read_text()))
close(corrected,json.loads((HERE/'CORRECTED_RECEIPT.json').read_text()))
rotor,outcomes=rotor_audit.execute(reexecute_source=True)
calls=rotor['source_calls_reexecuted'];rotor['source_calls_reexecuted']=0
close(rotor,json.loads((HERE/'ROTOR_AUDIT.json').read_text()))
close(outcomes,json.loads((HERE/'ROTOR_OUTCOMES.json').read_text()))
# Regenerate the original rotor products using its frozen integer inputs.
r=json.loads((HERE/'ROTOR_RECEIPT.json').read_text())
raw=json.loads((HERE/'ROTOR_RECORDS.json').read_text())
for coupling in('polar','apolar'):
    seq=[outcome_audit.maps_from_records(t)[0]for t in raw[coupling]]
    reference=seq[-1]
    for dwell in(1,2):
        got=run.interference([p for p in seq for _ in range(dwell)],reference)
        wanted=r['rotor_products'][coupling][str(dwell)]
        for key,value in got.items():close(value,wanted[key])
# Recheck the retained matched-histogram refinement, including failed dwell4.
base=[tuple(w)for w in receipt['weight_loop']]
base_maps={w:outcome_audit.maps_from_records(t)[0]for w,t in zip(base,records['source']['native']['loop'])}
for dwell in(4,8):
    wa=[w for w in base for _ in range(dwell)]
    wb=[w for w in base+base for _ in range(dwell//2)]
    got=outcome_audit.path_compare([base_maps[w]for w in wa],[base_maps[w]for w in wb])
    close(got,r['refinement'][str(dwell)]['comparison'])
    close(outcome_audit.exact_compare(wa,wb),r['refinement'][str(dwell)]['exact_check'])
print(json.dumps({'all_component_checks':True,'original_native_calls_reexecuted':receipt['native_calls'],
                  'rotor_native_calls_reexecuted':calls,'exact_source_branches_checked':check['checked_body_branches'],
                  'polar_visibility':r['rotor_products']['polar']['1']['relative_interference_visibility'],
                  'two_turn_visibility':rotor['two_turns']['relative_visibility'],
                  'physical_rotation_derived':False,'row_criterion_passed':False},sort_keys=True))
