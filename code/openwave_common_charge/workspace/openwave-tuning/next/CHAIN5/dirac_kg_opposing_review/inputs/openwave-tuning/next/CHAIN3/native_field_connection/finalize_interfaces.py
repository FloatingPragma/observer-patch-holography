"""Refresh small interface contracts without repeating completed field checks."""
from fractions import Fraction as F
import json
from providers import HERE, PINS, verify
import run
import native_oscillators as words

verify()
result=json.loads((HERE/'RESULT.json').read_text())
assert result['status']=='PASS_SCOPED_CONSTRUCTIVE_TENSOR_FIELD_INTERFACE'
result['complete_native_pointer']=run.pointer_control()
print('rational READY pointer all8 columns and both C/E fine effects verified',flush=True)
static=words.static_holding_contract(F(8))
run.write('STATIC_HOLDING.json',static)
for row in result['source_contracts']:
    run.complete_edge_bound(row)
result['static_holding_contract']='STATIC_HOLDING.json'
result['producer_input_files']=len(PINS)
run.write('RESULT.json',result)
run.freeze_outputs()
print(json.dumps({'status':result['status'],'static_caps':[static['coexact_cap'],static['harmonic_cap']],
                  'edge_TV_bounds':[row['complete_edge_Fourier_record_TV_bound_upper'] for row in result['source_contracts']]}),flush=True)
