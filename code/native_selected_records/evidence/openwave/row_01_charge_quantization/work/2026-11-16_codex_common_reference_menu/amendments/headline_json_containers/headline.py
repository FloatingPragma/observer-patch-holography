"""Bounded fresh primary-source replay, exact-path correction and sealed-grid custody."""
from pathlib import Path
import gzip
import hashlib
import json
import tempfile
import time
import experiment
from exact_path_amendment import correct
from reader import read

HERE=Path(__file__).resolve().parent
start=time.monotonic()
registration=json.loads((HERE/'EXACT_PATH_COUNT_AMENDMENT.json').read_text())
for name,digest in registration['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
corrected=json.loads((HERE/'EXACT_PATH_CORRECTION.json').read_text())
assert corrected['passed'] and corrected['cases']==224

def stable(value):
    if isinstance(value,dict):return {k:stable(v) for k,v in value.items() if k!='seconds'}
    if isinstance(value,list):return [stable(x) for x in value]
    return value

fresh=[]
with tempfile.TemporaryDirectory(prefix='common-reference-menu-') as folder:
    for mode in ('matched','stationary'):
        output=Path(folder)/mode
        receipt=experiment.execute(0,mode,output)
        actual=json.loads(gzip.decompress((output/'EVIDENCE.json.gz').read_bytes()))
        filed=json.loads(gzip.decompress((HERE/f'k0_{mode}/EVIDENCE.json.gz').read_bytes()))
        assert stable(actual)==stable(filed),'complete primary payload changed'
        repaired=correct(actual)
        saved=next(m for m in corrected['members'] if m['k']==0 and m['reference']==mode)
        expected={k:v for k,v in saved.items() if k!='primary_evidence_sha256'}
        assert repaired==expected and repaired['passed'],'source-derived correction changed'
        fresh.append({'mode':mode,'cases':16,'complete_primary_equal_except_runtime':True,
                      'original_primary_passed':receipt['passed'],'corrected_physical_law_passed':True})
finest=0
for mode in ('matched','stationary'):
    receipt=json.loads((HERE/f'k6_{mode}/RECEIPT.json').read_text())
    for case in receipt['cases']:
        records=case['pair_records']
        result=read(records['1'],records['8'],eta=.01,probability_error=2e-5)
        assert result==case['pair_reader'] and result['unique_component']
        assert experiment.SCALE*result['frequency_components'][0][1]<experiment.OLD_PAIR
        finest+=1
print(json.dumps({'candidate':experiment.CANDIDATE,'fresh_primary_cases':32,'fresh':fresh,
    'all_bound_cases':224,'all_corrected_cases_passed':True,'finest_record_only_pair_readouts':finest,
    'original_failed_raw_positivity_gates_preserved':True,'same_group':True,'external_G6':False,
    'seconds':time.monotonic()-start},indent=2))
