"""Bounded offline replay of the finest complete common-state case and pair domain."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import core
import independent
import pair_source
import pair_check

begin=time.perf_counter();here=Path(__file__).resolve().parent
fresh=core.run(6);stored=json.loads(gzip.decompress((here/'k6/EVIDENCE.json.gz').read_bytes()))['base_execution']
for result in (fresh,stored):result.pop('seconds')
assert json.loads(json.dumps(fresh))==stored
pair=pair_source.run();saved=json.loads(gzip.decompress((here/'PAIR_SOURCE.json.gz').read_bytes()))
for result in (pair,saved):result.pop('seconds')
assert json.loads(json.dumps(pair))==saved
small_audit=independent.run(here,0);pair_audit=pair_check.run(here)
receipts=[]
for k in range(7):
    r=json.loads((here/f'k{k}/RECEIPT.json').read_text())
    assert hashlib.sha256((here/f'k{k}/EVIDENCE.json.gz').read_bytes()).hexdigest()==r['evidence_sha256']
    assert r['passed'];receipts.append(r)
result={'candidate':'A-native-common-refinement-family-v1','fresh_finest_cases':16,
        'all_bound_family_cases':sum(len(r['cases']) for r in receipts),'all_bound_family_passed':True,
        'finest_complete_payload_equal_except_runtime':True,'pair_domain_equal_except_runtime':True,
        'representative_independent_audit':small_audit,'pair_independent_audit':pair_audit,
        'seconds':time.perf_counter()-begin,'same_group':True,'external_G6':False,'physical_award':False}
print(json.dumps(result,indent=2,sort_keys=True));assert result['seconds']<60
