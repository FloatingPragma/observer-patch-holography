"""Strict offline replay of final candidate plus separate source/selector checks."""
import contextlib
import gzip
import io
import json
from pathlib import Path
import time
from selectors_v2 import run
from audit import run as verify

start=time.perf_counter();here=Path(__file__).resolve().parent
fresh=run();recorded=json.loads(gzip.decompress((here/'v2_run/EVIDENCE.json.gz').read_bytes()))
for result in (fresh,recorded):
    result.pop('seconds');result['base_execution'].pop('seconds')
fresh=json.loads(json.dumps(fresh))  # canonicalize JSON object keys before comparing
assert json.dumps(fresh,sort_keys=True,separators=(',',':'))==json.dumps(recorded,sort_keys=True,separators=(',',':'))
audit=verify(here)
with contextlib.redirect_stdout(io.StringIO()):
    import selector_audit
summary={'candidate':fresh['candidate'],'strict_replay_excluding_runtime_only':True,'passed':fresh['passed'],
         'complete_source_audit':audit,'selector_audit_passed':selector_audit.out['passed'],
         'selector_first_failure_stages':sum(len(c['selector_first_failure_stages']) for c in selector_audit.out['all_selector_first_failure_ledgers']),
         'literal_native_calls':fresh['literal_all_calls'],'same_object_cases':len(fresh['cases']),
         'seconds':time.perf_counter()-start,'same_group':True,'external_G6':False,'physical_award':False}
print(json.dumps(summary,indent=2,sort_keys=True))
assert summary['passed'] and summary['seconds']<60
