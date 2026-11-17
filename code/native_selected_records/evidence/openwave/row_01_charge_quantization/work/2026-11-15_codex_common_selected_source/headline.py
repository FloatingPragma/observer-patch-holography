"""Recompute the common-source experiment and independently audit retained outcomes."""
import json
import time
from common import run
from audit import run as verify
from pathlib import Path

start=time.perf_counter();out=run();audit=verify(Path(__file__).resolve().parent)
summary={'candidate':out['candidate'],'gates':out['gates'],'fresh_execution_passed':out['passed'],
         'same_group_audit':audit,'literal_native_calls':out['literal_native_calls'],
         'native_visibilities':[c['visibility'] for c in out['cases'] if c['law']=='native'],
         'frequency_max_error':max(abs(c['reader']['frequency']-out['expected_frequency']) for c in out['cases']),
         'seconds':time.perf_counter()-start,'external_G6':False,'physical_row_award':False}
print(json.dumps(summary,indent=2,sort_keys=True))
assert summary['seconds']<60 and out['passed'] and audit['passed']
