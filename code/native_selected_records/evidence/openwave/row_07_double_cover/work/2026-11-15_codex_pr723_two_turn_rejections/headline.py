#!/usr/bin/env python3
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
for line in(HERE/'MANIFEST.sha256').read_text().splitlines():
    h,n=line.split('  ',1);assert hashlib.sha256((HERE/n).read_bytes()).hexdigest()==h,n
spec=importlib.util.spec_from_file_location('supplement',HERE/'run.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
result,outcomes=module.calculate()
def close(a,b):
    if isinstance(a,float):assert abs(a-b)<=2e-10
    elif isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:close(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):close(x,y)
    else:assert a==b
close(result,json.loads((HERE/'RECEIPT.json').read_text()))
close(outcomes,json.loads((HERE/'TWO_TURN_OUTCOMES.json').read_text()))
print(json.dumps(result,sort_keys=True))
