#!/usr/bin/env python3
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
from probe import run
r=run(); data=(json.dumps(r,sort_keys=True,indent=2)+'\n').encode()
same=data==(Path(__file__).with_name('RECEIPT.json')).read_bytes()
print(json.dumps({'native_calls':r['native_calls'],'green_target':r['green_target'],
                 'n512':r['native_moments']['512'],'all_component_checks':all(r['checks'].values()),
                 'byte_identical':same,'full_physical_rows_passed':False},sort_keys=True))
sys.exit(0 if same and all(r['checks'].values()) else 1)
