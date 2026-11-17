#!/usr/bin/env python3
from run_amended import experiment, HERE
import json
r,seconds=experiment()
assert r==json.loads((HERE/'RECEIPT.json').read_text())['scientific']
print(json.dumps({'native_calls':r['source']['native_calls'],'theta':r['clock']['theta'],
                  'position_envelope':r['localization']['position_envelope'],
                  'mass_clock_gates':r['energy']['mass_clock_gates'],
                  'autonomous_onset':r['onset']['autonomous_onset_established'],
                  'elapsed_seconds':seconds},sort_keys=True))
assert seconds<60
