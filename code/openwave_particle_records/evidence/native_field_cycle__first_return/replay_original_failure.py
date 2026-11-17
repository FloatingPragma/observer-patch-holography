#!/usr/bin/env python3
"""Reproduce the exact archived original failed gate after the amendment."""
import hashlib
import json
from pathlib import Path
import sys
import types
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
ARCHIVE=HERE/'attempts/01_neutral_control'
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def check():
    for name,h in json.loads((HERE/'ORIGINAL_REPLAY_PREREG.json').read_text())['files'].items():
        assert digest(HERE/name)==h,name
    for name,h in json.loads((ARCHIVE/'PREREG.json').read_text())['files'].items():
        assert digest(ARCHIVE/name)==h,name
    for name,h in json.loads((ARCHIVE/'SOURCE_PINS.json').read_text()).items():
        assert digest(ROOT/name)==h,name

if __name__=='__main__':
    if len(sys.argv)!=1:raise SystemExit('Read-only original-failure replay takes no options')
    check()
    primary=types.ModuleType('archived_first_return_failure')
    # Original __file__ preserves provider resolution; bytes come only from the archive.
    primary.__file__=str(HERE/'probe.py');sys.modules[primary.__name__]=primary
    exec(compile((ARCHIVE/'probe.py').read_text(),str(ARCHIVE/'probe.py'),'exec'),vars(primary))
    primary.check_pins=check
    found=None
    try:primary.main()
    except AssertionError as error:
        trace=error.__traceback__
        while trace:
            v=trace.tb_frame.f_locals
            if trace.tb_frame.f_code.co_name=='run_case' and 'changes' in v:
                found={'charge_number':v['q'],'reference_differences':v['changes'],
                       'matrix_max_error':max(x for e in v['allerrors'].values() for x in e.values()),
                       'failed_requirement':'omit_loop >1e-6 and zero_current >1e-6',
                       'original_gate_passed':False}
            trace=trace.tb_next
    assert found==json.loads((ARCHIVE/'FAILURE_GATE.log').read_text())
    check()
    print(json.dumps({'original_failure_reproduced':True,'diagnostic':found},sort_keys=True,indent=2))
