#!/usr/bin/env python3
"""Preserved failed-run replay, reporting the gate locals without changing it."""
import importlib.util
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('failed_return_primary',HERE/'probe.py')
p=importlib.util.module_from_spec(spec);sys.modules[spec.name]=p;spec.loader.exec_module(p)
try:
    p.main()
except AssertionError as error:
    trace=error.__traceback__;found=None
    while trace:
        values=trace.tb_frame.f_locals
        if trace.tb_frame.f_code.co_name=='run_case' and 'changes' in values:
            found={'charge_number':values['q'],'reference_differences':values['changes'],
                   'matrix_max_error':max(x for e in values['allerrors'].values() for x in e.values()),
                   'failed_requirement':'omit_loop >1e-6 and zero_current >1e-6',
                   'original_gate_passed':False}
        trace=trace.tb_next
    assert found is not None and found['charge_number']==0
    print(json.dumps(found,sort_keys=True,indent=2))
    raise SystemExit(0)
raise SystemExit('Expected preserved original gate to fail')
