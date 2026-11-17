#!/usr/bin/env python3
"""Offline custody check, exact record replay and independent finite check."""
import hashlib
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
for line in (HERE / 'MANIFEST.sha256').read_text().splitlines():
    digest, name = line.split('  ', 1)
    assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest, name
import run
import independent_check
result, records = run.run()
assert run.canonical(records) == (HERE / 'RECORDS.json').read_text()
saved = json.loads((HERE / 'RECEIPT.json').read_text())
# Exact records have byte custody. Numerical diagnostics use explicit tolerance
# across BLAS implementations; a bitwise numerical replay is not promised.
def close(a, b):
    if isinstance(a, float):
        assert abs(a - b) <= 2e-10
    elif isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a: close(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): close(x, y)
    else:
        assert a == b
close(result, saved)
independent = independent_check.run()
print(json.dumps({'diagnostic_checks': result['all_diagnostic_checks'],
                  'native_calls': result['native_endpoint_calls'],
                  'independent_branches': independent['checked_source_branches'],
                  'row_criterion_passed': result['row_criterion_passed'],
                  'native_candidate_window': result['dynamics']['native']['any_candidate_window']}, sort_keys=True))
