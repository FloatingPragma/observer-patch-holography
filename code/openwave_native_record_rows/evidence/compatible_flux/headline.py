#!/usr/bin/env python3
"""Custody and complete finite witness replay, normally a few seconds."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

P = Path(__file__).resolve().parent
for line in (P/'MANIFEST.sha256').read_text().splitlines():
    digest, name = line.split(maxsplit=1)
    assert hashlib.sha256((P/name).read_bytes()).hexdigest() == digest, name
names = ('RECEIPT.json', 'RECORDS.jsonl.gz', 'VERIFICATION.json')
before = {n: hashlib.sha256((P/'experiment'/n).read_bytes()).hexdigest() for n in names}
for script in ('run.py', 'verify.py'):
    subprocess.run([sys.executable, str(P/'experiment'/script)], check=True)
assert before == {n: hashlib.sha256((P/'experiment'/n).read_bytes()).hexdigest() for n in names}
print(json.dumps({'custody': True, 'full_witness_byte_replay': True,
                  'actual_native_calls': 98304, 'changed_law_events': 98304,
                  'full_Maxwell_execution': False, 'G6_pass': False}))
