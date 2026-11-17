#!/usr/bin/env python3
"""Verify sealed package, corrected execution freeze and source hashes."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

manifest=json.loads((HERE/'MANIFEST.json').read_text())
for name,value in manifest['files'].items():
    assert sha(HERE/name)==value,name
pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
for name,value in pins.items():assert sha(ROOT/name)==value,name
assert pins==json.loads((HERE/'FROZEN_SOURCE_PINS.json').read_text())
freeze=json.loads((HERE/'FROZEN_EXECUTION.json').read_text())
for filename,key in [('verify.py','code_sha256'),('preregistration.json','preregistration_sha256'),
                     ('FROZEN_SOURCE_PINS.json','source_pins_sha256')]:
    assert sha(HERE/filename)==freeze[key],filename
invalid=HERE/'attempts/01_invalid'
receipt=json.loads((invalid/'receipt.json').read_text())
for name,value in receipt['files'].items():assert sha(invalid/name)==value,name
print(json.dumps({'all_pass':True,'sealed_files':len(manifest['files']),
                  'source_inputs':len(pins),'manifest_sha256':sha(HERE/'MANIFEST.json')},indent=2))
