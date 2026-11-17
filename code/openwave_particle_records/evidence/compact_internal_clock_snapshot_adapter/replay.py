#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PACKET=HERE.parent/'compact_internal_clock'
sys.path.insert(0,str(PACKET))
import run_serialized as original
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
snapshots=json.loads((HERE/'SNAPSHOTS.json').read_text())
registration=json.loads((HERE/'PREREG.json').read_text())
for path,digest in registration['files'].items():assert sha(ROOT/path)==digest,path
def snapshot_bound_pins():
    for name in ('SOURCE_PINS.json','PREREG.json','PREREG_AMENDMENT.json','PREREG_SERIALIZATION.json'):
        for path,digest in json.loads((PACKET/name).read_text())['files'].items():
            target=HERE/snapshots[path]['snapshot'] if path in snapshots else ROOT/path
            assert sha(target)==digest,path
            if path in snapshots:assert snapshots[path]['sha256']==digest
original.pins=snapshot_bound_pins
result,seconds=original.experiment(False)
receipt=json.loads((PACKET/'RECEIPT.json').read_text())
assert result==receipt['scientific']
for path,digest in receipt['artifacts'].items():assert sha(PACKET/path)==digest
print(json.dumps({'scientific_receipt_unchanged':True,'snapshot_bindings':snapshots,'elapsed_seconds':seconds},sort_keys=True))
assert seconds<60
