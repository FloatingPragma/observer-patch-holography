#!/usr/bin/env python3
"""Read-only primary replay using the two original documentation snapshots."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def check_adapter():
    for name,h in json.loads((HERE/'REPLAY_PREREG.json').read_text())['files'].items():
        assert digest(HERE/name)==h,name

def check_pins():
    check_adapter()
    for name,h in json.loads((HERE/'PREREG.json').read_text())['files'].items():
        assert digest(HERE/name)==h,name
    mapping=json.loads((HERE/'FROZEN_INPUT_MAP.json').read_text())
    source=json.loads((HERE/'SOURCE_PINS.json').read_text())
    assert set(mapping)=={'openwave/row_05_magnetic_moment/THEOREMS.md','openwave/row_06_spin_half/THEOREMS.md'}
    for name,h in source.items():
        if name in mapping:
            entry=mapping[name];assert entry['sha256']==h,name
            path=HERE/entry['snapshot']
        else:path=ROOT/name
        assert digest(path)==h,name

if __name__=='__main__':
    if len(sys.argv)!=1:raise SystemExit('Read-only replay takes no options')
    check_pins()
    spec=importlib.util.spec_from_file_location('frozen_field_cycle_primary',HERE/'probe.py')
    primary=importlib.util.module_from_spec(spec);sys.modules[spec.name]=primary;spec.loader.exec_module(primary)
    primary.check_pins=check_pins
    primary.main()
    check_pins()
