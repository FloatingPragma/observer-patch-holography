"""Recompute the complete result and its repository-local source closure."""
import copy
import hashlib
import json
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from .model import calculate,source,PARENT,PARENT_SHA
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent


def validate_pins(pins):
    for name,expected in pins.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('source custody requires repository-relative paths')
        path=ROOT/relative
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('source custody mismatch: '+name)


def input_hashes():
    parent,_=source()
    inherited=dict(bound_pins(parent['input_sha256']));inherited[PARENT]=PARENT_SHA
    validate_pins(inherited)
    seeds={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(inherited);seeds.add('model/requirements.txt')
    names=import_closure(seeds)
    return bound_pins({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)})


def certificate():
    pins=input_hashes();value=copy.deepcopy(calculate())
    value['schema']='oph.interacting_block_proca.v1'
    value['input_sha256']=bound_pins(pins)
    value['immutable_parent_pins_validated']=[PARENT]
    return value


def write(path=None):
    value=certificate();path=Path(path) if path is not None else HERE/'results.json'
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    return value
