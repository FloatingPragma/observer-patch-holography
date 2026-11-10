"""Replay the source calculation and bind its full immutable input closure."""
import copy
import hashlib
import json
import math
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from .model import calculate, load_parents, PARENTS

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
    _,parent_pins=load_parents()
    validate_pins(parent_pins)
    seeds={str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(parent_pins)
    seeds.add('model/requirements.txt')
    names=import_closure(seeds)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def certificate():
    pins=input_hashes()
    value=copy.deepcopy(calculate())
    value['schema']='oph.interacting_block_propagation.v1'
    value['input_sha256']=pins
    value['immutable_parent_pins_validated']=PARENTS
    return value


def replay_differences(expected,actual):
    expected=json.loads(json.dumps(expected))
    actual=json.loads(json.dumps(actual))
    differences=[]
    def compare(left,right,path=''):
        if type(left)!=type(right):
            differences.append(path+' (type)');return
        if isinstance(left,dict):
            if set(left)!=set(right):
                differences.append(path+' (keys)');return
            for key in left:compare(left[key],right[key],path+'/'+key)
        elif isinstance(left,list):
            if len(left)!=len(right):
                differences.append(path+' (length)');return
            for i,(a,b) in enumerate(zip(left,right,strict=True)):compare(a,b,path+'/'+str(i))
        elif isinstance(left,float):
            if not math.isfinite(left) or not math.isfinite(right) or abs(left-right)>1e-10:
                differences.append(path+' (binary64 diagnostic)')
        elif left!=right:
            differences.append(path)
    compare(expected,actual)
    return differences


def write(path=None):
    value=certificate()
    path=Path(path) if path is not None else HERE/'results.json'
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    return value
