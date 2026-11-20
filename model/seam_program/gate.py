"""The eight-slot gate that the seam program compiles, and the archive helpers of its records."""
import json
from pathlib import Path

import numpy as np

from tools.strict_json import load

ROOT = Path(__file__).resolve().parents[2]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def multipliers(angles):
    return np.exp(-1j*np.repeat(angles, 2, axis=-1)*np.tile([1,-1],4))


def reflect(state, phase=1):
    require(state.ndim>=1 and state.shape[-1]==8,'the routing coin has eight slots')
    return (1+phase)*state.mean(axis=-1,keepdims=True)-state


def real_step(state, phase=1):
    """Literal four-index update of the target gate on the program's index grid."""
    require(state.ndim==5 and all(state.shape[:-1]),'four nonempty address coordinates required')
    coin = reflect(state,phase)
    out = np.empty_like(coin)
    for slot in range(8):
        out[...,slot^1] = np.roll(coin[...,slot], 1 if slot%2==0 else -1, axis=slot//2)
    return out


def counters(probabilities, shots):
    p = np.asarray(probabilities)
    require(isinstance(shots,(int,np.integer)) and not isinstance(shots,(bool,np.bool_))
            and 0<shots<=np.iinfo(np.int32).max,'exposures must fit a positive detector count')
    require(p.ndim==1 and p.size>0 and p.dtype.kind in 'fiu' and
            np.isfinite(p).all() and p.min()>=0 and abs(p.sum()-1)<2e-10,
            'incomplete detector law')
    edges = np.floor(shots*np.cumsum(p)+0.5).astype(np.int64)
    result = np.diff(edges,prepend=0)
    require(result.min()>=0 and result.sum()==shots,'lost a terminal outcome')
    return result.astype(np.int32)


def serialize(value):
    if isinstance(value,dict):
        return {k:serialize(v) for k,v in value.items()}
    if isinstance(value,list):
        return [serialize(v) for v in value]
    if isinstance(value,(float,np.floating)):
        return float(format(value,'.9g'))
    return value


def receipt(path,value,record):
    content = json.dumps(serialize(value),indent=1,sort_keys=True,allow_nan=False)+'\n'
    if record:
        with path.open('x',encoding='utf8',newline='\n') as out:
            out.write(content)
    stored = json.dumps(load(path),indent=1,sort_keys=True,allow_nan=False)+'\n'
    require(stored==content,f'changed result: {path.name}')


def archive(path,arrays,record):
    if record:
        with path.open('xb') as out:
            np.savez_compressed(out,**arrays)
    with np.load(path,allow_pickle=False) as stored:
        require(sorted(stored.files)==sorted(arrays),'changed detector menu')
        for k,a in arrays.items():
            require(stored[k].dtype==a.dtype and np.array_equal(stored[k],a),f'changed detector array: {k}')
        return dict(stored)
