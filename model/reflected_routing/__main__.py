"""Reproduce registered calibration, sealed predictions and subsequently released records."""
import argparse
import hashlib
import json
import secrets
from pathlib import Path

import numpy as np

from tools.strict_json import load
from . import clock,race
from .frontier import exact_tails
from .kernel import require

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECORDS = HERE/'records'


def source_pins():
    paths = sorted(HERE.glob('*.py'))
    paths = [p for p in paths if not p.name.startswith('test_')]
    paths.extend([HERE/'protocol.json',ROOT/'tools/strict_json.py',ROOT/'model/requirements.txt'])
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in paths}


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


def run(stage,root,record):
    require(source_pins()==load(HERE/'sources.json'),'construction or protocol changed')
    spec = load(HERE/'protocol.json')
    if record and stage=='seal':
        root.mkdir(exist_ok=False,parents=True)
    writes = record and stage=='seal'
    x = archive(root/'timing.npz',clock.acquire(spec,'timing'),writes)
    rest = archive(root/'rest.npz',clock.acquire(spec,'rest'),writes)
    forecast = clock.predict(spec,x,rest)
    targets = archive(root/'finite_forecast.npz',clock.finite_forecasts(spec,rest),writes)
    channels = race.forecast(spec['race'],forecast[f"{max(spec['sizes'])}_out"]['reference'])
    receipt(root/'forecast.json',dict(clock=forecast,channels=channels),writes)
    if stage=='seal':
        return
    if record:
        with (root/'release.json').open('x',encoding='utf8') as out:
            json.dump({'sender_entropy':secrets.token_hex(32)},out)
            out.write('\n')
    seed = bytes.fromhex(load(root/'release.json')['sender_entropy'])
    require(len(seed)==32,'invalid trial entropy')
    y = archive(root/'moving.npz',clock.acquire(spec,'moving'),record)
    actual = clock.open_pulse(spec,forecast,targets,y,rest)
    raw,result = race.acquire(spec['race'],channels,seed)
    archive(root/'channel.npz',raw,record)
    tails = exact_tails(spec['tail_waits'],spec['tail_speeds'])
    receipt(root/'result.json',dict(clock=actual,channels=result,tails=tails),record)
    print(json.dumps({k:v['decision'] for k,v in actual.items()}))
    print(json.dumps({k:[v['accepted'],v['ahead']] for k,v in result.items()}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['seal','open'],default='open')
    parser.add_argument('--record',type=Path,help='new receipt directory (seal first, then open)')
    args = parser.parse_args()
    run(args.stage,args.record or RECORDS,args.record is not None)


if __name__=='__main__':
    main()
