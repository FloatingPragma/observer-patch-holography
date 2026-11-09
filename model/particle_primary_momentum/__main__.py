"""Replay complete spatial Fourier records on the selected source."""
import argparse
import hashlib
import json
from pathlib import Path
from model.particle_common_contract import contract
from . import momentum

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def build():
    paths=set(momentum.SOURCE_INPUTS)
    paths.update(contract.member(0)['source']['implementation_sha256'])
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.is_file() and p.suffix in ('.py','.md'))
    pins={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}
    result=momentum.certificate(); result['input_sha256']=pins
    if any(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=value for name,value in pins.items()):
        raise ValueError('momentum prerequisite inputs changed during replay')
    return json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=build()
    if args.write:OUTPUT.write_text(result)
    elif not OUTPUT.is_file() or OUTPUT.read_text()!=result:raise ValueError('stale spatial Fourier receipt')
    print('PASS: complete spatial Fourier records; conserved primary momentum remains unassigned')


if __name__=='__main__':main()
