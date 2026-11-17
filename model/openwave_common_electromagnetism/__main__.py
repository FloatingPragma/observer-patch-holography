"""Replay the longitudinal common-source certificate."""
import argparse
import hashlib
import json
from pathlib import Path
from .source import certificate
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def pins():
    parent='model/openwave_common_mass_source/results.json'
    seeds=set(json.loads((ROOT/parent).read_text())['source_sha256'])|{parent}
    seeds|={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md','.tex')}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();before=pins();data=certificate();assert pins()==before
    data['source_sha256']=before;path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')
    else:assert equivalent(json.loads(path.read_text()),data)
    print(data['verdict'])


if __name__=='__main__':main()
