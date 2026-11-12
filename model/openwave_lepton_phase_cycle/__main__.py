"""Reproduce the exact phase-cycle certificate and complete mass records."""
import argparse
import hashlib
import json
from pathlib import Path
from model.lepton_mass_compiler.packet import bindings as original_bindings
from .compiler import run

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def produce():
    pins=original_bindings()
    pins.update({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in HERE.iterdir() if p.suffix in ('.py','.md')})
    result=run();result['source_sha256']=pins
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    choice=ap.add_mutually_exclusive_group(required=True)
    choice.add_argument('--write',action='store_true');choice.add_argument('--check',action='store_true')
    args=ap.parse_args();value=produce();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    elif json.loads(path.read_text())!=value:raise SystemExit('phase-cycle source result differs')
    print(value['verdict'])


if __name__=='__main__':main()
