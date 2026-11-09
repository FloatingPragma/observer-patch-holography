"""Replay populated primary finite responses with all failure records."""
import argparse
import json
from pathlib import Path
from model.particle_common_contract import contract,provenance
from . import model

HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def build():
    paths=set(model.SOURCE_INPUTS)
    paths.update(str(p.relative_to(contract.ROOT)) for p in HERE.iterdir()
                 if p.is_file() and p.suffix in ('.py','.md'))
    paths=provenance.definition_inputs(paths);pins=provenance.definition_hashes(paths)
    result=model.certificate();result['input_sha256']=pins
    if provenance.definition_hashes(paths)!=pins:raise ValueError('refinement inputs changed during replay')
    return json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();value=build()
    if args.write:OUTPUT.write_text(value)
    elif not OUTPUT.is_file() or OUTPUT.read_text()!=value:raise ValueError('stale primary refinement response')
    print('PASS: original populated finite primary responses at every capacity')


if __name__=='__main__':main()
