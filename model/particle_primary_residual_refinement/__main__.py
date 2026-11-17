"""Replay exact local moments and the full physical trial residual."""
import argparse
import json
from pathlib import Path

from model.particle_common_contract import contract,provenance
from model.particle_primary_variational_readout import variational
from . import model
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def build():
    paths=set(model.SOURCE_INPUTS)|set(variational.SOURCE_INPUTS)
    paths.update(bound_pins(contract.member(0)['source']['implementation_sha256']))
    paths.update(str(path.relative_to(ROOT)) for path in HERE.iterdir()
                 if path.is_file() and path.suffix in ('.py','.md'))
    paths=provenance.definition_inputs(paths);pins=provenance.definition_hashes(paths)
    result=model.certificate();result['input_sha256']=bound_pins(pins)
    if provenance.definition_hashes(paths)!=pins:raise ValueError('residual source changed during replay')
    return json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();result=build()
    if args.write:OUTPUT.write_text(result)
    elif not OUTPUT.is_file() or OUTPUT.read_text()!=result:raise ValueError('stale exact residual receipt')
    print('PASS: exact parent local moments and retained full-H residual')


if __name__=='__main__':main()
