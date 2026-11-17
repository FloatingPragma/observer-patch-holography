"""Replay complete joint energy and spatial-phase spectroscopy."""
import argparse
import json
from pathlib import Path
from model.particle_common_contract import contract, provenance
from model.particle_primary_endpoint_bounds import endpoint
from model.particle_primary_variational_readout import variational
from . import spectroscopy
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
    paths=set(spectroscopy.SOURCE_INPUTS)
    paths.update(endpoint.SOURCE_INPUTS)
    paths.update(variational.SOURCE_INPUTS)
    paths.update(bound_pins(contract.member(0)['source']['implementation_sha256']))
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.is_file() and p.suffix in ('.py','.md'))
    paths=provenance.definition_inputs(paths)
    pins=provenance.definition_hashes(paths)
    result=spectroscopy.certificate();result['input_sha256']=bound_pins(pins)
    if provenance.definition_hashes(paths)!=pins:
        raise ValueError('joint spectroscopy inputs changed during replay')
    return json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();result=build()
    if args.write:OUTPUT.write_text(result)
    elif not OUTPUT.is_file() or OUTPUT.read_text()!=result:raise ValueError('stale joint spectroscopy receipt')
    print('PASS: complete finite joint spectroscopy and scoped excitation intervals')


if __name__=='__main__':main()
