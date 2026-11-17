"""Replay source-bound excitation transport and conditional mass readouts."""
import argparse
import hashlib
import json
from pathlib import Path

from model.particle_common_contract import contract
from . import readout
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def build():
    paths = set(readout.SOURCE_INPUTS)
    paths.update(bound_pins(contract.member(0)['source']['implementation_sha256']))
    paths.update(str(path.relative_to(ROOT)) for path in HERE.iterdir()
                 if path.is_file() and path.suffix in ('.py', '.md'))
    pins = bound_pins({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)})
    result = readout.certificate()
    result['input_sha256'] = bound_pins(pins)
    if any(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != value
           for name, value in pins.items()):
        raise ValueError('inputs changed during primary mass readout replay')
    return json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.write:
        OUTPUT.write_text(result)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != result:
        raise ValueError('stale primary mass readout receipt')
    print('PASS: original excitation transport, conditional rest intervals and declared units')


if __name__ == '__main__':
    main()
