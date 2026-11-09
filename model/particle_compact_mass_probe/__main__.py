"""Replay compact four-outcome readout and its explicit resource boundary."""
import argparse
import hashlib
import json
from pathlib import Path

from model.particle_common_contract import contract
from . import compact

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def build():
    paths = set(compact.SOURCE_INPUTS)
    paths.update(contract.member(0)['source']['implementation_sha256'])
    paths.update(str(path.relative_to(ROOT)) for path in HERE.iterdir()
                 if path.is_file() and path.suffix in ('.py', '.md'))
    pins = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}
    result = compact.certificate()
    result['input_sha256'] = pins
    if any(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != value for name, value in pins.items()):
        raise ValueError('compact readout inputs changed during replay')
    return json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args(); result = build()
    if args.write:
        OUTPUT.write_text(result)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != result:
        raise ValueError('stale compact primary readout receipt')
    print('PASS: same original response with four complete outcomes; distinct failures')


if __name__ == '__main__':
    main()
