"""Replay spatial smearing and actual finite-source bulk transfer bounds."""
import argparse
import hashlib
import json
from pathlib import Path

from model.particle_common_contract import contract
from . import spatial

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def build():
    paths = set(spatial.SOURCE_INPUTS)
    paths.update(contract.member(0)['source']['implementation_sha256'])
    paths.update(str(path.relative_to(ROOT)) for path in HERE.iterdir()
                 if path.is_file() and path.suffix in ('.py', '.md'))
    pins = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)}
    result = spatial.certificate()
    result['input_sha256'] = pins
    if any(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != value
           for name, value in pins.items()):
        raise ValueError('inputs changed during spatial algebra replay')
    return json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.write:
        OUTPUT.write_text(result)
    elif args.check:
        if not OUTPUT.is_file() or OUTPUT.read_text() != result:
            raise ValueError('stale spatial algebra receipt')
    else:
        print(result, end='')
    if args.write or args.check:
        print('PASS: regular spatial translations and same-source bulk word comparison')


if __name__ == '__main__':
    main()
