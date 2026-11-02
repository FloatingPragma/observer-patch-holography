"""Replay the exact finite controls and source-bound continuum-energy report."""
import argparse
import json

from . import packet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = packet.HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(packet.build(), indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        packet.validate(json.loads(path.read_text()))
    print('PASS spin_canonical_energy '+('write' if args.write else 'check'))


if __name__ == '__main__':
    main()
