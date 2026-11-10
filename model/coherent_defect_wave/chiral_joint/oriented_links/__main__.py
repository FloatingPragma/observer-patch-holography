"""Regenerate or check the exact finite hardware receipt."""
import argparse
import json
from pathlib import Path

from .certificate import certificate


def main():
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = Path(__file__).with_name('results.json')
    rendered = json.dumps(certificate(), indent=2, sort_keys=True) + '\n'
    if args.write:
        path.write_text(rendered)
    elif path.read_text() != rendered:
        raise SystemExit('FAIL oriented-link receipt is stale')
    print('PASS full native-current links, contractive paths and proper cubic rotations')


if __name__ == '__main__':
    main()
