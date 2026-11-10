"""Recompute the full native neutral-odd spectral receipt."""
import argparse
import json
from pathlib import Path

from .certificate import certificate
from .model import generate_partition


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    for flag in ('write', 'check', 'generate-partition'):
        group.add_argument('--' + flag, action='store_true')
    args = parser.parse_args()
    if args.generate_partition:
        generate_partition()
    rendered = json.dumps(certificate(), indent=2, sort_keys=True) + '\n'
    path = Path(__file__).with_name('results.json')
    if args.check:
        if path.read_text() != rendered:
            raise SystemExit('FAIL stale neutral local spectral receipt')
    else:
        path.write_text(rendered)
    print('PASS exact full-native odd spectrum, centered charge and complete energy records')


if __name__ == '__main__':
    main()
