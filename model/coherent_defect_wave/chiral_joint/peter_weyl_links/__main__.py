"""Recompute every finite source and observer-record certificate."""
import argparse
import json
from pathlib import Path

from .certificate import certificate


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    rendered = json.dumps(certificate(), indent=2, sort_keys=True) + '\n'
    path = Path(__file__).with_name('results.json')
    if args.check:
        if path.read_text() != rendered:
            raise SystemExit('FAIL stale Peter-Weyl receipt')
    else:
        path.write_text(rendered)
    print('PASS complete Peter-Weyl links and finite physical neutral transfer')


if __name__ == '__main__':
    main()
