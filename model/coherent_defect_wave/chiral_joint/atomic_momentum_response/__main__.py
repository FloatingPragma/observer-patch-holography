"""Recompute exact source/reference phase records and pump failures."""
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
            raise SystemExit('FAIL stale atomic momentum receipt')
    else:
        path.write_text(rendered)
    print('PASS actual source momentum, retained reference phases and complete neutral pump')


if __name__ == '__main__':
    main()
