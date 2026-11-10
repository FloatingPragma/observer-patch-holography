"""Write or independently check the complete physical-arrival receipt."""
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
    path = Path(__file__).with_name('results.json')
    value = json.dumps(certificate(), indent=2, sort_keys=True) + '\n'
    if args.write:
        path.write_text(value)
    elif path.read_text() != value:
        raise SystemExit('FAIL neutral mobility receipt is stale')
    print('PASS full interacting source, physical odd arrival and complete outcome bounds')


if __name__ == '__main__':
    main()
