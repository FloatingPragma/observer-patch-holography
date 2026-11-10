"""Write or replay the complete periodic three-axis source receipt."""
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
    rendered = json.dumps(certificate(), indent=2, sort_keys=True) + '\n'
    if args.write:
        path.write_text(rendered)
    elif path.read_text() != rendered:
        raise SystemExit('FAIL stale periodic three-axis source receipt')
    print('PASS full periodic quantum source, three translations and physical phase records')


if __name__ == '__main__':
    main()
