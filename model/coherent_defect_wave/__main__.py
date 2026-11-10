"""Regenerate or replay the exact finite interacting defect-wave receipt."""
import argparse
import json
from pathlib import Path
from .certificate import certificate

OUTPUT = Path(__file__).with_name('results.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    rendered = json.dumps(certificate(), sort_keys=True, indent=2)+'\n'
    if args.write:
        OUTPUT.write_text(rendered)
    elif not args.check.exists() or args.check.read_text() != rendered:
        raise ValueError('stale joint defect-wave receipt')
    print('PASS exact joint defect-wave interaction, complete readout and same-clock reference')


if __name__ == '__main__':
    main()
