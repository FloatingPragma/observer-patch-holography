"""Reproduce the complete finite phase-spectroscopy receipt."""
import argparse
import json
from pathlib import Path
from .certificate import certificate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', nargs='?', const=str(Path(__file__).with_name('results.json')))
    args = parser.parse_args()
    text = json.dumps(certificate(), indent=2, sort_keys=True, allow_nan=False)+'\n'
    if args.write:
        Path(__file__).with_name('results.json').write_text(text)
    elif not Path(args.check).exists() or Path(args.check).read_text() != text:
        raise SystemExit('joint phase-spectroscopy receipt differs')
    print('PASS complete same-U phase readout, nonflat weight and energy-branch boundary')


if __name__ == '__main__':
    main()
