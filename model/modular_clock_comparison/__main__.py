"""Replay the source clock comparisons and exact modular frequency certificate."""
import argparse
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .comparison import report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args(argv)
    current = report()
    if args.write:
        args.write.write_text(canonical(current))
    else:
        compare_receipt(strict_json.load(args.check), current)
        print('Complete regional clock comparisons and exact nonresonance verified.')


if __name__ == '__main__':
    main()
