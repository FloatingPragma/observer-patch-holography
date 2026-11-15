"""Recompute the compiled native interface and exact binary error bounds."""
import argparse
from pathlib import Path

from tools import strict_json
from ..common_source.infinite_history import canonical
from ..common_source.output_certificate import compare_receipt
from .construction import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    paths = parser.add_mutually_exclusive_group()
    paths.add_argument('--write', type=Path)
    paths.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Compiled native interface receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
