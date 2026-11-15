"""Replay the finite native quantum transfer and coherent measurements."""
from __future__ import annotations

import argparse
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .construction import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if result['all_pass'] is not True:
        raise ArithmeticError('native quantum history failed')
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native quantum history receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
