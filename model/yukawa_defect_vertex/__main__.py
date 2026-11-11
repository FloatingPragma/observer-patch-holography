"""Write or replay the complete finite source receipt."""
import argparse
import json
from pathlib import Path

from .vertex import certificate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    path = Path(__file__).with_name('results.json')
    fresh = certificate()
    if args.check:
        if json.loads(path.read_text()) != fresh:
            raise ArithmeticError('the complete source receipt differs on replay')
        print('Complete scalar intervention and vertex receipt reproduced.')
    elif args.write:
        path.write_text(json.dumps(fresh, indent=2, sort_keys=True)+'\n')
        print('Wrote complete scalar intervention and vertex receipt.')
    else:
        print(json.dumps(fresh, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
