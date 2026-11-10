"""Replay the named coherent observer's complete finite certificate."""
import argparse
import json
from pathlib import Path
from .certificate import certificate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', nargs='?', const=str(Path(__file__).with_name('results.json')),
                       metavar='RECEIPT', help='recompute and compare with the optional receipt path')
    args = parser.parse_args()
    text = json.dumps(certificate(), indent=2, sort_keys=True, allow_nan=False) + '\n'
    path = Path(__file__).with_name('results.json')
    if args.write:
        temporary = path.with_suffix('.json.tmp')
        temporary.write_text(text)
        temporary.replace(path)
    elif not Path(args.check).exists() or Path(args.check).read_text() != text:
        raise SystemExit('coherent seam repair certificate differs')
    print('PASS coherent-seam-repair-v1: finite repair, reference, apparatus and static gauge response')


if __name__ == '__main__':
    main()
