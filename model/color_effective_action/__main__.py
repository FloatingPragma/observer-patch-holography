"""Build or replay the complete finite color effective action."""
import argparse
from pathlib import Path

from . import packet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        packet.verify(args.check)
        print('finite color effective action verified')
    else:
        args.write.write_text(packet.canonical(packet.build()), encoding='utf-8', newline='\n')
        print('finite color effective action written')


if __name__ == '__main__':
    main()
