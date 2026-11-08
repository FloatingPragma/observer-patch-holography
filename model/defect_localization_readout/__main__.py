"""Write or replay the localized source-excitation receipt at the given path."""
import argparse
from pathlib import Path

from . import packet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', nargs='?', const=packet.OUTPUT, type=Path)
    action.add_argument('--check', nargs='?', const=packet.OUTPUT, type=Path)
    args = parser.parse_args()
    if args.check is not None:
        packet.verify(args.check)
    else:
        args.write.write_text(packet.canonical(packet.build()))
    print('PASS defect_localization_readout')


if __name__ == '__main__':
    main()
