"""Replay the ordered-spin relaxed-response certificate."""
import argparse
from . import packet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', action='store_true')
    actions.add_argument('--check', nargs='?', const=str(packet.OUTPUT), metavar='PATH')
    args = parser.parse_args()
    if args.write:
        packet.OUTPUT.write_text(packet.canonical(packet.build()))
    else:
        packet.verify(args.check)
    print('PASS spin_relaxed_response '+('write' if args.write else 'check'))


if __name__ == '__main__':
    main()
