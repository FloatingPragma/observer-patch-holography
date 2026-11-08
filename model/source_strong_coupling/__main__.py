"""Write or replay the native coupling receipt at a requested path."""
import argparse
from pathlib import Path
from . import packet


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write',nargs='?',const=packet.OUTPUT,type=Path)
    actions.add_argument('--check',nargs='?',const=packet.OUTPUT,type=Path)
    args=parser.parse_args()
    if args.check is not None:packet.verify(args.check)
    else:args.write.write_text(packet.canonical(packet.build()))
    print('PASS source_strong_coupling')


if __name__=='__main__':main()
