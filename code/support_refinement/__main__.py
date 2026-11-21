"""Produce or exactly replay the finite rational-refinement certificate."""
import argparse
import json
from pathlib import Path

from .construction import certificate, need


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    payload = json.dumps(certificate(), indent=2, sort_keys=True)+'\n'
    if args.write:
        args.write.write_text(payload)
    else:
        need(args.check.read_text() == payload, 'support-refinement certificate differs')
        print('Exact equivariant rational support retracts and homotopies verified.')


if __name__ == '__main__':
    main()
