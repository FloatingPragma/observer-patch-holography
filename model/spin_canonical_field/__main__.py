import argparse
import json
from . import packet


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = packet.HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(packet.build(), indent=2, sort_keys=True)+'\n')
    else:
        packet.validate(json.loads(path.read_text()))
    print('spin_canonical_field: PASS')


if __name__ == '__main__':
    main()
