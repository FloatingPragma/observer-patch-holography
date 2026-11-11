import argparse
import json
from .closure import HERE, certificate, write


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        write()
    elif certificate() != json.loads((HERE/'results.json').read_text()):
        raise ArithmeticError('finite observable theorem packet differs')
    print('PASS: complete periodic Y, joint clocks and all-lag field recurrence')


if __name__ == '__main__':
    main()
