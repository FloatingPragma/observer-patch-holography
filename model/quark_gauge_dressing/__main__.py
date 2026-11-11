import argparse
import json
from .dressing import HERE, certificate


def main():
    parser = argparse.ArgumentParser()
    choices = parser.add_mutually_exclusive_group(required=True)
    choices.add_argument('--write', action='store_true')
    choices.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    elif result != json.loads(path.read_text()):
        raise ArithmeticError('the exact positive gauge dressing receipt differs')
    print('PASS: full-source quark symmetry and positive Gauss electric dressing')


if __name__ == '__main__':
    main()
