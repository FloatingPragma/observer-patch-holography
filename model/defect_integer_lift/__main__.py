"""Regenerate or replay the integer-lift receipt."""
import argparse
import json
from pathlib import Path

from .certificate import certificate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = certificate()
    failed = [name for name, ok in value['checks'].items() if not ok]
    if failed:
        raise SystemExit('failed checks: ' + ', '.join(failed))
    path = Path(__file__).with_name('results.json')
    text = json.dumps(value, indent=1, sort_keys=True) + '\n'
    if args.write:
        path.write_text(text)
    elif json.loads(path.read_text()) != json.loads(text):
        raise SystemExit('the integer lift receipt differs from its complete recomputation')
    print(f"PASS: {len(value['checks'])} checks: closed carriers, the exact gates of the lift, kinematics, "
          'readings and imported run counts')


if __name__ == '__main__':
    main()
