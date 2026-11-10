import argparse
import json
from pathlib import Path
from .certificate import certificate


def main():
    parser = argparse.ArgumentParser(description='Replay connected spatial defect-wave construction')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',nargs='?',const=str(Path(__file__).with_name('results.json')))
    args = parser.parse_args()
    content = json.dumps(certificate(),indent=2,sort_keys=True)+'\n'
    if args.write:
        Path(__file__).with_name('results.json').write_text(content)
    elif Path(args.check).read_text() != content:
        raise SystemExit('spatial defect-wave receipt differs')
    print('PASS connected local cell exchange, interior gauge constraints, causal cone and spectator band')


if __name__ == '__main__':
    main()
