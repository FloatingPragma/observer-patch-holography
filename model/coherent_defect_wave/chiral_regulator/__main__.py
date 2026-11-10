import argparse
import json
from pathlib import Path
from .certificate import certificate


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args();value=certificate();path=Path(__file__).with_name('results.json')
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    elif json.loads(path.read_text())!=value:raise SystemExit('FAIL: chiral spectral regulator receipt differs')
    print('PASS: complete left-Weyl spectrum, filled reference and finite-nonlocality boundary')


if __name__=='__main__':main()
