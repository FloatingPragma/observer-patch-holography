import argparse
import json
from pathlib import Path
from .certificate import certificate
from ..receipt_replay import check_receipt


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();rendered=json.dumps(certificate(),indent=2,sort_keys=True)+'\n'
    path=Path(__file__).with_name('results.json')
    if args.check:check_receipt(path,rendered,'same-H charge localization')
    else:path.write_text(rendered)
    print('PASS: unscreened string cost, full-bin localization and dynamical screening counterexample')


if __name__=='__main__':main()
