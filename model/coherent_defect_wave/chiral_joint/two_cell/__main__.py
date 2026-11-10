"""Write or independently replay the full two-cell certificate."""
import argparse
import json
from pathlib import Path
from .certificate import certificate
from ..receipt_replay import check_receipt


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    group.add_argument('--generate-witnesses',action='store_true',
                       help='explicitly replace dyadic eigenbasis proposals; then run --write')
    args=parser.parse_args()
    if args.generate_witnesses:
        from .spectrum import generate_witnesses
        print(f'PASS generated {generate_witnesses()} exact-checked eigenbasis witnesses')
        return
    path=Path(__file__).with_name('results.json')
    rendered=json.dumps(certificate(),indent=2,sort_keys=True)+'\n'
    if args.write:path.write_text(rendered)
    else:check_receipt(path,rendered,'complete two-cell')
    print('PASS full Gauss ground, complete charged records and finite propagation/refinement tests')


if __name__=='__main__':main()
