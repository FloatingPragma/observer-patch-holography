"""Write or check the finite common chiral source receipt."""
import argparse,json
from pathlib import Path
from .certificate import certificate


def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    args=parser.parse_args();path=Path(__file__).with_name('results.json')
    rendered=json.dumps(certificate(),indent=2,sort_keys=True)+'\n'
    if args.write:path.write_text(rendered)
    elif path.read_text()!=rendered:raise SystemExit('FAIL stale common chiral source receipt')
    print('PASS common finite chiral gauge/Higgs Hamiltonian, physical filter and encounter')


if __name__=='__main__':main()
