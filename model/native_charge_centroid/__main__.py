"""Write or replay the actual-preparation finite native centroid."""
import argparse,json
from pathlib import Path
from .certificate import certificate,validate

def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    path=Path(__file__).with_name('results.json')
    if args.check:
        expected=json.loads(path.read_text());validate(expected['input_sha256'])
    fresh=certificate()
    if args.write:path.write_text(json.dumps(fresh,indent=2,sort_keys=True)+'\n')
    elif fresh!=expected:raise ValueError('fresh native centroid replay differs')
    print('PASS: actual no-P0 preparation, full-native centroid and complete finite same-H phase records')

if __name__=='__main__':main()
