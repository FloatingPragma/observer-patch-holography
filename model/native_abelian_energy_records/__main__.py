"""Write or replay native-bin ground and ordered finite energy records."""
import argparse,json
from pathlib import Path
from .certificate import certificate,validate_pins

def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();path=Path(__file__).with_name('results.json')
    if args.check:
        expected=json.loads(path.read_text());validate_pins(expected['input_sha256'])
    fresh=certificate()
    if args.write:path.write_text(json.dumps(fresh,indent=2,sort_keys=True)+'\n')
    elif fresh!=expected:raise ValueError('fresh native energy record replay differs')
    print('PASS: actual native-bin ground, finite direct-H clock and complete ordered arrival/energy law')

if __name__=='__main__':main()
