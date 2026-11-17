"""Write or replay independent charged-spectrum validation."""
import argparse,json
from pathlib import Path
from .certificate import certificate,validate_pins

def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    path=Path(__file__).with_name('results.json')
    if a.check:
        expected=json.loads(path.read_text());validate_pins(expected['input_sha256'])
    fresh=certificate()
    if a.write:path.write_text(json.dumps(fresh,indent=2,sort_keys=True)+'\n')
    elif fresh!=expected:raise ValueError('fresh independent charged-spectrum replay differs')
    print('PASS: independent Taylor, rational inertia and actual charged-ground comparison')
if __name__=='__main__':main()
