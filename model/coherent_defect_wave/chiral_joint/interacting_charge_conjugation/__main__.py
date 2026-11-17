"""Write or replay the exact antiunitary source symmetry certificate."""
import argparse,json
from pathlib import Path
from .certificate import certificate,validate_pins
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[4]))
    from tools.pin_scope import bound_pins

def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    path=Path(__file__).with_name('results.json')
    if a.check:
        expected=json.loads(path.read_text());validate_pins(bound_pins(expected['input_sha256']))
    fresh=certificate()
    if a.write:path.write_text(json.dumps(fresh,indent=2,sort_keys=True)+'\n')
    elif fresh!=expected:raise ValueError('fresh antiunitary charge/reflection replay differs')
    print('PASS: exact native antiunitary symmetry, paired source spectra and explicit reference boundary')
if __name__=='__main__':main()
