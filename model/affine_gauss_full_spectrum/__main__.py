"""Verify a full-source residual; optionally write the fixed admitted receipt."""
import argparse,json
from pathlib import Path
from .certificate import certificate,replay
from model.affine_gauss_spectrum.model import validate_pins
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--check',action='store_true');group.add_argument('--write',action='store_true')
    args=parser.parse_args();path=Path(__file__).with_name('results.json')
    if args.check:
        saved=json.loads(path.read_text());validate_pins(bound_pins(saved['input_sha256']))
    fresh=certificate()
    if args.write:path.write_text(json.dumps(fresh,indent=2)+'\n')
    elif fresh!=saved:raise ValueError('fixed charged spectral certificate differs')
    print('PASS: complete70928-state charged-ground residual and readout bounds')
    print(json.dumps(replay()[3],sort_keys=True))


if __name__=='__main__':main()
