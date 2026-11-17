"""Replay the full native sign-response boundary certificate."""
import argparse
import hashlib
import json
from pathlib import Path

from .boundary import SOURCE_INPUTS, certificate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'results.json'


def build():
    paths = set(SOURCE_INPUTS)
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.suffix in ('.py', '.md'))
    def pins():
        return bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                for p in sorted(paths)})
    before = pins()
    answer = certificate()
    answer['input_sha256'] = bound_pins(before)
    if pins() != before:
        raise ValueError('source changed during the sign-response replay')
    return json.dumps(answer, indent=2, sort_keys=True, ensure_ascii=False)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = build()
    if args.write:
        OUTPUT.write_text(value)
    elif not OUTPUT.exists() or OUTPUT.read_text() != value:
        raise ValueError('stale sign-response boundary certificate')
    print('PASS: native sign sensitivity, full block census and original records')


if __name__ == '__main__':
    main()
