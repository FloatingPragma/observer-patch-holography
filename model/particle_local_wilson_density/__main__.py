"""Replay the exact native identity and same-state local spectral bounds."""
import argparse
import hashlib
import json
from pathlib import Path

from .density import SOURCE_INPUTS, certificate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'results.json'


def build():
    paths = set(SOURCE_INPUTS)
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.suffix in ('.py', '.md'))
    def pins():
        return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                for p in sorted(paths)}
    before = pins()
    value = certificate()
    value['input_sha256'] = before
    if pins() != before:
        raise ValueError('source changed during the local spectral replay')
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)+'\n'


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
        raise ValueError('stale local Wilson spectral certificate')
    print('PASS: native Wilson identity, local spectral weight and complete effects')


if __name__ == '__main__':
    main()
