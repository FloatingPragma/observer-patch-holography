"""Replay exact field-cone matching and restricted radial countercontrols."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import cone
from ..native_lorentz_transition.certificate import SOURCE_FILES as LORENTZ_SOURCES

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple(sorted(set((*cone.SOURCE_FILES, *LORENTZ_SOURCES,
                                'model/field_cone_compatibility/certificate.py'))))


def bindings():
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in SOURCE_FILES}


def build():
    before = bindings()
    result = cone.certificate()
    if bindings() != before:
        raise ValueError('field-cone source changed during replay')
    return dict(result, source_bindings=before, all_pass=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    mode.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    packet = build()
    if args.write is not None:
        args.write.write_text(json.dumps(packet, indent=2, sort_keys=True, allow_nan=False)+'\n')
        print('field-cone compatibility certificate written')
    else:
        if json.loads(args.check.read_text()) != packet:
            raise SystemExit('field-cone compatibility certificate differs')
        print('field-cone compatibility certificate matches')


if __name__ == '__main__':
    main()
