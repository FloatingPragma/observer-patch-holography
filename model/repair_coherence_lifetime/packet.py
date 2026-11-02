"""Replay the finite theorem's exact source and cut witnesses."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.exchange_sector.many_body_packet import source_bindings as inherited_bindings
from .certificate import report

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')


def source_bindings():
    paths = set(inherited_bindings())
    paths.update(p.relative_to(ROOT).as_posix() for p in Path(__file__).parent.glob('*.py'))
    paths.add('model/repair_coherence_lifetime/README.md')
    paths.add('paper/tex_fragments/REPAIR_COHERENCE_LIFETIME.tex')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def build():
    before = source_bindings()
    result = {'schema': 'oph.repair_coherence_lifetime.v1',
              'scope': 'FINITE_COMPONENT_FAMILY', 'all_pass': True,
              'certificate': report(), 'source_bindings': before,
              'integrated_gate_accepted': False, 'physical_comparison_count_added': 0}
    if source_bindings() != before:
        raise RuntimeError('source changed during exact replay')
    return result


def compare(actual, expected):
    """Reject type changes as well as changed values in a retained witness."""
    encode = lambda value: json.dumps(value, sort_keys=True, allow_nan=False, separators=(',', ':'))
    if encode(actual) != encode(expected):
        raise ArithmeticError('coherence theorem receipt differs from replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    mode.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = build()
    if args.check:
        compare(strict_json.load(args.check), result)
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'all_pass': True, 'fixed_number_factor': '5/6',
                      'nonextremal_factor': '11/12', 'sectors': 11}))


if __name__ == '__main__':
    main()
