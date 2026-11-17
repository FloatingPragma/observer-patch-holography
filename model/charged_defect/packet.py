"""Replay the charged-defect branch, spectrum and full-extension certificate."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from tools import strict_json
from ..particle_transport import packet as predecessor
from . import branches, dictionary


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
CONTRACT_SHA256 = 'f7e6248c1060ebbe8ab48341878b2ae56ef27ec5a758335cabf152f21db62068'
REPLAY_ATOL = 1e-13


def source_bindings():
    # The predecessor walks the transitive local-import graph, its eagerly
    # imported parents, pinned native source and explicit input artifacts.
    # Every new local source and its independent checks are additionally bound.
    result = predecessor.source_bindings()
    paths = [*HERE.glob('*.py'), HERE/'README.md', HERE/'contract.json',
             ROOT/'paper/tex_fragments/CHARGED_DEFECT_PROTECTED_ALGEBRA.tex']
    for path in paths:
        result[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(result.items()))


_IMPORTED_BINDINGS = source_bindings()


def compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual):
        raise ValueError(f'{path}: receipt type changed')
    if type(actual) is dict:
        if expected.keys() != actual.keys():
            raise ValueError(f'{path}: receipt keys changed')
        for key in actual:
            compare(expected[key], actual[key], path+'.'+key)
    elif type(actual) is list:
        if len(expected) != len(actual):
            raise ValueError(f'{path}: receipt length changed')
        for index, (left, right) in enumerate(zip(expected, actual, strict=True)):
            compare(left, right, f'{path}[{index}]')
    elif type(actual) is float:
        if (not math.isfinite(expected) or not math.isfinite(actual)
                or expected < 0 or actual < 0 or abs(expected-actual) > REPLAY_ATOL):
            raise ValueError(f'{path}: numerical replay changed')
    elif type(actual) not in (str, int, bool, type(None)) or expected != actual:
        raise ValueError(f'{path}: receipt value changed')


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('charged-defect sources changed; replay in a fresh interpreter')
    if pins['model/charged_defect/contract.json'] != CONTRACT_SHA256:
        raise ValueError('charged-defect contract changed')
    contract = strict_json.load(HERE/'contract.json')
    value = {'schema': 'oph.charged_defect.receipt.v1',
             'campaign_id': contract['campaign_id'], 'scope': contract['scope'],
             'contract_sha256': CONTRACT_SHA256,
             'branch_algebra': branches.report(),
             'same_candidate_dictionary': dictionary.report(),
             'integrated_gate_accepted': False,
             'physical_comparison_count_added': 0,
             'physical_particle_identification': False,
             'bare_native_law_modified': False,
             'full_field_source_execution_claimed': False,
             'source_bindings': pins,
             'replay_float_absolute_tolerance': REPLAY_ATOL}
    if source_bindings() != pins:
        raise ValueError('charged-defect sources changed during replay')
    json.dumps(value, allow_nan=False)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    value = report()
    if args.check:
        compare(strict_json.load(args.check), value)
    else:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print('Verified charged-defect complete branch algebra, same gap and joined-field obstruction.')


if __name__ == '__main__':
    main()
