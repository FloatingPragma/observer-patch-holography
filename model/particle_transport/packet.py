"""Reexecute the same-candidate transport and number-conserving detectors."""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import math
from pathlib import Path

from tools import strict_json
from ..observer.source_scalar import compiler
from ..common_source import slow_geometry
from . import algebra, detector, source


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'results.json'
CAMPAIGN = 'charged-particle-transport-20261031-v1'
# The frozen contract is checked separately from mutable source custody.
CONTRACT_SHA256 = '65337a2188bfcf7f6268108a9c06cda8593880770831e658f1f6d71458142945'
REPLAY_ATOL = 2e-10


def source_bindings():
    """Bind all local imports, eager initializers and explicit data inputs.

    Test code is included: an independent check is part of the retained
    evidence. External libraries are described by the pinned requirements.
    This static walk does not depend on which imports happen to be cached.
    """
    pending = [*HERE.glob('*.py'), HERE/'README.md', HERE/'contract.json',
               ROOT/'requirements.txt', ROOT/'model/requirements.txt',
               ROOT/compiler.RESPONSE, ROOT/compiler.CONTROL_RECEIPT,
               slow_geometry.SPEC, slow_geometry.ARCHIVE/'oph_exact/carrier.py',
               slow_geometry.ARCHIVE/'oph_fpe/dynamics/self_readback_repair_closure.py']
    visited = set()
    while pending:
        path = pending.pop().resolve()
        if path in visited:
            continue
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError(f'missing or escaping source artifact: {path}')
        visited.add(path)
        if path.suffix != '.py':
            continue
        import_root = slow_geometry.ARCHIVE if path.is_relative_to(slow_geometry.ARCHIVE) else ROOT
        parts = path.relative_to(import_root).with_suffix('').parts
        package = '.'.join(parts[:-1])
        for end in range(1, len(parts)):
            initializer = import_root.joinpath(*parts[:end], '__init__.py')
            if initializer.is_file():
                pending.append(initializer)
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'))):
            names = []
            if isinstance(node, ast.Import):
                names.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                name = '.'*node.level+(node.module or '')
                name = importlib.util.resolve_name(name, package) if node.level else name
                names.append(name)
                names.extend(name+'.'+alias.name for alias in node.names if alias.name != '*')
            for name in names:
                top = name.split('.')[0]
                if top in ('model', 'tools', 'code'):
                    module_root = ROOT
                elif top in ('oph_exact', 'oph_fpe'):
                    module_root = slow_geometry.ARCHIVE
                else:
                    continue
                stem = module_root.joinpath(*name.split('.'))
                pending.extend(candidate for candidate in (stem.with_suffix('.py'), stem/'__init__.py')
                               if candidate.is_file())
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(visited)}


_IMPORTED_BINDINGS = source_bindings()


def compare(expected, actual, path='receipt'):
    """Every key, type, count and scope is mandatory; floats have one fixed atol."""
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
                or abs(expected-actual) > REPLAY_ATOL):
            raise ValueError(f'{path}: numerical replay changed')
        field = path.rsplit('.', 1)[-1]
        if (field in ('preparation_abort_probability', 'continuation_abort_probability')
                or field.startswith('pointer_probabilities[')):
            if not 0 <= expected <= 1 or not 0 <= actual <= 1:
                raise ValueError(f'{path}: outcome probability outside [0, 1]')
    elif type(actual) not in (str, int, bool, type(None)) or expected != actual:
        raise ValueError(f'{path}: receipt value changed')


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('sources changed since packet import; replay in a fresh interpreter')
    if pins['model/particle_transport/contract.json'] != CONTRACT_SHA256:
        raise ValueError('charged-particle transport contract changed')
    contract = strict_json.load(HERE/'contract.json')
    compiler._source.cache_clear()
    compiler._source()
    slow_geometry.source.cache_clear()
    slow_geometry.source()
    result = {
        'schema': 'oph.charged_particle_transport.receipt.v1',
        'campaign_id': CAMPAIGN,
        'scope': contract['scope'],
        'contract_sha256': CONTRACT_SHA256,
        'same_candidate_algebra': algebra.certificate(),
        'source_transport': source.report(),
        'number_conserving_detector': detector.report(),
        'integrated_gate_accepted': False,
        'physical_particle_identification': False,
        'physical_comparison_count_added': 0,
        'failed_dispersion_benchmark_repaired': False,
        'native_topological_protection_established': False,
        'replay_float_absolute_tolerance': REPLAY_ATOL,
        'source_bindings': pins,
    }
    if source_bindings() != pins:
        raise ArithmeticError('sources changed during charged-particle replay')
    # Reject NaN and unsupported payloads even on the write path.
    json.dumps(result, allow_nan=False)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = report()
    if args.check:
        compare(strict_json.load(args.check), result)
    else:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print('Verified same-candidate transport, complete exchange and number-conserving readouts.')


if __name__ == '__main__':
    main()
