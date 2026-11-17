"""Replay the unchanged-action charged dynamics, proof and complete source history."""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import math
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import compare_receipt
from . import certification, observables, source

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
CAMPAIGN = 'charged-action-dynamics-20261031-v1'
CONTRACT_SHA256 = 'af2a56e8af97c68dcfe986203d80950ca382b0db1417de64475dcbc67f286b24'


def contract():
    path = HERE/'contract.json'
    if hashlib.sha256(path.read_bytes()).hexdigest() != CONTRACT_SHA256:
        raise ValueError('charged dynamics contract changed')
    return strict_json.load(path)


def source_bindings():
    """All owning files and recursive repository imports, including eager parents."""
    pending = list(HERE.glob('*.py'))+[HERE/'README.md', HERE/'contract.json',
                                     ROOT/'model/requirements.txt', ROOT/'requirements.txt',
                                     ROOT/source.compiler.RESPONSE, ROOT/source.compiler.CONTROL_RECEIPT]
    seen = set()

    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        stem = name.replace('.', '/')
        candidates = [ROOT/(stem+'.py'), ROOT/stem/'__init__.py']
        parts = stem.split('/')
        candidates.extend(ROOT.joinpath(*parts[:end], '__init__.py')
                          for end in range(1, len(parts)))
        pending.extend(path for path in candidates if path.is_file() and path not in seen)

    while pending:
        path = pending.pop()
        if path in seen:
            continue
        if not path.is_file() or not path.resolve().is_relative_to(ROOT):
            raise ValueError('missing or escaping source dependency')
        seen.add(path)
        if path.suffix != '.py':
            continue
        module = path.relative_to(ROOT).with_suffix('').as_posix().replace('/', '.')
        package = (module.removesuffix('.__init__') if module.endswith('.__init__')
                   else module.rpartition('.')[0])
        add_module(module)
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'))):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add_module(target)
                for alias in node.names:
                    add_module(target+'.'+alias.name)
    return {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(seen)}


# A warm interpreter executes its imported code, even if the files have since
# changed. Bind that import lifetime as well as the start/end of each replay.
_IMPORTED_SOURCE_BINDINGS = tuple(source_bindings().items())


def report():
    declared = contract()
    pins = source_bindings()
    if pins != dict(_IMPORTED_SOURCE_BINDINGS):
        raise ArithmeticError('source files changed since packet import; restart the replay process')
    result = dict(schema='oph.charged_action_dynamics.receipt.v1', campaign_id=CAMPAIGN,
                  scope=declared['scope'], contract_sha256=CONTRACT_SHA256,
                  band_and_transfer=certification.certificate(),
                  field_observables=observables.certificate(), source_history=source.certificate(),
                  integrated_gate_accepted=False, physical_comparison_count_added=0,
                  physical_electron_identification=False, nontrivial_topological_defect_assigned=False,
                  source_bindings=pins)
    if source_bindings() != pins:
        raise ArithmeticError('source files changed during charged dynamics replay')
    compare(result, result)
    return result


def compare(actual, expected):
    """Exact schema/counts/proofs and finite floating diagnostics, never negative mass."""
    def probabilities(value, path=()):
        if isinstance(value, dict):
            for key, child in value.items():
                probabilities(child, (*path, key))
        elif isinstance(value, list):
            for child in value:
                probabilities(child, path)
        elif type(value) is float:
            if not math.isfinite(value):
                raise ValueError('nonfinite charged dynamics receipt')
            if any('probabilit' in key or key == 'block_number_bins'
                   or key.endswith(('weights', '_weight')) for key in path):
                if value < 0 or value > 1+2e-10:
                    raise ValueError('invalid charged dynamics probability')
    probabilities(actual)
    probabilities(expected)
    compare_receipt(actual, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = report()
    if args.check:
        compare(strict_json.load(args.check), result)
    else:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n',
                              encoding='utf-8', newline='\n')
    print('Verified complete charged-band continuation, physical currents and source history; no physical species assigned.')


if __name__ == '__main__':
    main()
