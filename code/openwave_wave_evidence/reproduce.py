#!/usr/bin/env python3
"""Replay the bundled wave evidence without Git or network access.

Scientific computations and comparisons are the sealed amendment's code.
The tests and saved-result verifier receive a local copy of the sealed source
through an explicit replacement of their Git extraction function.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / 'evidence'
BUNDLE = EVIDENCE / 'shared/unified_wave_submission_20261115'
AMENDMENT = EVIDENCE / 'shared/unified_wave_replay_amendment_20261115'
REVIEW = EVIDENCE / 'notes/wave_review_20261115'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sealed_directory(root, anchor):
    manifest = root / 'MANIFEST.sha256'
    require(digest(manifest) == anchor['manifest_sha256'], 'Manifest digest mismatch')
    expected = {'MANIFEST.sha256'}
    for line in manifest.read_text(encoding='utf-8').splitlines():
        wanted, name = line.split(maxsplit=1)
        name = name.lstrip('*')
        require(name not in expected, 'Duplicate manifest entry')
        require(name and '\\' not in name and ':' not in name, 'Invalid manifest path')
        path = root / name
        require(path.resolve().is_relative_to(root.resolve()) and path.resolve() != root.resolve(),
                'Path escapes evidence')
        expected.add(name)
        require(digest(path) == wanted, 'Evidence digest mismatch: ' + name)
    actual = set()
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Symlink in evidence')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
    require(actual == expected, 'Evidence inventory mismatch')
    require(len(expected) - 1 == anchor['entries'], 'Evidence count mismatch')
    return anchor['entries']


def verify():
    require(sys.flags.optimize == 0, 'Optimized Python is unsupported')
    seals = json.loads((HERE / 'SEALS.json').read_bytes())
    counts = {}
    for relative, anchor in seals.items():
        require(relative in ('shared/unified_wave_submission_20261115',
                             'shared/unified_wave_replay_amendment_20261115',
                             'notes/wave_review_20261115'), 'Unknown seal')
        counts[relative] = sealed_directory(EVIDENCE / relative, anchor)
    require(len(counts) == 3, 'Missing evidence seal')
    return counts


def local_source(destination, policy):
    """Supply byte-identical input to frozen tests/verifier; no science changes."""
    require(policy['bundle'] == 'shared/unified_wave_submission_20261115', 'Wrong source')
    target = destination / policy['bundle']
    shutil.copytree(BUNDLE, target)
    return target


def load_amendment():
    sys.path.insert(0, str(AMENDMENT))
    return importlib.import_module('replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('custody')
    for name in ('replay', 'independent', 'tests', 'verify'):
        command = commands.add_parser(name)
        command.add_argument('--output', type=Path, required=True)
        if name in ('tests', 'verify'):
            command.add_argument('--run', type=Path, required=True,
                                 nargs='+' if name == 'verify' else None)
    args = parser.parse_args()
    counts = verify()
    if args.command == 'custody':
        print(json.dumps({'status': 'pass', 'entries': counts}, indent=2))
        return
    output = args.output.resolve()
    require(not output.exists(), 'Choose a fresh output path')
    require(not output.is_relative_to(HERE), 'Outputs must be outside this package')
    if args.command in ('replay', 'independent'):
        script = AMENDMENT / 'replay.py' if args.command == 'replay' else REVIEW / 'independent.py'
        completed = subprocess.run([sys.executable, '-B', str(script), '--bundle', str(BUNDLE),
                                    '--output', str(output)], check=False)
        verify()
        raise SystemExit(completed.returncode)
    r = load_amendment()
    with patch.object(r, 'extract_git', local_source):
        if args.command == 'tests':
            test_module = importlib.import_module('test_replay')
            test_module.RUN = args.run.resolve()
            suite = unittest.defaultTestLoader.loadTestsFromModule(test_module)
            result = unittest.TextTestRunner(verbosity=2).run(suite)
            passed = result.wasSuccessful() and result.testsRun == 21 and not result.skipped
            r.save(output, {'status': 'pass' if passed else 'fail', 'tests_run': result.testsRun,
                            'skipped': len(result.skipped), 'failures': len(result.failures),
                            'errors': len(result.errors), 'source_fixture': 'bundled exact bytes',
                            'sealed_test_sha256': digest(AMENDMENT / 'test_replay.py')})
            require(passed, 'Frozen regression suite failed')
        else:
            verifier = importlib.import_module('verify_results')
            with patch.object(sys, 'argv', ['verify_results.py',
                              *[str(run.resolve()) for run in args.run], '--output', str(output)]):
                verifier.main()
    verify()


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'Failed: {type(error).__name__}: {error}', file=sys.stderr)
        raise SystemExit(1)
