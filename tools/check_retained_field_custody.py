#!/usr/bin/env python3
"""Run the unchanged 31 field-custody tests in their actual retained Git checkout.

The current worktree and pytest configuration are never copied into that view.
This verifies the dated custody contract; it neither validates today's source
pins against old receipts nor relabels old evidence as a current replay.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'c7b09261853f5a3c5a67e21294c3d9511ed59fc7'
TESTS = ('model/completion/test_field_custody.py', 'model/completion/test_field_navigation.py')
TEST_COUNT = 31
ORIGINAL_SHA256 = {'model/completion/test_field_custody.py': 'f8828194c39da3968a5d4ca970dd1cf8823346b658eda4ec0b9d8f776cf2c22a', 'model/completion/test_field_navigation.py': '3e0c38e1b2fb945257824502831b3094ad74c81b723b915397d69ca8a4d6b7ec', 'model/completion/field_custody.py': '3708ccbaa0bb2e10e15d0d43a388b7db0a7dbcbfffa7877e710b2eb13a14e930', 'model/completion/field_navigation.py': '5e3cd590e56ebbd2cb5b1e9b75114aad5eb0cfe61829bb7f42961ce52fd35d64', 'model/completion/field_custody_authority.json': '2417e7a60e23516d44603838d4f248b20fd221e47c4732d5ba53a401d8e5b49b', 'model/completion/field_navigation_authority.json': '0b8c5169a4396a62fd324c6109a4a677331b9a6ac4fec943321d77afbd47987a', 'pytest.ini': '38ee02e44262432ef0ec05ae311b3d05236687723fbbb6c77cba9edec32bbc9e'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def authenticate(view):
    for name, expected in ORIGINAL_SHA256.items():
        actual = hashlib.sha256((Path(view) / name).read_bytes()).hexdigest()
        require(actual == expected, 'retained checker, test, authority or configuration changed: ' + name)


def extract_checkout(destination, root=ROOT):
    actual = subprocess.check_output(['git', 'rev-parse', CHECKPOINT + '^{commit}'], cwd=root, text=True).strip()
    require(actual == CHECKPOINT, 'retained checkpoint identity changed')
    with subprocess.Popen(['git', 'archive', '--format=tar', CHECKPOINT], cwd=root,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
        with tarfile.open(fileobj=process.stdout, mode='r|') as archive:
            archive.extractall(destination, filter='data')
        errors = process.stderr.read().decode()
        require(process.wait() == 0, 'retained checkout extraction failed: ' + errors)
    authenticate(destination)


def execute(*, root=ROOT):
    with tempfile.TemporaryDirectory(prefix='oph-retained-field-custody-c7-') as folder:
        view = Path(folder)
        extract_checkout(view, root)
        environment = dict(os.environ, PYTHONPATH=str(view), PYTHONDONTWRITEBYTECODE='1',
                           PYTEST_DISABLE_PLUGIN_AUTOLOAD='1')
        command = [sys.executable, '-m', 'pytest', '-q', '-Werror', *TESTS]
        collected = subprocess.run(command[:4] + ['--collect-only'] + command[4:], cwd=view,
                                   env=environment, capture_output=True, text=True)
        require(collected.returncode == 0, 'retained collection failed:\n' + collected.stdout + collected.stderr)
        require(re.search(rf'\b{TEST_COUNT} tests collected\b', collected.stdout) is not None,
                'retained test inventory changed:\n' + collected.stdout)
        completed = subprocess.run(command, cwd=view, env=environment, capture_output=True, text=True)
        print(completed.stdout, end='')
        if completed.stderr:
            print(completed.stderr, file=sys.stderr, end='')
        require(completed.returncode == 0, 'retained custody regression failed')
        require(re.search(rf'\b{TEST_COUNT} passed\b', completed.stdout) is not None,
                'retained tests did not all execute')
        authenticate(view)
        return {'schema': 'oph.retained_field_custody_test_execution.v1', 'checkpoint': CHECKPOINT,
                'original_source_sha256': ORIGINAL_SHA256, 'tests_passed': TEST_COUNT,
                'original_pytest_configuration_used': True,
                'current_worktree_sources_injected': False,
                'current_live_source_custody_validated_by_this_historical_run': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = execute()
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
