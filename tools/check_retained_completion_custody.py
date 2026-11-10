#!/usr/bin/env python3
"""Run the unchanged 53 custody tests in their actual retained Git checkout.

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
CHECKPOINT = '0194fc685be3f7ba5f732681713c9cfd04e5290f'
TEST = 'model/completion/test_custody_renewal.py'
TEST_COUNT = 53
ORIGINAL_SHA256 = {
    TEST: '8ff269ae43229ba3f0f7e41663d090314400b42ff33eb6b461415acd59dd7e0d',
    'model/completion/r2063_custody_renewal.py': 'a19fc6ab5e88a9b75a5a2221f06e61d61acde2f6faeb9ad57688408337649fcc',
    'model/completion/custody_authority_r2063.json': '0f47735b000df58b670f9ae4b217f5dcfe8e404868adb3b502f9aa8b37d41f81',
    'pytest.ini': '38ee02e44262432ef0ec05ae311b3d05236687723fbbb6c77cba9edec32bbc9e',
}


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
    with tempfile.TemporaryDirectory(prefix='oph-retained-custody-0194-') as folder:
        view = Path(folder)
        extract_checkout(view, root)
        environment = dict(os.environ, PYTHONPATH=str(view), PYTHONDONTWRITEBYTECODE='1',
                           PYTEST_DISABLE_PLUGIN_AUTOLOAD='1')
        command = [sys.executable, '-m', 'pytest', '-q', '-Werror', TEST]
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
        return {'schema': 'oph.retained_completion_test_execution.v1', 'checkpoint': CHECKPOINT,
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
