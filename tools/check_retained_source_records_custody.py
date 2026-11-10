#!/usr/bin/env python3
"""Run the unchanged 73 source-record custody and routing tests in their actual retained Git checkout.

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
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = '7fc48ae176316cb6d2b7ef6decb2060dab9f7db8'
TESTS = ('model/completion/test_source_records_navigation.py', 'model/completion/test_source_records_custody.py', 'tools/test_source_records_atlas.py', 'tools/test_source_records_comparison.py', 'tools/test_source_record_readouts.py', 'tools/test_source_record_manifests.py', 'tools/test_retained_custody_routing.py', 'tools/test_retained_field_custody.py')
TEST_COUNT = 73
ORIGINAL_SHA256 = {'model/completion/source_records_custody.py': 'e303faa262771249421a5775487ccb03d47ad0305980f0a64a819e5124a74d8a', 'model/completion/source_records_custody_authority.json': '8594f608e26e2086432ed417de1737f7aaf9e6f7bd5830e8c1c5656de0f9ba0e', 'model/completion/source_records_navigation.py': '423a6676c5c6462f2b6c22ff89f5a18dac68ffec4d5ff72d3686596148dd1288', 'model/completion/source_records_navigation_authority.json': 'e52da9fe4c7dd9e612351af36f1bf61e0168a46f1d054fc25e3f0e9fac2a8be0', 'model/completion/test_source_records_custody.py': '78426b8abadae6e02f48485ed52fc325339b74e31d2e17116544788115702be6', 'model/completion/test_source_records_navigation.py': 'f4f45318b831d86a79b7b18cb6293458ac3aa8aba2b57d3083b6cba0b6176cb1', 'pytest.ini': 'e2e34163f479f213fca5a6d39fdde600e884e020f56b1387f2d288c441e1a695', 'tools/renew_source_record_manifests.py': '52e031e50b5d5f395d647c40850c8bd9c27ab3be01cae054e767d3faf5956e38', 'tools/renew_source_record_readouts.py': '428e0e33a5ee3e2a785ae7050e7819a1c77c35894b6cc0667e9d251fe535cbf8', 'tools/renew_source_records_atlas.py': '2ff18d131cd59d4276c44e2833bef36b8bcf1e8270b4eeb2d16dc0f09bcde3d0', 'tools/renew_source_records_comparison.py': 'd6da3661305bf60410e6f4a3b3d72efbf660735f8abe2c3eb3e9a7fa34315f66', 'tools/source_record_manifest_authority.json': 'a7cb0c314049fc18212f3785e4b348b5816fb06b3b4562bcd28e69be539a8992', 'tools/source_record_readout_owners.json': '5ce926726f445fba143bd4b18936596fa72f9e12481185031e58184bb73b87dc', 'tools/source_records_atlas_authority.json': '95bb023c369a5e108992b624be65803b95dada499f06bbfccfcb55e7c84b9e82', 'tools/test_source_record_manifests.py': '38fd6c5e90da3370ee70d4299953b9adc73f313d50df2391f687c2852a2d17a6', 'tools/test_source_record_readouts.py': '8b4e94e96b75053a50c8b5040d3ff4f0c5f187f55f50af1e862d42d06fc57103', 'tools/test_source_records_atlas.py': 'c88d6fc9c3b9a237183de1c9a6fd79a18fd2f19acde7ce220d73dc8f3c682ceb', 'tools/test_source_records_comparison.py': 'd5e1d3fd55a414b946a6e3faf3c88b2383bea4b7883b00fbc493e04530acb302', 'tools/test_retained_custody_routing.py': 'bb33e38470adc63c78aac466090d1529dc34afb62470b4a884b55c9689a6b59b', 'tools/test_retained_field_custody.py': '5e619f2e61be5015dcb090e2a89ecd633c19590632a38600bd625b2ac8b8d162'}


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
    # The unchanged historical routing controls invoke git show on their own
    # authenticated ancestor commits. A private shared clone supplies those
    # immutable objects without injecting any current worktree source.
    subprocess.run(['git', 'clone', '--shared', '--no-checkout', str(Path(root).resolve()), str(destination)],
                   check=True, capture_output=True, text=True)
    subprocess.run(['git', 'checkout', '--detach', CHECKPOINT], cwd=destination,
                   check=True, capture_output=True, text=True)
    authenticate(destination)


def execute(*, root=ROOT):
    with tempfile.TemporaryDirectory(prefix='oph-retained-source-record-custody-7fc-') as folder:
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
        return {'schema': 'oph.retained_source_records_custody_test_execution.v1', 'checkpoint': CHECKPOINT,
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
