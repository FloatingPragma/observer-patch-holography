#!/usr/bin/env python3
"""Replay the exact reviewed finite chain in a disposable, isolated workspace."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    manifest = json.loads((HERE / 'PUBLIC_MANIFEST.json').read_text())
    for relative, expected in manifest['files_sha256'].items():
        path = HERE / relative
        if path.is_symlink() or not path.is_file() or digest(path) != expected:
            raise ValueError('public input mismatch: ' + relative)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Write replay receipt to this path')
    args = parser.parse_args()
    manifest = verify()
    expected = json.loads((HERE / 'provenance/RESULT.json').read_text())
    with tempfile.TemporaryDirectory(prefix='oph-chain-public-') as temporary:
        bundle = Path(temporary).resolve() / 'bundle'
        shutil.copytree(HERE / 'provenance', bundle)
        shutil.copytree(HERE / 'workspace', bundle / 'frozen_workspace')
        environment = dict(os.environ)
        environment.pop('PYTHONPATH', None)
        environment['PYTHONDONTWRITEBYTECODE'] = '1'
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
            environment[name] = '1'
        subprocess.run([sys.executable, '-I', str(bundle / 'run.py')],
                       cwd=bundle, env=environment, check=True)
        actual = json.loads((bundle / 'RESULT.json').read_text())
        runtime = actual['snapshot_execution'].pop('runtime_isolation')
        expected_runtime = expected['snapshot_execution'].pop('runtime_isolation')
        if actual != expected:
            raise ValueError('replayed scientific result differs from reviewed result')
        for name in ('historical_workspace', 'all_successful_project_reads_inside_snapshot',
                     'mutable_source_read_probe_rejected', 'site_startup_disabled',
                     'environment_pythonpath_ignored'):
            if runtime[name] is not True:
                raise ValueError('isolation condition failed: ' + name)
        if runtime['snapshot_module_files'] != expected_runtime['snapshot_module_files']:
            raise ValueError('local imported-module closure differs from reviewed execution')
        if runtime['blocked_mutable_reads'] != expected_runtime['blocked_mutable_reads']:
            raise ValueError('unexpected attempt to read mutable project source')
        execution = json.loads((bundle / 'SNAPSHOT_EXECUTION.json').read_text())
        if execution['all_snapshot_bytes_unchanged'] is not True:
            raise ValueError('snapshot custody failed')
        receipt = {
            'schema': 'oph.chain.public_replay.v1',
            'status': 'PASS_SCOPED_FINITE_COMPOSITIONS',
            'public_input_count': len(manifest['files_sha256']),
            'original_result_sha256': digest(HERE / 'provenance/RESULT.json'),
            'original_snapshot_manifest_sha256': digest(HERE / 'provenance/SNAPSHOT_MANIFEST.json'),
            'scientific_result_exactly_reproduced': True,
            'runtime_metadata_compared_separately': runtime,
            'snapshot_file_count': execution['snapshot_file_count'],
            'all_snapshot_bytes_unchanged': True,
            'rejected_mutations': actual['mutation_count'],
            'official_rows_changed': [],
            'all31_physical_complete': actual['all31_physical_complete'],
            'common_physics_complete': actual['common_physics_complete'],
        }
    verify()
    if args.output:
        args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items()
                      if key != 'runtime_metadata_compared_separately'}, indent=2))


if __name__ == '__main__':
    main()
