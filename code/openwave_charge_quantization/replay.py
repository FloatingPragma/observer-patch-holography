#!/usr/bin/env python3
"""Replay the complete native charge criterion and independent reviews."""
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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    public = HERE / 'PUBLIC_MANIFEST.json'
    if not public.is_file():
        raise ValueError('Public manifest required')
    for name, expected in json.loads(public.read_text())['files_sha256'].items():
        p = HERE / name
        if p.is_symlink() or sha(p) != expected:
            raise ValueError('Public package mismatch: ' + name)
    snapshot = json.loads((HERE / 'SNAPSHOT_MANIFEST.json').read_text())
    actual = {str(p.relative_to(HERE / 'workspace'))
              for p in (HERE / 'workspace').rglob('*') if p.is_file()}
    if actual != set(snapshot['files_sha256']):
        raise ValueError('Snapshot file set differs')
    for name, expected in snapshot['files_sha256'].items():
        p = HERE / 'workspace' / name
        if p.is_symlink() or sha(p) != expected:
            raise ValueError('Source mismatch: ' + name)
    return snapshot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('all', 'candidate', 'auditor_root', 'auditor_dyn17', 'codex'), default='all')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    snapshot = verify()
    modes = ['candidate', 'auditor_root', 'auditor_dyn17', 'codex'] if args.mode == 'all' else [args.mode]
    results = {}
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    with tempfile.TemporaryDirectory(prefix='oph-native-charge-') as temporary:
        root = Path(temporary).resolve()
        for mode in modes:
            home = root / mode
            shutil.copytree(HERE / 'workspace', home / 'workspace')
            shutil.copyfile(HERE / 'worker.py', home / 'worker.py')
            shutil.copyfile(HERE / 'SNAPSHOT_MANIFEST.json', home / 'SNAPSHOT_MANIFEST.json')
            child = subprocess.run([sys.executable, '-I', '-S', str(home / 'worker.py'), mode],
                                   cwd=home, env=env, capture_output=True, text=True)
            if child.returncode:
                raise RuntimeError(mode + ':\n' + child.stdout + child.stderr)
            results[mode] = json.loads((home / 'REPLAY.json').read_text())
            print(mode + ': PASS', flush=True)
    verify()
    result = {'schema': 'oph.native_charge.public_replay.v1',
              'status': 'PASS_DECLARED_FULL_STATIC_CRITERION',
              'candidate_id': 'literal-native-relative-spatial-charge-AU5-record-only-v1',
              'snapshot_manifest_sha256': sha(HERE / 'SNAPSHOT_MANIFEST.json'),
              'executions': results,
              'scope': 'Signed relative topological charge; total native load one or two; true G; faithful closing records; complete retained failures.',
              'official_rows_changed': [], 'all31_complete': False}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'executions'}, indent=2))


if __name__ == '__main__':
    main()
