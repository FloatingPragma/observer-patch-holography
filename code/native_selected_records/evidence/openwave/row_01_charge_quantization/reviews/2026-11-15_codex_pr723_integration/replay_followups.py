#!/usr/bin/env python3
"""Offline, read-only replay of the scoped PR723 follow-ups.

Requires Python 3 and already-installed NumPy. No download, install or archive
execution. --record creates a new receipt exclusively; default only prints it.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
R1 = 'row_01_charge_quantization/work/2026-11-15_codex_pr723_unit_selection'
R2 = 'row_02_electron_rest_energy/work/2026-11-15_codex_pr723_rest_selection'
R7 = 'row_07_double_cover/work/2026-11-15_codex_pr723_rotation'
R7A = 'row_07_double_cover/work/2026-11-15_codex_pr723_two_turn_rejections'
COMMON = 'row_01_charge_quantization/reviews/2026-11-15_codex_pr723_common_instrument'
BINDINGS = {
    R1: '35b67994df0c89d8478624815a95585293b99249da46b4dc90ba5fd27b72d5ae',
    R2: '8c6d86aecb8c516afeca52bbe5db8f6014271e4c00c1c50af81cfdb2db117366',
    R7: '05af1f898eb98a68413084b3b01586e600c79aa231d2f0318f0f90b887158f32',
    R7 + '/berry': '8add06ee1425a84e00997a7311f37ba222e7f77199487040506838cc21f00b7f',
    R7A: '1861b420260787509cc944f291919af901a87a7f5dc7ca8c86c91c7ab4cfe512',
    'row_01_charge_quantization/reviews/2026-11-15_codex_pr723_cross_review': '949b0ef26f557f96b582c7860be4e11f2d27bb07a96432eb7294f7de89384e5a',
    'row_02_electron_rest_energy/reviews/2026-11-15_codex_pr723_adversarial': 'b59a9c0e8d6f8b80a634563ec9e19d54fc5d2d16fb699c48cd8b5cedb417da16',
    'row_07_double_cover/reviews/2026-11-15_codex_pr723_adversarial': '760155ed834df197e359cd209580122f7b10f9c5ae8359af789c9a488f0a0caa',
    COMMON: 'd1d26e907932a04ad4ed3bb19cb556fd3e800161d56bcd9b0bd7987bbb85e10b',
}
COMMANDS = [(R1, 'headline_strict.py'), (R2, 'headline.py'),
            (R2, 'independent_check.py'), (R2, 'trajectory_check.py'),
            (R2, 'tunability.py'), (R2, 'tunability_check.py'),
            (R7, 'replay_all.py'), (R7A, 'headline.py'),
            (COMMON, 'test_common_instrument.py')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_manifests():
    result = {}
    for directory, digest in BINDINGS.items():
        base = ROOT / directory
        manifest = base / 'MANIFEST.sha256'
        if sha(manifest) != digest:
            raise ValueError(f'Manifest changed: {manifest}')
        leaves = []
        for line in manifest.read_text().splitlines():
            expected, name = line.split('  ', 1)
            path = (base / name).resolve()
            if not path.is_relative_to(base.resolve()):
                raise ValueError(f'Out-of-directory manifest entry: {name}')
            if sha(path) != expected:
                raise ValueError(f'Leaf changed: {path}')
            leaves.append(name)
        result[directory] = {'sha256': digest, 'files': len(leaves)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, help='Create a new receipt; never overwrite')
    args = parser.parse_args()
    if args.record and args.record.exists():
        parser.error('Receipt already exists; omit --record for a read-only replay')
    before = check_manifests()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    runs = []
    for directory, script in COMMANDS:
        start = time.monotonic()
        completed = subprocess.run([sys.executable, '-B', script], cwd=ROOT/directory,
                                   env=env, capture_output=True, text=True, timeout=60)
        if completed.returncode:
            raise RuntimeError(f'{directory}/{script}: {completed.stdout}\n{completed.stderr}')
        runs.append({'directory': directory, 'command': f'python3 -B {script}',
                     'exit_code': completed.returncode,
                     'seconds': round(time.monotonic()-start, 6),
                     'stdout_sha256': hashlib.sha256(completed.stdout.encode()).hexdigest(),
                     'stderr': completed.stderr})
    if before != check_manifests():
        raise ValueError('Custody changed during replay')
    result = {'passed': True, 'same_group': True, 'external_G6': False,
              'physical_rows_awarded': False, 'python': sys.version.split()[0],
              'manifests': before, 'runs': runs,
              'scope': 'Finite follow-up experiments, independent readers and complete outcome checks; no physical claim or external award inferred.'}
    content = json.dumps(result, indent=2) + '\n'
    if args.record:
        with args.record.open('x') as output:
            output.write(content)
    print(content, end='')


if __name__ == '__main__':
    main()
