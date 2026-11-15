#!/usr/bin/env python3
"""Replay frozen common-source instruments in independent isolated interpreters."""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json
import os
import shutil
import site
import subprocess
import sys
import sysconfig
import tempfile

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    public = json.loads((HERE / 'PUBLIC_MANIFEST.json').read_text())
    for relative, expected in public['files_sha256'].items():
        path = HERE / relative
        if path.is_symlink() or not path.is_file() or sha(path) != expected:
            raise ValueError('public package input mismatch: ' + relative)
    snapshot = json.loads((HERE / 'SNAPSHOT_MANIFEST.json').read_text())
    for relative, expected in snapshot['files_sha256'].items():
        path = HERE / 'workspace' / relative
        if path.is_symlink() or not path.is_file() or sha(path) != expected:
            raise ValueError('frozen source mismatch: ' + relative)
    actual = {str(p.relative_to(HERE / 'workspace')) for p in (HERE / 'workspace').rglob('*') if p.is_file()}
    if actual != set(snapshot['files_sha256']):
        raise ValueError('unmanifested files in workspace')
    return snapshot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Receipt destination')
    parser.add_argument('--mode', choices=('all', 'integration', 'quantum', 'reference', 'two_charge'), default='all')
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 13):
        raise RuntimeError('Use Python 3.13 with the dependencies in requirements.txt; see README.md.')
    snapshot = verify()
    versions = {}
    for line in (HERE / 'requirements.txt').read_text().splitlines():
        name, expected = line.split('==')
        versions[name] = importlib.metadata.version(name)
        if versions[name] != expected:
            raise RuntimeError('Numerical replay requires ' + line + '; found ' + versions[name])
    modes = ['integration', 'quantum', 'reference', 'two_charge'] if args.mode == 'all' else [args.mode]
    environment = dict(os.environ)
    environment.pop('PYTHONPATH', None)
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        environment[name] = '1'
    sites = sorted({str(p) for p in [sysconfig.get_path('purelib'), sysconfig.get_path('platlib'),
                                   site.getusersitepackages(), *site.getsitepackages()] if p and Path(p).is_dir()})
    records = {}
    with tempfile.TemporaryDirectory(prefix='oph-common-public-') as name:
        temporary = Path(name).resolve()
        forbidden = temporary / 'outside_workspace.py'
        forbidden.write_text('raise RuntimeError("mutable source must not execute")\n')
        for mode in modes:
            control = temporary / mode
            control.mkdir()
            shutil.copytree(HERE / 'workspace', control / 'workspace')
            for filename in ('worker.py', 'factors.py', 'SNAPSHOT_MANIFEST.json'):
                shutil.copyfile(HERE / filename, control / filename)
            aliases = {}
            if mode == 'two_charge':
                consumer = 'openwave-tuning/next/CHAIN2/two_charge_response/'
                for filename in ('INPUT_SNAPSHOT.json', 'IMPORT_ANNEX.json'):
                    data = json.loads((control / 'workspace' / consumer / filename).read_text())
                    for relative, expected in data['files'].items():
                        if snapshot['files_sha256'].get(relative) != expected:
                            raise ValueError('two-charge input disagrees with shared snapshot: ' + relative)
                        alias = consumer + 'inputs/' + relative
                        dest = control / 'workspace' / alias
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copyfile(control / 'workspace' / relative, dest)
                        aliases[alias] = relative
            (control / 'runtime.json').write_text(json.dumps({'site_directories': sites,
                'forbidden_probe': str(forbidden), 'runtime_aliases': aliases,
                'dependency_versions': versions}))
            generated = {'integration': ['openwave-tuning/common_dynamics/integration/INTEGRATION.json',
                                        'openwave-tuning/common_dynamics/integration/BLOCK_PROGRAM.json'],
                         'two_charge': ['openwave-tuning/next/CHAIN2/two_charge_independent/RESULT.json']}.get(mode, [])
            for relative in generated:
                (control / 'workspace' / relative).unlink()
            # Keep the selected interpreter spelling; resolving a venv symlink changes its environment.
            child = subprocess.run([sys.executable, '-I', '-S', str(control / 'worker.py'), mode],
                                   cwd=control, env=environment, text=True, capture_output=True)
            if child.returncode:
                raise RuntimeError(mode + ' replay failed:\n' + child.stdout + child.stderr)
            print(child.stdout.strip(), flush=True)
            record = json.loads((control / (mode + '.json')).read_text())
            for relative, expected in snapshot['files_sha256'].items():
                if sha(control / 'workspace' / relative) != expected:
                    raise ValueError('regenerated output or source differs: ' + relative)
            for alias, relative in aliases.items():
                if sha(control / 'workspace' / alias) != snapshot['files_sha256'][relative]:
                    raise ValueError('runtime input materialization changed')
            record['all_source_and_regenerated_output_bytes_match'] = True
            record['outputs_deleted_before_execution'] = generated
            record['runtime_aliases_materialized_from_flat_snapshot'] = len(aliases)
            records[mode] = record
    verify()
    result = {'schema': 'oph.common_dynamics.public_replay.v1', 'status': 'PASS_SCOPED_FINITE_INSTRUMENTS',
              'snapshot_manifest_sha256': sha(HERE / 'SNAPSHOT_MANIFEST.json'),
              'snapshot_file_count': snapshot['file_count'], 'executions': records,
              'official_rows_changed': [], 'official_row_awards': 0, 'all31_physical_complete': False,
              'one_interacting_continuum_construction_established': False}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'executions'}, indent=2))


if __name__ == '__main__':
    main()
