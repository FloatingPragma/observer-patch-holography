#!/usr/bin/env python3
"""Replay frozen common-source instruments in independent isolated interpreters."""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json
import math
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


def compare_diagnostics(expected, observed, allowed, path=''):
    """Allow roundoff only in named diagnostics; certificates stay exact."""
    if type(expected) is not type(observed):
        raise ValueError('diagnostic type changed: ' + path)
    if isinstance(expected, dict):
        if expected.keys() != observed.keys():
            raise ValueError('diagnostic keys changed: ' + path)
        return sum((compare_diagnostics(value, observed[key], allowed, path+'/'+key)
                    for key, value in expected.items()), [])
    if isinstance(expected, list):
        if len(expected) != len(observed):
            raise ValueError('diagnostic length changed: ' + path)
        return sum((compare_diagnostics(value, observed[i], allowed, path+'/'+str(i))
                    for i, value in enumerate(expected)), [])
    if isinstance(expected, float) and path.split('/')[-1] in allowed:
        if not math.isfinite(expected) or not math.isfinite(observed) or abs(expected-observed) > 5e-13:
            raise ValueError('diagnostic exceeds roundoff allowance: ' + path)
        return [] if expected == observed else [{'path': path, 'archived': expected,
                    'replayed': observed, 'absolute_difference': abs(expected-observed)}]
    if expected != observed:
        raise ValueError('exact scientific value changed: ' + path)
    return []


def normalize_source_paths(value, equivalences, snapshot):
    """Normalize only declared identical frozen copies in the read-path receipt."""
    changes = []
    if not equivalences:
        return value, changes
    value = dict(value)
    references = {}
    for name, digest in value['source_files_read'].items():
        target = equivalences.get(name, name)
        if snapshot.get(name) != digest or snapshot.get(target) != digest:
            raise ValueError('source-path equivalence changes bytes: ' + name)
        if target in references:
            raise ValueError('source-path equivalence merges distinct reads')
        references[target] = digest
        if target != name:
            changes.append({'replayed_path': name, 'archived_path': target,
                            'identical_source_sha256': digest})
    value['source_files_read'] = references
    return value, changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Receipt destination')
    parser.add_argument('--mode', choices=('all', 'topology', 'native', 'terminal', 'terminal_opposing', 'terminal_channel'), default='all')
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
    plans = json.loads((HERE / 'PLANS.json').read_text())
    modes = list(plans) if args.mode == 'all' else [args.mode]
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
            for filename in ('worker.py', 'PLANS.json', 'SNAPSHOT_MANIFEST.json'):
                shutil.copyfile(HERE / filename, control / filename)
            aliases = {}
            if 'alias_inputs' in plans[mode]:
                data = json.loads((control / 'workspace' / plans[mode]['alias_inputs']).read_text())
                for relative, expected in data[plans[mode].get('alias_group', 'files_sha256')].items():
                    if snapshot['files_sha256'].get(relative) != expected:
                        raise ValueError('review input differs from snapshot: ' + relative)
                    alias = plans[mode]['alias_directory'] + relative
                    dest = control / 'workspace' / alias
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(control / 'workspace' / relative, dest)
                    aliases[alias] = relative
            (control / 'runtime.json').write_text(json.dumps({'site_directories': sites,
                'forbidden_probe': str(forbidden), 'runtime_aliases': aliases,
                'dependency_versions': versions}))
            generated = plans[mode]['outputs']
            for relative in generated:
                (control / 'workspace' / relative).unlink()
            withheld = plans[mode].get('withheld_review_freezes', [])
            for relative in withheld:
                # A frozen reviewer refuses to overwrite its own freeze. Only
                # this temporary verified copy omits that admission marker.
                assert relative in snapshot['files_sha256'] and relative not in generated
                (control / 'workspace' / relative).unlink()
            # Keep the selected interpreter spelling; resolving a venv symlink changes its environment.
            child = subprocess.run([sys.executable, '-I', '-S', str(control / 'worker.py'), mode],
                                   cwd=control, env=environment, text=True, capture_output=True)
            if child.returncode:
                raise RuntimeError(mode + ' replay failed:\n' + child.stdout + child.stderr)
            print(child.stdout.strip(), flush=True)
            record = json.loads((control / (mode + '.json')).read_text())
            for relative in withheld:
                shutil.copyfile(HERE / 'workspace' / relative, control / 'workspace' / relative)
            differences = {}
            for relative, expected in snapshot['files_sha256'].items():
                if sha(control / 'workspace' / relative) != expected:
                    allowed = plans[mode].get('roundoff_diagnostics', {}).get(relative)
                    if relative not in generated or not allowed:
                        raise ValueError('regenerated output or source differs: ' + relative)
                    observed, path_changes = normalize_source_paths(
                        json.loads((control / 'workspace' / relative).read_text()),
                        plans[mode].get('identical_source_path_equivalences', {}).get(relative, {}),
                        snapshot['files_sha256'])
                    delta = compare_diagnostics(json.loads((HERE / 'workspace' / relative).read_text()),
                        observed, allowed)
                    differences[relative] = {'regenerated_sha256': sha(control / 'workspace' / relative),
                        'archive_sha256': expected, 'absolute_roundoff_allowance': 5e-13,
                        'differences': delta, 'identical_source_path_changes': path_changes,
                        'all_other_values_exact': True}
            for alias, relative in aliases.items():
                if sha(control / 'workspace' / alias) != snapshot['files_sha256'][relative]:
                    raise ValueError('runtime input materialization changed')
            record['all_source_and_regenerated_output_bytes_match'] = not differences
            record['all_source_bytes_match'] = True
            record['scientific_replay_matches_with_declared_roundoff'] = True
            record['numeric_roundoff_differences'] = differences
            record['outputs_deleted_before_execution'] = generated
            record['runtime_aliases_materialized_from_flat_snapshot'] = len(aliases)
            record['review_freezes_withheld_during_reproduction'] = {
                relative: snapshot['files_sha256'][relative] for relative in withheld}
            records[mode] = record
    verify()
    result = {'schema': 'oph.common_charge.public_replay.v1', 'status': 'PASS_COMMON_CHARGE_SOURCE_AND_TERMINAL_REVIEWS',
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
