#!/usr/bin/env python3
"""Verify portable source-repair snapshots or replay one bounded packet offline."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
RUNNER_HOME = HERE


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON key: ' + key)
            result[key] = value
        return result

    def number(value):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError('nonfinite JSON number')
        return result

    def constant(value):
        raise ValueError('nonfinite JSON constant: ' + value)

    return json.loads(text, object_pairs_hook=pairs, parse_float=number,
                      parse_constant=constant)


def snapshot(root):
    """Detect deletion, creation, rewriting and replacement during verification."""
    result = {}
    for path in sorted(root.rglob('*')):
        if path.is_symlink():
            raise ValueError('symlink in isolated evidence: ' + str(path))
        if path.is_file():
            stat = path.stat()
            result[str(path.relative_to(root))] = (
                digest(path), stat.st_size, stat.st_mtime_ns, stat.st_ino)
    return result


def output_verdicts(stdout):
    """Reject explicit failed verdicts; substantive checks remain in bound scripts.

    Some preserved programs print prose around pretty-printed JSON. Decode
    complete objects in that stream without treating arbitrary false-valued
    scientific observations (for example a control outcome) as verdicts.
    """
    decoder = json.JSONDecoder()
    found = 0

    def inspect(value, path='output'):
        nonlocal found
        if isinstance(value, dict):
            for key, item in value.items():
                location = path + '.' + key
                if key in ('all_pass', 'all_passed', 'all_gates_pass', 'all_checks_pass'):
                    found += 1
                    if item is not True:
                        raise ValueError('failed or nonboolean verdict: ' + location)
                elif key in ('gates', 'checks'):
                    # A scalar checks count is an archival summary, not a verdict.
                    if key == 'checks' and type(item) is int and item >= 0:
                        pass
                    elif isinstance(item, dict) and item:
                        for name, verdict in item.items():
                            found += 1
                            if verdict is not True:
                                raise ValueError('failed or nonboolean verdict: ' + location + '.' + name)
                    else:
                        raise ValueError('invalid or empty verdict container: ' + location)
                inspect(item, location)
        elif isinstance(value, list):
            for i, item in enumerate(value):
                inspect(item, path + '[' + str(i) + ']')
        elif isinstance(value, str) and path.endswith(('.stdout', '.stderr')):
            found += output_verdicts(value)

    position = 0
    while position < len(stdout):
        start = stdout.find('{', position)
        if start < 0:
            break
        try:
            _, length = decoder.raw_decode(stdout[start:])
        except json.JSONDecodeError:
            position = start + 1
            continue
        inspect(strict_json(stdout[start:start + length]))
        position = start + length
    return found


def run_checked(command, cwd, env, evidence_root):
    before = snapshot(evidence_root)
    process = subprocess.run(command, cwd=cwd, env=env, text=True,
                             capture_output=True, timeout=60)
    print(process.stdout, end='')
    print(process.stderr, end='', file=sys.stderr)
    after = snapshot(evidence_root)
    if before != after:
        changed = sorted(k for k in before.keys() | after.keys()
                         if before.get(k) != after.get(k))
        raise ValueError('verification changed isolated evidence: ' + ', '.join(changed))
    if process.returncode:
        raise RuntimeError('scientific child exited ' + str(process.returncode))
    return output_verdicts(process.stdout) + output_verdicts(process.stderr)


def custody():
    origins = strict_json((HERE / 'ORIGINS.json').read_text())
    counts = {}
    for name, record in origins['packets'].items():
        directory = HERE / name
        manifest = directory / 'MANIFEST.sha256'
        if digest(manifest) != record['manifest_sha256']:
            raise ValueError('manifest changed: ' + name)
        count = 0
        seen = set()
        for line in manifest.read_text().splitlines():
            expected, relative = line.split('  ', 1)
            path = directory / relative
            if relative in seen or path.is_symlink() or not path.resolve().is_relative_to(directory.resolve()):
                raise ValueError('path outside snapshot: ' + relative)
            seen.add(relative)
            if digest(path) != expected:
                raise ValueError('source changed: ' + name + '/' + relative)
            count += 1
        counts[name] = count
    for relative, expected in origins.get('review_files_sha256', {}).items():
        path = HERE / relative
        if not path.resolve().is_relative_to((HERE / 'reviews').resolve()):
            raise ValueError('review path outside archive: ' + relative)
        if digest(path) != expected:
            raise ValueError('review changed: ' + relative)
    failed = strict_json((HERE / 'contact_failed/RECEIPT.json').read_text())
    actual = sorted(k for k, value in failed['checks'].items() if not value)
    if actual != origins['packets']['contact_failed']['expected_failed_checks']:
        raise ValueError('failed predecessor has been changed or hidden')
    return counts


def main():
    global HERE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=HERE,
                        help='Preserved package to verify; supports a separate reviewed entrypoint')
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--packet', choices=('charge', 'rotation', 'contact', 'scalar', 'bridge', 'clock', 'apolar_join', 'pair_loader', 'pair_route', 'flow', 'energy', 'common', 'unknown_total', 'aligned_pair', 'source_local_flow', 'common_v2'))
    parser.add_argument('--independent', action='store_true')
    args = parser.parse_args()
    HERE = args.source_root.resolve()
    counts = custody()
    if args.verify_only:
        print(json.dumps({'custody_pass': True, 'files': counts,
                          'scientific_execution_performed': False}, sort_keys=True))
        return
    if args.packet is None:
        parser.error('choose --packet for scientific replay or --verify-only for custody')
    scripts = {'charge': 'verify.py', 'rotation': 'independent_check.py',
               'contact': 'independent_check.py', 'scalar': 'independent_check.py',
               'bridge': 'verify.py', 'clock': 'verify.py', 'apolar_join': 'headline.py',
               'pair_loader': 'verify.py', 'pair_route': 'verify.py', 'flow': 'headline.py',
               'energy': 'portable_replay.py', 'common': 'headline.py',
               'unknown_total': 'check.py', 'aligned_pair': 'headline.py',
               'source_local_flow': 'portable_replay.py', 'common_v2': 'headline.py'}
    script = scripts[args.packet] if args.independent else (
        'audit.py' if args.packet == 'scalar' else 'headline.py')
    if args.packet == 'flow':
        script = 'workspace/openwave/shared/compact_logarithmic_flow/headline.py'
    if args.packet in ('energy', 'source_local_flow'):
        script = 'portable_replay.py'
    with tempfile.TemporaryDirectory(prefix='source_repair_replay_') as scratch:
        target = Path(scratch) / ('selected_compact_source_join' if args.packet == 'common' else args.packet)
        shutil.copytree(HERE / args.packet, target,
                        ignore=shutil.ignore_patterns('__pycache__'))
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1',
                   PYTHONINTMAXSTRDIGITS='0', PYTHONOPTIMIZE='0',
                   OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        started = time.perf_counter()
        command = [sys.executable, '-B', str(target / script)]
        if args.packet == 'flow':
            command += ['--attempt', 'stage2_002']
        if args.packet in ('energy', 'source_local_flow'):
            command += ['independent' if args.independent else 'headline']
        if args.packet == 'apolar_join' and args.independent:
            command = [sys.executable, '-B', str(HERE / 'reviews/apolar_join/check_join.py'), str(target)]
            script = 'reviews/apolar_join/check_join.py'
        if args.packet == 'common' and args.independent:
            review = Path(scratch) / 'selected_compact_source_join_review_rotation'
            shutil.copytree(HERE / 'reviews/common', review,
                            ignore=shutil.ignore_patterns('__pycache__'))
            command = [sys.executable, '-B', str(review / 'check.py')]
            script = 'reviews/common/check.py'
        if args.packet in ('unknown_total', 'aligned_pair') and args.independent:
            review = Path(scratch) / 'reviews' / args.packet
            shutil.copytree(HERE / 'reviews' / args.packet, review,
                            ignore=shutil.ignore_patterns('__pycache__'))
            if args.packet == 'unknown_total':
                # Preserve the sealed reviewer's two-level relative source path.
                source = Path(scratch) / 'work/2026-11-16_codex_unknown_total_extraction'
                source.parent.mkdir()
                shutil.move(str(target), str(source))
                target = source
                command = [sys.executable, '-B', str(review / 'check.py')]
            else:
                # The sealed ancestor writes its RESULT on every invocation.
                # Use the separate reviewed read-only successor, retaining the
                # ancestor as evidence rather than treating regeneration as replay.
                adapter = RUNNER_HOME / 'verification'
                origin = strict_json((adapter / 'ORIGIN.json').read_text())
                for name, expected in origin['files'].items():
                    if Path(name).name != name or digest(adapter / name) != expected:
                        raise ValueError('verification successor binding changed: ' + name)
                copied = Path(scratch) / 'verification'
                shutil.copytree(adapter, copied)
                common = Path(scratch) / 'common_v2'
                shutil.copytree(HERE / 'common_v2', common,
                                ignore=shutil.ignore_patterns('__pycache__'))
                command = [sys.executable, '-B', str(copied / 'verify_spatial_v2.py'),
                           '--source', str(target), '--common', str(common),
                           '--reference', str(copied / 'RESULT.json')]
            script = 'reviews/' + args.packet + '/check.py'
            if args.packet == 'aligned_pair':
                script = 'verification/verify_spatial_v2.py (read-only successor)'
        if args.packet == 'common_v2' and args.independent:
            review = Path(scratch) / 'common_v2_review'
            shutil.copytree(HERE / 'reviews/common_v2', review,
                            ignore=shutil.ignore_patterns('__pycache__'))
            command = [sys.executable, '-B', str(review / 'check.py'), str(target)]
            script = 'reviews/common_v2/check.py + check_calibration.py'
        verdict_count = run_checked(command, target, env, Path(scratch))
        if args.packet == 'common_v2' and args.independent:
            verdict_count += run_checked(
                [sys.executable, '-B', str(review / 'check_calibration.py'), str(target)],
                target, env, Path(scratch))
        # Recheck the preserved provider, in addition to the actual isolated copy.
        if custody() != counts:
            raise ValueError('provider custody changed during replay')
        print(json.dumps({'packet': args.packet, 'program': script,
                          'exit_code': 0, 'explicit_verdicts_checked': verdict_count,
                          'saved_evidence_unchanged': True,
                          'child_python_optimize': 0,
                          'seconds': time.perf_counter()-started,
                          'review_group': 'Codex same-group'}, sort_keys=True))


if __name__ == '__main__':
    main()
