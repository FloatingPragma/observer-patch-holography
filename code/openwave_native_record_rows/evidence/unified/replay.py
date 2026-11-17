#!/usr/bin/env python3
"""Replay sealed ancestors and the unified candidate in isolated directories."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--custody-only', action='store_true')
    args = parser.parse_args()
    pins = json.loads((ROOT / 'SOURCE_PINS.json').read_text())
    for relative, expected in pins.items():
        if digest(ROOT / relative) != expected:
            raise AssertionError(relative)
    origins = json.loads((ROOT / 'ORIGINS.json').read_text())
    for name, origin in origins.items():
        base = ROOT / 'evidence' / name
        manifest = base / 'MANIFEST.sha256'
        assert digest(manifest) == origin['manifest_sha256'], name
        for line in manifest.read_text().splitlines():
            expected, relative = line.split(maxsplit=1)
            assert digest(base / relative.lstrip('*')) == expected, (name, relative)
    report = dict(candidate='MFR-unified-native-field-v1',
                  official='untested', copied_evidence_leaves=len(pins),
                  custody='pass', new_scientific_thresholds=False)
    if args.custody_only:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='oph_unified_submission_') as directory:
        work = Path(directory)
        # Every executable package is self-contained. No workspace import is used.
        for name in ('free', 'orbital', 'unified', 'boost'):
            shutil.copytree(ROOT / 'evidence' / name, work / name)
        results = {}
        for name, script in (('free', 'scripts/task_common.py'),
                             ('orbital', 'scripts/task_orbital.py'),
                             ('unified', 'headline.py'),
                             ('boost', 'headline.py')):
            output = subprocess.run([sys.executable, str(work / name / script)],
                                    cwd=work / name, text=True, check=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            results[name] = json.loads(output.stdout)
        comparisons = {}
        for name, files in {'unified': ('RECEIPT.json', 'CONTROL_RECORDS.json.gz'),
                            'boost': ('RECEIPT.json',)}.items():
            for filename in files:
                key = f'{name}/{filename}'
                comparisons[key] = digest(work / key) == digest(ROOT / 'evidence' / key)
        assert all(comparisons.values()), comparisons
        unified = json.loads((work / 'unified/RECEIPT.json').read_text())
        assert unified['all_registered_checks'] == 'pass'
        assert unified['total_native_calls'] == 492837
        assert results['free']['actual_native_calls'] == 169364
        assert results['orbital']['actual_native_calls'] == 347533
        assert results['boost']['all_registered_analysis_checks'] == 'pass'
        report.update(replays=results, byte_replay=comparisons,
                      actual_native_calls_in_this_bundle_replay=1009734,
                      unified_candidate_native_calls=492837,
                      ancestor_replay_native_calls=516897,
                      reused_ancestor_histories_are_independent_evidence=False,
                      boost_new_native_calls=0,
                      elapsed_seconds=time.monotonic()-started)
    # Inputs must remain sealed even if a child script writes its outputs.
    for relative, expected in pins.items():
        assert digest(ROOT / relative) == expected, relative
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
