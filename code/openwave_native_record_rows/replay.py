#!/usr/bin/env python3
"""Verify frozen native-record evidence and replay selected complete sources."""
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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(root=ROOT):
    pins = json.loads((root / 'SNAPSHOT_MANIFEST.json').read_text())
    for relative, expected in pins.items():
        path = root / relative
        if path.is_symlink() or sha(path) != expected:
            raise ValueError('Changed native-record evidence: ' + relative)
    origins = json.loads((root / 'ORIGINS.json').read_text())
    for name, entry in origins.items():
        base = root / 'evidence' / name
        manifest = base / 'MANIFEST.sha256'
        if sha(manifest) != entry['manifest_sha256']:
            raise ValueError('Changed original manifest: ' + name)
        for line in manifest.read_text().splitlines():
            expected, relative = line.split(maxsplit=1)
            if sha(base / relative.lstrip('*')) != expected:
                raise ValueError('Changed original leaf: ' + name + '/' + relative)
    return len(pins)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('custody', 'representative'), default='custody')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    result = dict(custody='pass', evidence_leaves=verify(), mode=args.mode,
                  official_awards_added=0, executions={}, drafts_executed=False)
    if args.mode == 'representative':
        with tempfile.TemporaryDirectory(prefix='oph_native_record_rows_') as td:
            work = Path(td)
            for name in ('unified', 'confinement', 'fixed_charge'):
                shutil.copytree(ROOT / 'evidence' / name, work / name)
            commands = {
                'unified': ('replay.py',),
                'confinement': ('check.py',),
                'fixed_charge': ('headline.py', 'verify_records.py'),
            }
            # This reader intentionally creates its output exclusively.
            (work / 'confinement/SELF_CHECK.json').unlink()
            outputs = {'unified': (), 'confinement': ('SELF_CHECK.json',),
                       'fixed_charge': ('RAW_RECORDS.jsonl.gz', 'RECEIPT.json', 'VERIFICATION.json')}
            for name, scripts in commands.items():
                began = time.monotonic()
                reports = []
                for script in scripts:
                    proc = subprocess.run([sys.executable, '-B', str(work / name / script)],
                                          cwd=work / name, check=True, text=True,
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    reports.append(json.loads(proc.stdout))
                matches = {f: sha(work / name / f) == sha(ROOT / 'evidence' / name / f)
                           for f in outputs[name]}
                if not all(matches.values()):
                    raise AssertionError((name, matches))
                result['executions'][name] = dict(reports=reports, byte_replay=matches,
                                                  seconds=time.monotonic()-began)
        # The fixed-charge producer and reader each call the native primitive.
        result['actual_native_calls_in_replay'] = 1009734 + 663600 + 2 * 2454480
        result['replayed_ancestry_is_independent_evidence'] = False
    verify()
    result['seconds'] = time.monotonic()-started
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
