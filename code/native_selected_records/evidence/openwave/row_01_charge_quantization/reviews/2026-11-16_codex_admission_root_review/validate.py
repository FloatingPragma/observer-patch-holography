"""Read-only held-out audit, exact custody, and actual rejection checks."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]/'work'
SOURCES = {
    '2026-11-16_codex_record_only_admission_v2': 'ba4fc9e36bf41592b272feee134c4bdd54d67c51416c8f992416535d536c9f2f',
    '2026-11-16_codex_admission_coverage_amendment': 'f697fd2c0fa1784baf1df69e28c5a6bcc8334530eeb178dbb68cb60a54a936ac',
}


def require(value, why):
    if not value:
        raise ValueError(why)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(paths):
    return {str(p): [sha(p), p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ino] for p in paths}


def run():
    protected = [HERE/'AUDIT.json']
    leaves = {}
    for name, expected in SOURCES.items():
        folder = WORK/name
        manifest = folder/'MANIFEST.sha256'
        require(sha(manifest) == expected, 'source manifest changed')
        lines = manifest.read_text().splitlines()
        for line in lines:
            digest, relative = line.split('  ', 1)
            path = (folder/relative).resolve()
            require(path.is_relative_to(folder.resolve()) and sha(path) == digest, 'unbound source leaf')
            protected.append(path)
        protected.append(manifest)
        leaves[name] = len(lines)
    for filename, expected in json.loads((HERE/'PREREG.json').read_text())['files'].items():
        require(sha(HERE/filename) == expected, 'held-out preregistration changed')
    before = snapshot(protected)
    checks = []

    def call(label, args, success, predicate=None):
        result = subprocess.run([sys.executable, '-B', *args], cwd=HERE, capture_output=True, text=True, timeout=60)
        require((result.returncode == 0) == success, label+' wrong exit')
        require(predicate is None or predicate in result.stderr+result.stdout, label+' wrong rejection')
        checks.append({'case': label, 'exit': result.returncode, 'intended_predicate': predicate})

    call('held_out_default', ['audit.py'], True)
    call('held_out_optimized', ['-O', 'audit.py'], True)
    call('coverage_amendment', [str(WORK/'2026-11-16_codex_admission_coverage_amendment/check.py')], True)
    with tempfile.TemporaryDirectory(prefix='admission_root_review_') as scratch:
        changed = json.loads((HERE/'AUDIT.json').read_text())
        changed['signed_extractions'] += 1
        target = Path(scratch)/'changed.json'
        target.write_text(json.dumps(changed))
        call('changed_frozen_evidence', ['-O', 'audit.py', '--verify', str(target)], False, 'saved audit differs')
        call('exclusive_record_creation', ['audit.py', '--record-to', str(target)], False, 'refusing to overwrite evidence')
    mutant = "import audit; a=audit.load_controller();\ndef corrupt(request,graph):\n r=a(request,graph);r['initial_core'][0]+=1;return r\naudit.load_controller=lambda:corrupt\naudit.run()"
    call('actual_census_guard_optimized', ['-O', '-c', mutant], False, 'hidden initial-state reconstruction')
    require(snapshot(protected) == before, 'verification changed protected evidence')
    return {'passed': True, 'checks': checks, 'sealed_source_leaves': leaves,
            'held_out_cases': 480, 'source_events': 59852, 'saved_bytes_and_metadata_unchanged': True,
            'review_group': 'Codex', 'external_G6': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    if args.record_to:
        require(not args.record_to.exists(), 'refusing to overwrite validation')
    result = run()
    if args.record_to:
        with args.record_to.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    else:
        require(result == json.loads((HERE/'VALIDATION.json').read_text()), 'validation record differs')
    print(json.dumps(result))
