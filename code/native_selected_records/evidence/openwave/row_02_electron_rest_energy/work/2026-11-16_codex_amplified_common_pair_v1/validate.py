"""Actual rejection and evidence-preservation checks for the common join."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def require(value, why):
    if not value:
        raise ValueError(why)


def identity(path):
    stat = path.stat()
    return [hashlib.sha256(path.read_bytes()).hexdigest(), stat.st_size, stat.st_mtime_ns, stat.st_ino]


def run():
    saved = HERE/'SOURCE_JOIN.json.gz'
    before = identity(saved)
    checks = []
    def call(name, args, success, predicate=None):
        result = subprocess.run([sys.executable, '-B', *args], cwd=HERE,
                                capture_output=True, text=True, timeout=60)
        require((result.returncode == 0) == success, name+' exit mismatch')
        require(predicate is None or predicate in result.stdout+result.stderr, name+' wrong rejection')
        checks.append({'name': name, 'exit': result.returncode, 'predicate': predicate})
    call('optimized_saved_replay', ['-O', 'join.py'], True)
    code = ("import join\noriginal=join.native\n"
            "def corrupt(q,e,c):\n out,record=original(q,e,c);record[3]+=1;return out,record\n"
            "join.native=corrupt\njoin.run()")
    call('actual_global_source_guard_optimized', ['-O', '-c', code], False,
         'global current differs from actual calibration')
    with tempfile.TemporaryDirectory(prefix='common_pair_mutation_') as scratch:
        with gzip.open(saved, 'rt') as stream:
            altered = json.load(stream)
        altered['literal_continuing_source_calls'] += 1
        path = Path(scratch)/'changed.json.gz'
        with gzip.open(path, 'wt') as stream:
            json.dump(altered, stream)
        call('changed_saved_verdict', ['-O', 'join.py', '--verify', str(path)], False,
             'saved common evidence differs')
        call('exclusive_first_record', ['join.py', '--record-to', str(path)], False,
             'refusing to overwrite common evidence')
    require(identity(saved) == before, 'saved source bytes or metadata changed')
    return {'passed': True, 'checks': checks, 'saved_evidence_unchanged': True,
            'source_sha256': before[0], 'mutation_scope': 'Integrity and actual scientific guard; native competing-law evidence is separate.',
            'external_G6': False}


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
        require(result == json.loads((HERE/'VALIDATION.json').read_text()), 'validation changed')
    print(json.dumps(result))
