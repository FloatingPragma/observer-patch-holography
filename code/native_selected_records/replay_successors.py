"""Read-only replay of explicitly bound successor experiments.

This full-suite driver is separate from the inspectable compact headlines.
Every child has its own sealed saved-evidence comparison. The driver checks
all emitted scientific verdicts, explicit output fields, source custody and
file metadata. It never creates or replaces scientific evidence.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def snapshot():
    return {str(p.relative_to(HERE)): (hashlib.sha256(p.read_bytes()).hexdigest(),
            p.stat().st_mtime_ns, p.stat().st_ino)
            for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts}


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON key '+key)
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError('nonfinite JSON constant '+value)


def finite_float(value):
    result = float(value)
    require(math.isfinite(result), 'nonfinite JSON number '+value)
    return result


def load_json(path):
    return json.loads(path.read_text(), object_pairs_hook=unique_object,
                      parse_constant=reject_constant, parse_float=finite_float)


def objects(text):
    decoder = json.JSONDecoder(object_pairs_hook=unique_object,
                               parse_constant=reject_constant, parse_float=finite_float)
    result = []
    remaining = text.lstrip()
    while remaining:
        value, end = decoder.raw_decode(remaining)
        require(isinstance(value, dict), 'non-object child output')
        result.append(value)
        remaining = remaining[end:].lstrip()
    require(bool(result), 'empty child output')
    return result


def verdicts(value):
    if isinstance(value, dict):
        for key, child in value.items():
            if key in ('passed', 'gates_pass', 'scientific_gates_pass', 'verified'):
                require(type(child) is bool and child, 'child scientific failure')
                yield child
            yield from verdicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from verdicts(child)


def same_typed(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same_typed(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(same_typed(x, y) for x, y in zip(a, b))
    return a == b


def one(job):
    folder = (HERE/'evidence/openwave'/job['directory']).resolve()
    require(folder.is_relative_to((HERE/'evidence/openwave').resolve()), 'job escapes evidence')
    script = folder/job['script']
    require(hashlib.sha256(script.read_bytes()).hexdigest() == job['script_sha256'], 'unbound replay script')
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONOPTIMIZE='0',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    start = time.perf_counter()
    child = subprocess.run([sys.executable, '-B', job['script'], *job['arguments']],
        cwd=folder, env=env, capture_output=True, text=True, timeout=job['timeout_seconds'])
    require(child.returncode == 0, job['id']+' failed: '+child.stderr[-2000:]+child.stdout[-2000:])
    output = objects(child.stdout)
    require(bool(list(verdicts(output))), 'child emitted no scientific verdict')
    if 'saved_stdout' in job:
        saved = load_json(folder/job['saved_stdout'])
        require(len(output) == 1 and type(saved) is dict and same_typed(output[0], saved),
                'compact output differs from sealed evidence')
    final = output[-1]
    for key, expected in job['expected_final_fields'].items():
        actual = final[key]
        require(type(actual) is type(expected) and actual == expected, 'unexpected final field '+key)
    return {'id': job['id'], 'passed': True, 'exit_status': child.returncode,
            'seconds': time.perf_counter()-start, 'stdout': output, 'stderr': child.stderr}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    args = parser.parse_args()
    import verify
    before = snapshot()
    custody = verify.verify()
    contract_path = HERE/'SUCCESSOR_JOBS.json'
    contracts = load_json(contract_path)
    jobs = [j for j in contracts['jobs'] if args.full or j['headline']]
    require(bool(jobs) and len({j['id'] for j in jobs}) == len(jobs), 'invalid job set')
    started = time.perf_counter()
    def recorded_one(job):
        try:
            return one(job)
        except Exception as error:
            return {'id': job['id'], 'passed': False, 'error': repr(error)}
    with ThreadPoolExecutor(max_workers=2) as pool:
        result = list(pool.map(recorded_one, jobs))
    after_custody = verify.verify()
    require(snapshot() == before, 'replay changed saved bytes or metadata')
    passed = all(job['passed'] is True for job in result)
    print(json.dumps({'passed': passed, 'scope': contracts['scope'], 'jobs': result,
        'seconds': time.perf_counter()-started, 'saved_evidence_unchanged': True,
        'custody_before': custody, 'custody_after': after_custody,
        'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'contracts_sha256': hashlib.sha256(contract_path.read_bytes()).hexdigest(),
        'external_G6': False, 'official_awards_added': 0}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
