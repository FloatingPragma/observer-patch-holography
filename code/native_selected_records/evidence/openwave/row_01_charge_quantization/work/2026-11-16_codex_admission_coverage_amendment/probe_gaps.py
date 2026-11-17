"""Verification-only counterexamples to the prior independent reader's coverage."""
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/'2026-11-16_codex_record_only_admission_v2'


def need(ok, message):
    if not ok:
        raise AssertionError(message)


def mutations(data):
    first = next(i for i, r in enumerate(data['rows']) if 'consumer_global_columns' in r)
    x = copy.deepcopy(data)
    x['rows'][first]['consumer_global_columns'].pop()
    yield 'drop_column', x, 0
    x = copy.deepcopy(data)
    columns = x['rows'][first]['consumer_global_columns']
    columns[-1] = copy.deepcopy(columns[0])
    yield 'duplicate_column', x, 0
    x = copy.deepcopy(data)
    x['rows'].pop()
    yield 'drop_input_without_replacement', x, 1
    # Select copies by audited vector/current records, not the claimed total label.
    choices = [i for i, r in enumerate(data['rows']) if r['B'] == 2 and r['M'] == 4
               and not r['control'] and sum(r['controller']['initial_core']) == 2]
    need(len(choices) == 2, 'Expected the two public stream histories')
    x = copy.deepcopy(data)
    x['rows'][choices[1]] = copy.deepcopy(x['rows'][choices[0]])
    yield 'duplicate_stream_input', x, 0
    other = next(i for i, r in enumerate(data['rows']) if r['B'] == 2 and r['M'] == 4
                 and not r['control'] and sum(r['controller']['initial_core']) == 3)
    x = copy.deepcopy(data)
    x['rows'][other] = copy.deepcopy(x['rows'][choices[0]])
    yield 'duplicate_other_total_input', x, 0


def calculate():
    evidence = OLD/'records/EVIDENCE.json.gz'
    before = evidence.read_bytes()
    data = json.loads(gzip.decompress(before))
    rows = []
    with tempfile.TemporaryDirectory(prefix='admission-coverage-probe-') as tmp:
        for label, modified, expected in mutations(data):
            path = Path(tmp)/(label+'.json.gz')
            path.write_bytes(gzip.compress(json.dumps(modified).encode(), mtime=0))
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            p = subprocess.run([sys.executable, '-B', str(OLD/'check.py'), '--evidence', str(path)],
                               capture_output=True, text=True, timeout=60)
            need(p.returncode == expected, 'Unexpected original-reader result: '+label)
            need(hashlib.sha256(path.read_bytes()).hexdigest() == digest, 'Mutation was overwritten')
            result = json.loads(p.stdout)
            result.pop('seconds', None)
            rows.append({'mutation': label, 'original_exit': p.returncode,
                         'original_result': result, 'saved_mutant_unchanged': True})
    need(evidence.read_bytes() == before, 'Original evidence changed')
    return {'passed': True, 'verification_only': True, 'external_G6': False,
            'gap_demonstrated': True, 'source_evidence_sha256': hashlib.sha256(before).hexdigest(),
            'mutations': rows, 'original_evidence_unchanged': True}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    if args.record_to:
        need(not args.record_to.exists(), 'First-record destination exists')
    start = time.perf_counter()
    result = calculate()
    if args.record_to:
        with args.record_to.open('x') as stream:
            stream.write(json.dumps(result, indent=2)+'\n')
    else:
        need(result == json.loads((HERE/'GAP_PROBE.json').read_text()), 'Saved gap probe differs')
    print(json.dumps({'passed': True, 'cases': len(result['mutations']), 'seconds': time.perf_counter()-start}))
