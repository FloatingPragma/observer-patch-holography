"""Same-group verification-interface audit; no new scientific assay or award."""
import argparse
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROW = HERE.parents[1]
PACKET = ROW/'work/2026-11-16_codex_verification_contract_v3'
INTEGRATION = ROW/'reviews/2026-11-16_codex_resubmission_integration'
ORIGINAL = INTEGRATION/'outgoing_v2/openwave/xperiments/m10_oph/research'
SOURCE = INTEGRATION/'source_export_v2/source'
MANIFEST = 'b25c9c07328cc62fe7852261af837ddf5d33a03d3031dd74b54d56c74ac1dd42'


def need(ok, message):
    if not ok:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def snapshot(folder):
    return {str(p.relative_to(folder)): (sha(p), p.stat().st_size,
            p.stat().st_mtime_ns, p.stat().st_ino)
            for p in folder.rglob('*') if p.is_file()}


class ExpectedChecks(ast.NodeTransformer):
    """Transform ORIGINAL predicates, independently of the producer's reverse map."""
    def visit_Assert(self, node):
        message = node.msg or ast.Constant('Scientific check at original line '+str(node.lineno))
        return ast.Expr(ast.Call(ast.Name('_require', ast.Load()), [node.test, message], []))

    def visit_Expr(self, node):
        if (isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name)
                and node.value.func.id == 'print'):
            need(len(node.value.args) == 1 and not node.value.keywords, 'Unexpected original print')
            return ast.Assign([ast.Name('_verification_payload', ast.Store())], node.value.args[0])
        return self.generic_visit(node)


def is_main(node):
    return isinstance(node, ast.If) and ast.unparse(node.test) == "__name__ == '__main__'"


def normalized(path, original, direct=False):
    tree = ast.parse(path.read_text())
    if direct:
        tree.body = [n for n in tree.body if not is_main(n)]
    if original:
        tree = ExpectedChecks().visit(tree)
    else:
        tree.body = [n for n in tree.body if not
                     (isinstance(n, ast.FunctionDef) and n.name in ('_require', '_verify_or_record'))]
        if not direct:
            need(is_main(tree.body[-1]), 'Missing output-only main guard')
            need(ast.unparse(tree.body[-1].body[0]) == '_verify_or_record(_verification_payload)',
                 'Unexpected output-only main body')
            need(len(tree.body[-1].body) == 1 and not tree.body[-1].orelse, 'Hidden main content')
            tree.body.pop()
    return ast.dump(tree, include_attributes=False)


def pointer(value, path):
    for key in path:
        value = value[key]
    return value


def canonical(value):
    return json.dumps(value, sort_keys=True, allow_nan=False, separators=(',', ':'))


def without_times(value, paths):
    result = copy.deepcopy(value)
    for path in paths:
        need(path[-1] in ('seconds', 'total_seconds'), 'A scientific field is excluded')
        need(type(pointer(result, path)) in (int, float), 'Runtime is not numeric')
        pointer(result, path[:-1])[path[-1]] = 0
    return canonical(result)


def calculate():
    before = snapshot(PACKET)
    need(sha(PACKET/'MANIFEST.sha256') == MANIFEST, 'Wrong reviewed seal')
    entries = (PACKET/'MANIFEST.sha256').read_text().splitlines()
    for line in entries:
        digest, name = line.split('  ', 1)
        need(sha(PACKET/name) == digest, 'Changed producer leaf '+name)
    calculation_rows = []
    for path in sorted((PACKET/'descendants/scripts').glob('*.py')):
        old = ORIGINAL/'scripts'/path.name
        need(normalized(old, True) == normalized(path, False), 'Changed full calculation '+path.name)
        need(not any(isinstance(n, ast.Assert) for n in ast.walk(ast.parse(path.read_text()))),
             'Optimization-sensitive descendant assertion')
        need((ORIGINAL/'data'/(path.stem+'_summary.json')).read_bytes() ==
             (PACKET/'descendants/data'/(path.stem+'_summary.json')).read_bytes(), 'Changed saved summary')
        calculation_rows.append({'file': path.name, 'source': sha(old), 'descendant': sha(path)})
    producer_rows = []
    for item in load(PACKET/'PRODUCER_BUILD.json'):
        old, new = INTEGRATION/item['original'], PACKET/item['descendant']
        need(normalized(old, True, True) == normalized(new, False, True),
             'Changed direct producer calculation '+item['descendant'])
        producer_rows.append({'file': item['descendant'], 'source': sha(old), 'descendant': sha(new)})
    contracts = load(PACKET/'descendants/SOURCE_JOB_CONTRACTS.json')['jobs']
    old_replay = load(INTEGRATION/'SOURCE_EXPORT_CLEAN_REPLAY_V2.json')['result']['current_replay']['jobs']
    replay = load(PACKET/'validation_harness/SOURCE_REPLAY.json')
    need(len(contracts) == len(old_replay) == len(replay['jobs']) == 28, 'Wrong archive coverage')
    runtime_count = verdict_count = 0
    for i, (contract, old, fresh) in enumerate(zip(contracts, old_replay, replay['jobs'])):
        need(contract['index'] == i and old['directory'] == contract['directory'], 'Wrong historical job')
        need(old['command'][2:] == [contract['script'], *contract['arguments']], 'Changed command')
        need(canonical(old['output']) == canonical(contract['expected_output']), 'Baseline was regenerated')
        need(fresh['exit_code'] == 0 and fresh['stderr'] == '' and fresh['passed'], 'Failed retained replay')
        need(without_times(fresh['output'], contract['runtime_paths']) ==
             without_times(old['output'], contract['runtime_paths']), 'Changed replay science')
        for p in contract['required_verdict_paths']:
            need(pointer(fresh['output'], p) in (True, 'PASS'), 'Failed required current verdict')
        runtime_count += len(contract['runtime_paths'])
        verdict_count += len(contract['required_verdict_paths'])
    need(replay['passed'] and replay['saved_evidence_unchanged'], 'Archive no-write audit failed')
    # Every retained compact mutant must fail specifically at its changed predicate.
    lines = dict(charge=132, energy=46, spin=139, common=353, pair=77)
    for name, line in lines.items():
        for mode in ('normal', 'optimized'):
            r = load(PACKET/f'validation_payloads/{name}_mutant_{mode}.json')
            need(r['exit_code'] == 1 and not r.get('timeout') and
                 r['stderr'].strip().endswith('AssertionError: Scientific check at original line '+str(line)),
                 'Wrong retained predicate failure')
    wrappers = load(PACKET/'WRAPPER_AUDIT.json')['results']
    for item in wrappers:
        r = load(PACKET/'wrapper_processes'/(item['label']+'.json'))
        need(r['exit_code'] == 1 and not r.get('timeout') and 'AssertionError' in r['stderr'],
             'Wrong historical predicate rejection')
        if 'mutation_fragment' in item:
            need(item['mutation_fragment'] in r['stderr'] or item['label'] == 'job_02_scientific_predicate',
                 'Trace does not identify intended predicate')
    # Execute interfaces from disposable copies, retaining source/saved metadata.
    processes = []
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONOPTIMIZE='2',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')

    def invoke(label, command, expected, cause=None):
        p = subprocess.run(command, env=env, capture_output=True, text=True, timeout=60)
        need(p.returncode == expected, label+' exit status')
        if cause:
            need(cause in p.stderr, label+' failed for wrong reason')
        processes.append({'label': label, 'exit_code': p.returncode, 'required_stderr': cause})
        return p

    with tempfile.TemporaryDirectory(prefix='independent-verification-') as temp:
        work = Path(temp)
        shutil.copytree(PACKET/'descendants/scripts', work/'scripts')
        shutil.copytree(PACKET/'descendants/data', work/'data')
        saved = snapshot(work/'data')
        sites = {'charge': ('entries == {"native": 120, "swap_hold": 0}', 'entries == {"native": 121, "swap_hold": 0}'),
                 'energy': ('np.max(abs(R @ R - np.eye(33))) == 0', 'np.max(abs(R @ R - np.eye(33))) > 1'),
                 'spin': ('visibility < 0 if law == "native" and turns == 1 else visibility > 0',
                          'visibility > 0 if law == "native" and turns == 1 else visibility < 0'),
                 'common': ('impurity < .002 and pair_impurity < .01', 'impurity < -1 and pair_impurity < -1'),
                 'pair': ('energy[1]<OLD_PAIR', 'energy[1]<1000')}
        for name, (old, new) in sites.items():
            path = work/'scripts'/('m10_4_'+name+'.py')
            text = path.read_text()
            need(text.count(old) == 1, 'Ambiguous scientific mutation')
            path.write_text(text.replace(old, new))
            invoke(name+'_predicate_optimized', [sys.executable, '-B', '-O', str(path)], 1,
                   'AssertionError: Scientific check at original line '+str(lines[name]))
            path.write_text(text)
        charge = work/'scripts/m10_4_charge.py'
        invoke('charge_optimized_clean', [sys.executable, '-B', '-O', str(charge)], 0)
        target = work/'first.json'
        invoke('explicit_first_record', [sys.executable, '-B', str(charge), '--record-to', str(target)], 0)
        target_before = (sha(target), target.stat().st_mtime_ns, target.stat().st_ino)
        invoke('exclusive_overwrite', [sys.executable, '-B', str(charge), '--record-to', str(target)], 1,
               'FileExistsError')
        need(target_before == (sha(target), target.stat().st_mtime_ns, target.stat().st_ino), 'Record overwritten')
        need(saved == snapshot(work/'data'), 'Verification changed saved records')
        # Earlier/nested failures and legitimate historical negatives use actual parser.
        parser = runpy.run_path(str(PACKET/'descendants/verify_source.py'))
        cases = [([{'passed': False}, {'passed': True}], [{'passed': True}, {'passed': True}], False),
                 ([{'passed': True, 'checks': {'inner': False}}], [{'passed': True, 'checks': {'inner': True}}], False),
                 ([{'passed': True, 'history': {'passed': False}, 'external_G6': False}],
                  [{'passed': True, 'history': {'passed': False}, 'external_G6': False}], True)]
        parser_checks = []
        for actual, expected, success in cases:
            try:
                parser['check_output'](actual, {'expected_output': expected, 'runtime_paths': [],
                                              'required_verdict_paths': [[0, 'passed']]})
                passed = True
            except AssertionError:
                passed = False
            need(passed == success, 'Wrong parser disposition')
            parser_checks.append({'actual': actual, 'accepted': passed})
        # A real job exercises clean custody, full expected-output comparison and -O parent.
        r = invoke('historical_job_23_optimized_parent', [sys.executable, '-B', '-O',
                   str(PACKET/'descendants/verify_source.py'), '--source', str(SOURCE), '--job', '23'], 0)
        replay_one = json.loads(r.stdout)
        need(replay_one['passed'] and replay_one['saved_evidence_unchanged'] and
             len(replay_one['jobs']) == 1, 'Focused archive replay failed')
        invoke('saved_adjudication_optimized', [sys.executable, '-B', '-O', str(PACKET/'adjudicate.py')], 0)
    need(before == snapshot(PACKET), 'Sealed producer packet changed')
    return {'passed': True, 'same_group': True, 'external_G6': False,
            'producer_manifest': MANIFEST, 'verified_leaves': len(entries),
            'payload_calculations': calculation_rows, 'direct_producer_calculations': producer_rows,
            'retained_full_archive_jobs': 28, 'runtime_only_paths': runtime_count,
            'required_current_verdict_paths': verdict_count, 'retained_scientific_predicate_mutants': 10,
            'retained_historical_entry_mutants': len(wrappers), 'fresh_processes': processes,
            'parser_cases': parser_checks, 'producer_packet_unchanged': True,
            'historical_focused_replay_unchanged': True, 'saved_copy_bytes_metadata_unchanged': True}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    if args.record_to:
        need(not args.record_to.exists(), 'Exclusive first-record destination exists')
    start = time.monotonic()
    result = calculate()
    if args.record_to:
        with args.record_to.open('x') as stream:
            stream.write(json.dumps(result, indent=2)+'\n')
    else:
        need(result == load(HERE/'AUDIT.json'), 'Saved independent audit differs')
    print(json.dumps({'passed': result['passed'], 'fresh_processes': len(result['fresh_processes']),
                      'seconds': time.monotonic()-start, 'same_group': True, 'external_G6': False}))
