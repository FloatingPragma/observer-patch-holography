"""Replay retained finite candidates with explicit public-root path adapters.

All scientific producer syntax is unchanged. Original bytes are authenticated;
only the enumerated path constants are adapted in disposable executable copies.
"""
from __future__ import annotations

import argparse
import ast
import copy
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RETAINED = HERE / 'retained'
RETAINED_SHA256 = '264f191e2d8290de65072102e31db0b052d26207c3ecb14fe10c59076eb20ba7'
OLD_ROOT = '/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality'
MOVED_ROOT = '/Users/muellerberndt/Projects/oph-meta/repair-law/reverse-engineering-reality'
PATH_ADAPTERS = {
    'CM/cm_r2_compare.py': {MOVED_ROOT: 1},
    'FR/fr_common.py': {OLD_ROOT: 1},
    'SIM/derive/lnat/lnat_common.py': {OLD_ROOT: 1},
    'SIM/derive/lnat3/lnat3_common.py': {OLD_ROOT: 1, MOVED_ROOT: 1},
}
OUTPUTS = {
    'CM': ['r2.json', 'x8.json', 'r2_compare.json'],
    'DR': ['register.json'],
    'SIM/derive/lnat3': ['lnat3_t1.json', 'lnat3_t2.json', 'lnat3_t3.json', 'lnat3_t4.json'],
}

CHILD = r'''
import hashlib,json,os,pathlib,runpy,sys
program=pathlib.Path(sys.argv[1]).resolve(); audit=pathlib.Path(sys.argv[2]); args=sys.argv[3:]
seen=set()
def hook(event,args):
    if event=='open' and isinstance(args[0],str):
        p=pathlib.Path(args[0]).resolve()
        if p.suffix in {'.py','.json','.md','.npy','.npz'}: seen.add(str(p))
sys.addaudithook(hook)
sys.path.insert(0,str(program.parent));sys.argv=[str(program),*args]
runpy.run_path(str(program),run_name='__main__')
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path and str(path).endswith('.py'):seen.add(str(pathlib.Path(path).resolve()))
paths=sorted(seen)
rows={p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest() for p in paths if pathlib.Path(p).is_file()}
audit.write_text(json.dumps(rows,sort_keys=True))
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def differences(old, new, path=()):
    """Report every leaf difference, without a permissive numerical comparator."""
    if type(old) is not type(new):
        return [{'path': list(path), 'old': old, 'fresh': new}]
    if isinstance(old, float) and math.isnan(old) and math.isnan(new):
        return []  # Recorded separately as undefined source numeric leaves.
    if isinstance(old, dict):
        if set(old) != set(new):
            return [{'path': list(path), 'old_keys': sorted(old), 'fresh_keys': sorted(new)}]
        return sum((differences(old[k], new[k], (*path, k)) for k in sorted(old)), [])
    if isinstance(old, list):
        if len(old) != len(new):
            return [{'path': list(path), 'old_length': len(old), 'fresh_length': len(new)}]
        return sum((differences(a, b, (*path, i)) for i, (a, b) in enumerate(zip(old, new))), [])
    return [] if old == new else [{'path': list(path), 'old': old, 'fresh': new}]


def nonfinite_leaves(value, path=()):
    if isinstance(value, float) and not math.isfinite(value):
        return [{'path': list(path), 'value': str(value)}]
    if isinstance(value, dict):
        return sum((nonfinite_leaves(v, (*path, k)) for k, v in value.items()), [])
    if isinstance(value, list):
        return sum((nonfinite_leaves(v, (*path, i)) for i, v in enumerate(value)), [])
    return []


def project_path_outputs(name, original, fresh, public_root):
    """Normalize only four authenticated execution-location output fields."""
    paths = {
        'CM/r2.json': ('FR_register_check', 'repository_root_used_by_FR_code'),
        'DR/register.json': ('provenance_FR', 'repository_root_used'),
        'SIM/derive/lnat3/lnat3_t1.json': ('RER_path_used',),
        'SIM/derive/lnat3/lnat3_t2.json': ('RER_path_used',),
    }
    before, after = copy.deepcopy(original), copy.deepcopy(fresh)
    records = []
    if name in paths:
        path = paths[name]
        old, new = before, after
        for key in path[:-1]:
            old, new = old[key], new[key]
        if new[path[-1]] != str(public_root):
            raise ValueError('producer escaped the public source root')
        records.append({'path': list(path), 'retained_value': old[path[-1]],
                        'current_value': '<PUBLIC_REPOSITORY_ROOT>'})
        old[path[-1]] = new[path[-1]] = '<PUBLIC_REPOSITORY_ROOT>'
    if name == 'DR/register.json':
        old, new = before['provenance_FR'], after['provenance_FR']
        if old['canonical_checkout_present'] is not False or new['canonical_checkout_present'] is not True:
            raise ValueError('public source presence transition differs')
        records.append({'path': ['provenance_FR', 'canonical_checkout_present'],
                        'retained_value': False, 'current_value': True})
        old['canonical_checkout_present'] = new['canonical_checkout_present'] = True
    return before, after, records


def retained_manifest():
    path = HERE/'retained.json'
    if sha(path) != RETAINED_SHA256:
        raise ValueError('retained source authority differs')
    data = json.loads(path.read_text())
    for name, row in data['files'].items():
        path = HERE/row['retained_path']
        if not path.resolve().is_relative_to(RETAINED.resolve()) or sha(path) != row['sha256']:
            raise ValueError('retained original differs: '+name)
    return data


def adapt_source(name, raw, public_root):
    """An exact AST constant adapter; executable structure is conserved."""
    allowed = PATH_ADAPTERS.get(name, {})
    tree = ast.parse(raw)
    before = ast.dump(tree, include_attributes=False)
    changed = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value in allowed:
            changed.append((node, node.value))
            node.value = str(public_root)
    if {old: sum(value == old for _, value in changed) for old in allowed} != allowed:
        raise ValueError('path metadata literal count changed: '+name)
    transformed = ast.unparse(tree)+'\n'
    for node, old in changed:
        node.value = old
    if ast.dump(tree, include_attributes=False) != before:
        raise ValueError('scientific producer syntax changed')
    return transformed if changed else raw, [old for _, old in changed]


def stage(destination, manifest):
    adaptations = {}
    for name, row in manifest['files'].items():
        if not name.endswith('.py'):
            continue
        target = destination/name
        target.parent.mkdir(parents=True, exist_ok=True)
        original = (HERE/row['retained_path']).read_text()
        adapted, replaced = adapt_source(name, original, ROOT)
        target.write_text(adapted)
        adaptations[name] = {'original_sha256': row['sha256'],
                             'path_constants_replaced': replaced,
                             'replacement': '<PUBLIC_REPOSITORY_ROOT>' if replaced else None,
                             'all_other_AST_nodes_equal': True}
    return adaptations


def execute_group(name, destination):
    jobs = {
        'CM': [('cm_r2.py', []), ('cm_r2_x8.py', []), ('cm_r2_compare.py', [])],
        'DR': [('dr_register.py', ['register.json'])],
        'SIM/derive/lnat3': [(f't{i}_{suffix}.py', []) for i, suffix in
                           [(1, 'source'), (2, 'analyzer'), (3, 'assemble'), (4, 'fr_multitick')]],
    }
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', FR_RER=str(ROOT),
               PYTHONPATH=os.pathsep.join([str(ROOT), str(destination / 'FR')]),
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    custody = {}
    for script, args in jobs[name]:
        cwd = destination / name
        audit = cwd / (script + '.inputs.json')
        command = [sys.executable, '-B', '-c', CHILD, str(cwd / script), str(audit), *args]
        run = subprocess.run(command, cwd=cwd, env=env, text=True, capture_output=True, timeout=360)
        if run.returncode:
            raise RuntimeError(f'{name}/{script}:\n{run.stdout}\n{run.stderr}')
        custody.update(json.loads(audit.read_text()))
    return custody


def run():
    import mpmath, numpy, scipy, sympy
    manifest = retained_manifest()
    explicit = {str((HERE/row['retained_path']).relative_to(ROOT)): row['sha256']
                for row in manifest['files'].values()}
    explicit[str((HERE/'retained.json').relative_to(ROOT))] = RETAINED_SHA256
    for path in (*HERE.glob('*.py'), *HERE.glob('*.md')):
        explicit[str(path.relative_to(ROOT))] = sha(path)
    with tempfile.TemporaryDirectory(prefix='oph-fable-audit-') as td:
        destination = Path(td)
        adaptations = stage(destination, manifest)
        with ThreadPoolExecutor(max_workers=3) as pool:
            records = list(pool.map(lambda name: execute_group(name, destination), OUTPUTS))
        for custody in records:
            for path, digest in custody.items():
                p = Path(path)
                if p.is_relative_to(destination):
                    relative = p.relative_to(destination)
                    # Execute only the authenticated source or its typed
                    # path-constant projection. Generated outputs are below.
                    if p.suffix == '.py':
                        name = str(relative)
                        original = (HERE/manifest['files'][name]['retained_path']).read_text()
                        adapted, _ = adapt_source(name, original, ROOT)
                        assert hashlib.sha256(adapted.encode()).hexdigest() == digest
                elif p.is_relative_to(ROOT) and not p.is_relative_to(ROOT / '.venv'):
                    explicit[str(p.relative_to(ROOT))] = digest
        fresh, comparisons = {}, {}
        for name, outputs in OUTPUTS.items():
            for output in outputs:
                key = name + '/' + output
                raw = (destination / key).read_bytes()
                value = json.loads(raw)
                retained = HERE/manifest['files'][key]['retained_path']
                original = json.loads(retained.read_text())
                fresh[key] = value
                old_projected, fresh_projected, output_projection = project_path_outputs(key, original, value, ROOT)
                comparisons[key] = {'original_sha256': sha(retained),
                                    'fresh_byte_sha256': hashlib.sha256(raw).hexdigest(),
                                    'explicit_path_metadata_output_projection': output_projection,
                                    'retained_nonfinite_numeric_leaves': nonfinite_leaves(original),
                                    'fresh_nonfinite_numeric_leaves': nonfinite_leaves(value),
                                    'leaf_differences': differences(old_projected, fresh_projected)}
        from .exact_controls import controls, summarize
        exact = controls(ROOT, RETAINED)
        for name in exact.pop('input_paths'):
            explicit[str(Path(name).relative_to(ROOT))] = sha(name)
        result = {'schema': 'oph.particle_construction_scope.v1',
                  'scope': 'Retained scientific producer replay with typed path-constant adapters; no physical-status promotion or empirical prediction.',
                  'path_metadata_adaptations': adaptations,
                  'environment': {'python': sys.version.split()[0], 'numpy': numpy.__version__,
                                  'scipy': scipy.__version__, 'sympy': sympy.__version__, 'mpmath': mpmath.__version__},
                  'producer_comparisons': comparisons, 'findings': summarize(fresh, exact),
                  'exact_controls': exact, 'input_sha256': dict(sorted(explicit.items()))}
    for name, digest in explicit.items():
        assert sha(ROOT / name) == digest, 'source changed during replay: ' + name
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write == args.check:
        parser.error('choose exactly one of --write/--check')
    result = run()
    path = HERE / 'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    else:
        saved = json.loads(path.read_text())
        # Timings and fresh JSON byte hashes may vary. Every scientific leaf
        # difference is retained for inspection rather than silently admitted.
        for value in (result, saved):
            for row in value['producer_comparisons'].values():
                row.pop('fresh_byte_sha256')
                row['leaf_differences'] = [d for d in row['leaf_differences']
                    if d['path'] != ['seconds']]
        if canonical(result) != canonical(saved):
            raise ValueError('audit differs; inspect exact source/environment/result differences')
    print('PASS: material Fable replay and exact finite-scope controls')


if __name__ == '__main__':
    main()
