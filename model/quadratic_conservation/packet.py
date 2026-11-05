"""Cold replay of the complete quadratic-conservation certificate."""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from . import algebra


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'


def source_bindings():
    paths = {p.relative_to(ROOT).as_posix() for p in PACKAGE.iterdir()
             if p.suffix in ('.py', '.md')}
    paths.add('model/requirements.txt')
    pending, seen = list(paths), set()
    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.')
        candidates = [module.replace('.', '/')+'.py', module.replace('.', '/')+'/__init__.py']
        candidates += ['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for path in candidates:
            if (ROOT/path).is_file() and path not in paths:
                paths.add(path)
                pending.append(path)
    while pending:
        path = pending.pop()
        if path in seen or not path.endswith('.py'):
            continue
        seen.add(path)
        module = path[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/path).read_text(encoding='utf-8'))):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                if node.level:
                    target = importlib.util.resolve_name(target, package)
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def produce():
    pins = source_bindings()
    result = algebra.certificate()
    if pins != source_bindings():
        raise ValueError('quadratic conservation sources changed during replay')
    return {'schema': 'oph.quadratic_conservation.v1', 'issues': [1140, 1142],
            'all_pass': True, 'certificate': result,
            'canonical_physical_premises_discharged': [],
            'physical_comparison_count_added': 0,
            'integrated_gate_accepted': False,
            'source_law_changed': False, 'source_bindings': pins}


def validate(value, expected=None):
    expected = produce() if expected is None else expected
    if canonical(value) != canonical(expected):
        raise ValueError('quadratic conservation receipt differs from exact reconstruction')


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-quadratic-conservation-') as folder:
        path = Path(folder)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(folder)/'pycache'),
                        '-m', 'model.quadratic_conservation.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.load(path)
    if pins != source_bindings() or value['source_bindings'] != pins:
        raise ValueError('quadratic conservation sources changed across cold replay')
    return value


def verify(path=OUTPUT):
    validate(strict_json.load(Path(path)), report())
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    value = produce() if args.worker else report()
    if args.check:
        validate(strict_json.load(args.check), value)
    else:
        (args.worker or args.write).write_text(canonical(value), encoding='utf-8')
    if not args.worker:
        print('Complete quadratic conservation certificate OK')


if __name__ == '__main__':
    main()
