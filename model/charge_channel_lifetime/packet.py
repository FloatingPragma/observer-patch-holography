"""Full replay of complete source-charge retention and detector probabilities."""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
from tools import strict_json
from . import survival

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
PROOFS = ('model/native_ordered_ground/README.md',
          'model/protected_defect_catalog/PROOF.md',
          'model/native_spin_observer/README.md',
          'model/defect_localization_readout/PROOF.md',
          'model/protected_spin_limit/README.md',
          'model/requirements.txt',
          'code/source_selection_model/response.json',
          'code/m1_source_realization/receipt.json',
          'evidence/native_geometric_source_20260925/spec.json')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    legacy = 'evidence/observer_dynamics_20260925/oph-physics-sim/'
    names = set(PROOFS)|{survival.WITNESS,
                        legacy+'oph_exact/carrier.py', legacy+'oph_exact/__init__.py'}
    names.update(path.relative_to(ROOT).as_posix() for path in HERE.iterdir()
                 if path.suffix in ('.py', '.md'))
    pending = list(names)
    visited = set()

    def add(module):
        parts = module.split('.')
        prefix = legacy if parts[0] in ('oph_exact', 'oph_fpe') else ''
        stem = prefix+'/'.join(parts)
        candidates = [stem+'.py', stem+'/__init__.py']
        candidates += [prefix+'/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)
    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removeprefix(legacy)[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def build():
    before = bindings()
    survival.source_data.cache_clear()
    value = {'schema':'oph.charge_channel_lifetime.v1',
             'source_bindings':before, 'complete_readout':survival.build()}
    if bindings()!=before:
        raise RuntimeError('source changed during retention certificate')
    return value


def validate(value):
    if canonical(value)!=canonical(build()):
        raise ValueError('retention receipt, full outcomes, source or physical scope differs')
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args()
    if args.write:
        args.write.write_text(canonical(build()))
    else:
        validate(strict_json.load(args.check))
    print('PASS: complete conserved-charge contrast and same-clock retention')
