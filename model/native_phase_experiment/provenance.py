"""Resolve the experiment's local import closure and declared source data."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = 'evidence/observer_dynamics_20260925/oph-physics-sim'
SPEC = 'evidence/native_geometric_source_20260925/spec.json'
DATA = ('code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json', SPEC)
CLASSIFIER = 'Lean/ObservableNormalForms/tools/verify_protected_obstruction_models.py'


def source_graph():
    """Follow Python imports, parent packages and the pinned geometry manifest."""
    paths, pending, imports = set(DATA), [], {}
    modules = ('model.native_phase_experiment.experiment', 'model.native_phase_experiment.independent',
               'model.native_phase_experiment.verify', 'model.native_phase_experiment.lineage',
               'model.native_phase_experiment.__main__')

    def add_path(path, module):
        path = path.resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError('source dependency leaves the repository')
        name = path.relative_to(ROOT).as_posix()
        if name not in paths:
            paths.add(name)
            pending.append((name, module))
        return name

    def add_module(module):
        prefix = module.split('.')[0]
        if prefix in ('model', 'tools'):
            base = ROOT
        elif prefix in ('oph_exact', 'oph_fpe'):
            base = ROOT / ARCHIVE
        else:
            return []
        parts = module.split('.')
        added = []
        for length in range(1, len(parts) + 1):
            package = base.joinpath(*parts[:length], '__init__.py')
            if package.is_file():
                added.append(add_path(package, '.'.join(parts[:length]) + '.__init__'))
        source = base.joinpath(*parts).with_suffix('.py')
        if source.is_file():
            added.append(add_path(source, module))
        return added

    for module in modules:
        add_module(module)
    add_path(ROOT / CLASSIFIER, 'oph_protected_classifier')
    spec = json.loads((ROOT / SPEC).read_text(encoding='utf-8'))
    for name in spec['source_sha256']:
        if name.endswith('.py'):
            add_module(name[:-3].replace('/', '.'))
    while pending:
        name, module = pending.pop()
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        linked = set()
        for node in ast.walk(ast.parse((ROOT / name).read_text(encoding='utf-8'))):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    linked.update(add_module(alias.name))
            elif isinstance(node, ast.ImportFrom):
                relative = '.' * node.level + (node.module or '')
                target = importlib.util.resolve_name(relative, package) if node.level else relative
                linked.update(add_module(target))
                for alias in node.names:
                    linked.update(add_module(target + '.' + alias.name))
        imports[name] = sorted(linked)
    return {'files': sorted(paths), 'imports': imports, 'declared_data': list(DATA),
            'geometry_manifest': SPEC}


def source_hashes():
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in source_graph()['files']}
