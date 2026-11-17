"""Project import closure and immutable data ancestry for the field program.

The AST walk follows the repository's completion-program source-custody
algorithm, including imports inside functions and package initializers.
Frozen comparison inputs provide conservative coverage of dynamic archive
imports and source data; including their bytes does not transfer their laws.
"""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
try:
    from tools.pin_scope import bound, bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound, bound_pins

ROOT=Path(__file__).resolve().parents[2]
BASE='model/particle_calibration_selection/results.json'


def import_closure(seeds):
    """Pinned repository paths: the seeds and their Python import closure,
    without prose and test modules (tools/pin_scope.py), which are neither
    walked nor returned."""
    paths=set(bound(seeds));pending=list(paths);visited=set()
    def add_module(name):
        stem=name.replace('.','/')
        parts=stem.split('/')
        candidates=[stem+'.py',stem+'/__init__.py']
        candidates += ['/'.join(parts[:end])+'/__init__.py' for end in range(1,len(parts))]
        for candidate in candidates:
            if (ROOT/candidate).is_file() and candidate not in paths:
                paths.add(candidate);pending.append(candidate)
    while pending:
        name=pending.pop()
        if name in visited or not name.endswith('.py'):continue
        visited.add(name)
        module=name.removesuffix('.py').replace('/','.')
        package=module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add_module(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:add_module(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target='.'*node.level+(node.module or '')
                if node.level:target=importlib.util.resolve_name(target,package)
                add_module(target)
                for alias in node.names:add_module(target+'.'+alias.name)
    return set(bound(paths))


def definition_inputs(seeds):
    inherited=bound_pins(json.loads((ROOT/BASE).read_text())['input_sha256'])
    return tuple(sorted(import_closure(set(seeds)|set(inherited)|{BASE,
        'model/particle_common_contract/provenance.py','model/requirements.txt'})))


def definition_hashes(paths):
    hashes=bound_pins({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    inherited=bound_pins(json.loads((ROOT/BASE).read_text())['input_sha256'])
    for path,expected in inherited.items():
        if hashes.get(path)!=expected:raise ValueError('immutable source ancestry changed: '+path)
    return hashes
