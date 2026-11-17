"""The supported native-function scaffold and source dependency graph."""
from __future__ import annotations

import ast
import hashlib
import importlib.util
from pathlib import Path
import re


def require(condition, message):
    if not condition:
        raise ValueError(message)


def native_text(root):
    text = (Path(root)/'model/observer/native.py').read_bytes().decode('utf-8')
    functions = [n for n in ast.parse(text).body
                 if isinstance(n, ast.FunctionDef) and n.name == 'native']
    require(len(functions) == 1, 'source:ambiguous-native-function')
    return ast.get_source_segment(text, functions[0])


def _expression(node, allowed):
    if isinstance(node, ast.Constant):
        return type(node.value) is int and node.value >= 0
    if isinstance(node, ast.Name):
        return node.id in allowed
    if isinstance(node, ast.Subscript):
        return (isinstance(node.value, ast.Name) and node.value.id == 'loads'
                and isinstance(node.slice, ast.Name) and node.slice.id in ('a', 'b'))
    if isinstance(node, ast.UnaryOp):
        return isinstance(node.op, ast.USub) and _expression(node.operand, allowed)
    if isinstance(node, ast.BinOp):
        return (isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Mod, ast.FloorDiv))
                and _expression(node.left, allowed) and _expression(node.right, allowed))
    return (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id == 'abs' and len(node.args) == 1 and not node.keywords
            and _expression(node.args[0], allowed))


def parse_source(text):
    require(type(text) is str and text.isascii() and not any(c in text for c in '\r\t'),
            'grammar:source-character-profile')
    lines = [line for line in text.split('\n') if line]
    require(len(lines) == 8 and lines[0] == 'def native(loads, seam, coin):'
            and lines[1] == '    a, b = seam' and lines[4] == '    out = list(loads)'
            and lines[2].startswith('    d = ') and lines[3].startswith('    current = ')
            and re.fullmatch(r'    out\[a\] [+-]= .+', lines[5]) is not None
            and re.fullmatch(r'    out\[b\] [+-]= .+', lines[6]) is not None
            and lines[7].startswith('    return tuple(out), (') and lines[7].endswith(')'),
            'grammar:source-scaffold')
    try:
        tree = ast.parse(text)
    except SyntaxError as error:
        raise ValueError('grammar:source-expression') from error
    require(len(tree.body) == 1 and isinstance(tree.body[0], ast.FunctionDef), 'grammar:function')
    body = tree.body[0].body
    require(len(body) == 7 and isinstance(body[1], ast.Assign) and isinstance(body[2], ast.Assign)
            and isinstance(body[4], ast.AugAssign) and isinstance(body[5], ast.AugAssign)
            and isinstance(body[6], ast.Return), 'grammar:statement-shape')
    require(_expression(body[1].value, {'coin'})
            and _expression(body[2].value, {'coin', 'd'})
            and _expression(body[4].value, {'coin', 'd', 'current'})
            and _expression(body[5].value, {'coin', 'd', 'current'}), 'grammar:expression-domain')
    result = body[6].value
    require(isinstance(result, ast.Tuple) and len(result.elts) == 2
            and isinstance(result.elts[1], ast.Tuple) and len(result.elts[1].elts) == 2
            and all(_expression(n, {'coin', 'd', 'current'}) for n in result.elts[1].elts),
            'grammar:record-domain')
    return tree


def compile_native(text):
    tree = parse_source(text)
    namespace = {'__builtins__': {'abs': abs, 'list': list, 'tuple': tuple}}
    exec(compile(tree, '<native-source-slice>', 'exec'), namespace)
    return namespace['native']


def grammar_cases(text):
    return {
        'missing_return_boundary': text.replace('return tuple', 'returntuple'),
        'split_identifier': text.replace('abs(d)', 'a bs(d)'),
        'split_floor_token': text.replace('// 2', '/ / 2'),
        'split_update_token': text.replace('-= current', '- = current'),
        'wrong_indentation': text.replace('    out[a]', 'out[a]'),
        'difference_self_reference': text.replace('d = loads[a] - loads[b]', 'd = d + 1'),
        'current_self_reference': text.replace('(d - (2 * coin', '(current - (2 * coin'),
        'leading_zero': text.replace('// 2', '// 02'),
    }


def _python_file(root, module):
    path = root/Path(*module.split('.'))
    if path.with_suffix('.py').is_file():
        return path.with_suffix('.py')
    if (path/'__init__.py').is_file():
        return path/'__init__.py'
    return None


def python_closure(root, seeds):
    root = Path(root).resolve()
    pending, visited, paths = list(seeds), set(), set()
    while pending:
        module = pending.pop()
        if module in visited:
            continue
        visited.add(module)
        path = _python_file(root, module)
        if path is None:
            continue
        paths.add(path)
        package = module if path.name == '__init__.py' else module.rpartition('.')[0]
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'))):
            if isinstance(node, ast.Import):
                pending.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                name = '.'*node.level+(node.module or '')
                base = importlib.util.resolve_name(name, package) if node.level else name
                if base:
                    pending.append(base)
                    pending.extend(base+'.'+alias.name for alias in node.names)
        # Importing a module executes each package initializer too.
        parts = module.split('.')
        for count in range(1, len(parts)):
            parent = '.'.join(parts[:count])
            parent_file = _python_file(root, parent)
            if parent_file and parent_file.name == '__init__.py':
                pending.append(parent)
    return paths


def lean_closure(root, seeds):
    root = Path(root)
    pending, visited, paths = list(seeds), set(), set()
    while pending:
        module = pending.pop()
        if module in visited:
            continue
        visited.add(module)
        path = root/'Lean'/Path(*module.split('.')).with_suffix('.lean')
        if not path.is_file():
            continue
        paths.add(path)
        for line in path.read_text(encoding='utf-8').splitlines():
            if line.startswith('import '):
                pending.extend(line[7:].split('--', 1)[0].split())
    return paths


def source_bindings(root):
    root = Path(root).resolve()
    package = 'model.native_binary_source_refinement'
    paths = python_closure(root, [package+'.__main__', package+'.certificate', package+'.native_validation',
                                 'model.native_quantum_reconstruction.certificate',
                                 'model.native_gravity_response.certificate'])
    paths |= lean_closure(root, ['EventAlgebra.NativeBinarySourceWitnesses',
                                'EventAlgebra.NativeBinarySourceReplay',
                                'EventAlgebra.NativeBinarySourceGrammarCheck'])
    paths |= {root/'Lean/lakefile.lean', root/'Lean/lake-manifest.json', root/'Lean/lean-toolchain'}
    # The carrier is loaded dynamically through slow_geometry.source(). Its
    # public specification supplies the archive roots; follow their imports.
    import json
    spec_path = root/'evidence/native_geometric_source_20260925/spec.json'
    archive = root/'evidence/observer_dynamics_20260925/oph-physics-sim'
    declared = json.loads(spec_path.read_bytes())['source_sha256']
    seeds = []
    for name in declared:
        target = (archive/name).resolve()
        require(target.is_relative_to(archive.resolve()) and target.is_file(), 'source:archive-path')
        require(hashlib.sha256(target.read_bytes()).hexdigest() == declared[name],
                'source:archive-specification-binding')
        paths.add(target)
        if name.endswith('.py'):
            seeds.append('.'.join(Path(name).with_suffix('').parts))
    paths |= python_closure(archive, seeds)
    paths.add(spec_path)
    return {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}
