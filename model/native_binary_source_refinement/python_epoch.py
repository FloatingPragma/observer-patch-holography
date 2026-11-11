"""Bind receiving Python modules and their caches to one source epoch."""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from functools import _lru_cache_wrapper, wraps
import builtins
import hashlib
import importlib
from importlib.abc import MetaPathFinder
from importlib.machinery import PathFinder, SourceFileLoader
from importlib.util import resolve_name, spec_from_file_location
from pathlib import Path
import sys
from types import CodeType, FunctionType, MappingProxyType, ModuleType


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# The receiving adapter and this helper are the interpreter bootstrap. This
# value belongs to the helper's import, rather than its first verification.
_IMPORTED_HELPER = (Path(__file__).resolve(), _digest(__file__))
_ARCHIVE = 'evidence/observer_dynamics_20260925/oph-physics-sim/'
_SPEC = 'evidence/native_geometric_source_20260925/spec.json'
_MISSING = object()


class _EpochLoader(SourceFileLoader):
    def __init__(self, name, path, epoch):
        super().__init__(name, str(path))
        self.epoch = epoch
        self.raw = None

    def get_code(self, fullname):
        # Compile the bytes enrolled for this import; do not adopt a cached
        # module from another receiving context as a source witness.
        self.raw = Path(self.path).read_bytes()
        self.epoch._accept_import(self.name, Path(self.path), self.raw)
        code = compile(self.raw, self.path, 'exec', dont_inherit=True)
        self.epoch._remember_code(self.name, code)
        return code

    def exec_module(self, module):
        module.__dict__['__builtins__'] = self.epoch._builtins
        self.epoch._record_context(module, self.name, Path(self.path))
        self.epoch._executing[self.name] = module
        try:
            super().exec_module(module)
            if self.raw is None or Path(self.path).read_bytes() != self.raw:
                raise ValueError('source:python-changed-during-import')
            self.epoch._enrolled[self.name] = (module, Path(self.path).resolve(),
                                              hashlib.sha256(self.raw).hexdigest())
            self.epoch._capture_bindings(module, module.__dict__, self.name)
        finally:
            self.epoch._executing.pop(self.name, None)


class _EpochFinder(MetaPathFinder):
    def __init__(self, epoch):
        self.epoch = epoch

    def find_spec(self, fullname, path=None, target=None):
        if not self.epoch._imports_active.get():
            return None
        spec = PathFinder.find_spec(fullname, path, target)
        if (spec is None or not spec.origin
                or not isinstance(spec.loader, SourceFileLoader)):
            return None
        source = Path(spec.origin).resolve()
        # The active interpreter can live in the repository's ignored temp
        # directory. Its installed dependencies belong to the runtime TCB.
        if any(source.is_relative_to(Path(prefix).resolve())
               for prefix in (sys.prefix, sys.base_prefix)):
            return None
        if not source.is_relative_to(self.epoch.root):
            return None
        if target is not None:
            previous = self.epoch._enrolled.get(fullname)
            if previous is None or previous[0] is not target:
                raise ValueError('source:unbound-python-context')
        return spec_from_file_location(fullname, source,
            loader=_EpochLoader(fullname, source, self.epoch),
            submodule_search_locations=spec.submodule_search_locations)


class PythonSourceEpoch:
    """An immutable source context; unknown preloaded producers fail at use."""

    def __init__(self, root, owner):
        self.root = Path(root).resolve()
        self._bindings = None
        self._enrolled = {}
        self._codes = {}
        self._captured = {}
        self._references = []
        self._cache_targets = {}
        self._local_modules = ()
        self._module_names = set()
        self._executing = {}
        self._import_records = {}
        self._possible_modules = {}
        self._imports_active = ContextVar('native_source_imports', default=False)
        self._original_import = builtins.__import__
        self._original_import_module = importlib.import_module
        self._builtins = dict(vars(builtins), __import__=self._source_import)
        self._importlib = ModuleType(importlib.__name__)
        self._importlib.__dict__.update(vars(importlib))
        self._importlib.import_module = self._source_import_module
        self._module_contexts = {}
        self._namespace = None
        self._problem = None
        self._bootstrap = {}
        names = (__name__, owner, 'model', 'model.native_binary_source_refinement', '__main__')
        for name in names:
            module = sys.modules.get(name)
            path = getattr(module, '__file__', None)
            if path is None:
                continue
            path = Path(path).resolve()
            if not path.is_relative_to(self.root):
                continue
            # Only the declared adapter entry can be a receiving __main__.
            if name == '__main__' and path != self.root/'model/native_binary_source_refinement/__main__.py':
                continue
            value = _IMPORTED_HELPER[1] if name == __name__ else _digest(path)
            self._bootstrap[name] = (module, path, value)
            self._record_context(module, name, path, bootstrap=True)
        # This is the sole non-Python/non-Lean input in the source profile.
        # Read it before consumer imports can initialize data or caches.
        self._data = {_SPEC: _digest(self.root/_SPEC)}

    @contextmanager
    def imports(self):
        if self._imports_active.get():
            yield
            return
        finder = _EpochFinder(self)
        token = self._imports_active.set(True)
        sys.meta_path.insert(0, finder)
        try:
            yield
        finally:
            self._imports_active.reset(token)
            sys.meta_path.remove(finder)

    def _possible_module(self, name):
        if self._bindings is not None:
            return name in self._module_names
        if name not in self._possible_modules:
            relative = Path(*name.split('.'))
            self._possible_modules[name] = any(
                (base/relative).with_suffix('.py').is_file()
                or (base/relative/'__init__.py').is_file()
                for base in (self.root, self.root/_ARCHIVE))
        return self._possible_modules[name]

    def _record_import(self, name, module):
        if not self._possible_module(name):
            return
        # Other consumers use ordinary imports. The receiving check still sees
        # the objects they acquired, including imports that populated a cache.
        self._import_records[name, id(module)] = module
        if self._bindings is not None and self._imports_active.get():
            self._require_import(name, module, during_import=True)

    def _record_context(self, module, name, path, bootstrap=False):
        package = (module.__package__ if bootstrap else
                   name if Path(path).name == '__init__.py' else name.rpartition('.')[0])
        self._module_contexts[id(module.__dict__)] = (module, name, package)

    def _import_package(self, namespace):
        context = self._module_contexts.get(id(namespace))
        if context is None or context[0].__dict__ is not namespace:
            return namespace.get('__package__')
        _, name, package = context
        if (type(namespace.get('__name__')) is not str or namespace['__name__'] != name
                or type(namespace.get('__package__')) is not str
                or namespace['__package__'] != package):
            raise ValueError('source:python-import-context')
        return package

    def _require_import(self, name, module, during_import=False):
        record = self._enrolled.get(name) or self._bootstrap.get(name)
        actual = record[0] if record else None
        if actual is None and during_import:
            actual = self._executing.get(name)
        if actual is not module:
            raise ValueError('source:unbound-python-context')

    def _source_import(self, name, globals=None, locals=None, fromlist=(), level=0):
        package = self._import_package(globals) if globals is not None else None
        absolute = resolve_name('.'*level+name, package) if level else name
        module = self._original_import(name, globals, locals, fromlist, level)
        if fromlist:
            self._record_import(absolute, module)
            for item in fromlist:
                child = absolute+'.'+item
                if item != '*' and self._possible_module(child):
                    value = getattr(module, item, _MISSING)
                    if value is _MISSING:
                        # IMPORT_FROM can use the canonical module entry when
                        # a preloaded child is absent from its package attrs.
                        value = sys.modules.get(child)
                    if isinstance(value, ModuleType):
                        self._record_import(child, value)
        elif self._possible_module(absolute):
            value = module
            for item in absolute.split('.')[1:]:
                value = getattr(value, item, None)
            self._record_import(absolute, value)
        return self._importlib if module is importlib else module

    def _source_import_module(self, name, package=None):
        self._import_package(sys._getframe(1).f_globals)
        absolute = resolve_name(name, package) if name.startswith('.') else name
        module = self._original_import_module(name, package)
        self._record_import(absolute, module)
        return module

    def _accept_import(self, name, path, raw):
        if self._bindings is not None:
            relative = path.resolve().relative_to(self.root).as_posix()
            if self._bindings.get(relative) != hashlib.sha256(raw).hexdigest():
                raise ValueError('source:python-import-outside-epoch')
        previous = self._enrolled.get(name)
        if previous is not None and sys.modules.get(name) is not previous[0]:
            raise ValueError('source:python-module-replaced')

    def _remember_code(self, name, code):
        # Use identity, not CodeType's structural equality. An earlier import
        # of identical bytes still does not belong to this receiving epoch.
        self._codes[id(code)] = (code, name)
        for value in code.co_consts:
            if isinstance(value, CodeType):
                self._remember_code(name, value)

    def _capture_bindings(self, owner, namespace, module_name):
        if id(owner) in self._captured:
            return
        self._captured[id(owner)] = owner
        for name, value in tuple(namespace.items()):
            if (isinstance(value, (ModuleType, staticmethod, classmethod, property))
                    or callable(value)):
                self._references.append((namespace, name, value))
            if isinstance(value, _lru_cache_wrapper):
                # This association is captured when the source import creates
                # the wrapper, before callers can change inspection metadata.
                self._cache_targets[id(value)] = (value, value.__wrapped__)
            if isinstance(value, type) and value.__module__ == module_name:
                self._capture_bindings(value, value.__dict__, module_name)

    def _check_references(self):
        contexts = {**self._bootstrap, **self._enrolled}
        modules = {id(row[0]): (name, row) for name, row in contexts.items()}
        namespaces = {id(row[0].__dict__): (name, row) for name, row in contexts.items()}
        seen = set()
        paths = {}

        def local(path):
            if path is None:
                return False
            key = str(path)
            if key in paths:
                return paths[key]
            path = Path(path).resolve()
            result = (path.is_relative_to(self.root)
                      and path.relative_to(self.root).as_posix() in self._bindings)
            paths[key] = result
            return result

        def local_code(value):
            if isinstance(value, FunctionType):
                return local(value.__code__.co_filename)
            if isinstance(value, _lru_cache_wrapper):
                record = self._cache_targets.get(id(value))
                target = record[1] if record else getattr(value, '__wrapped__', None)
                # This only identifies a source-owned reference to reject;
                # mutable metadata never enrolls an unknown wrapper.
                return isinstance(target, FunctionType) and local(target.__code__.co_filename)
            if isinstance(value, (staticmethod, classmethod)):
                return local_code(value.__func__)
            if isinstance(value, property):
                return any(local_code(f) for f in (value.fget, value.fset, value.fdel))
            if isinstance(value, type):
                return any(local_code(f) for f in value.__dict__.values()
                    if isinstance(f, (FunctionType, _lru_cache_wrapper,
                                      staticmethod, classmethod, property)))
            return False

        def namespace(values):
            for key, value in tuple(values.items()):
                if not key.startswith('__'):
                    visit(value)

        def visit(value):
            if id(value) in seen:
                return
            seen.add(id(value))
            if isinstance(value, ModuleType):
                record = modules.get(id(value))
                if record is None or record[1][0] is not value:
                    if (local(getattr(value, '__file__', None))
                            or value.__name__ in self._module_names
                            or any(local_code(v) for v in value.__dict__.values())):
                        raise ValueError('source:unbound-python-context')
                    return
                namespace(value.__dict__)
            elif isinstance(value, FunctionType):
                values = value.__globals__
                record = namespaces.get(id(values))
                if record is None or record[1][0].__dict__ is not values:
                    if local(values.get('__file__')) or local_code(value):
                        raise ValueError('source:unbound-python-context')
                    return
                name = record[0]
                if name not in self._bootstrap:
                    code = self._codes.get(id(value.__code__))
                    if code is None or code[0] is not value.__code__ or code[1] != name:
                        raise ValueError('source:unbound-python-function')
                namespace(values)
                visit(value.__defaults__)
                visit(value.__kwdefaults__)
                for cell in value.__closure__ or ():
                    try:
                        captured = cell.cell_contents
                    except ValueError:
                        continue
                    visit(captured)
            elif isinstance(value, _lru_cache_wrapper):
                record = self._cache_targets.get(id(value))
                if record is None or record[0] is not value:
                    raise ValueError('source:unbound-python-cache')
                visit(record[1])
            elif isinstance(value, (staticmethod, classmethod)):
                visit(value.__func__)
            elif isinstance(value, property):
                for function in (value.fget, value.fset, value.fdel):
                    visit(function)
            elif isinstance(value, type):
                namespace(value.__dict__)
            elif isinstance(value, (tuple, list)):
                for item in value:
                    visit(item)
            elif isinstance(value, dict):
                for item in tuple(value.values()):
                    visit(item)

        namespace(self._namespace)

    def seal(self, bindings, namespace, local_modules):
        if self._bindings is not None:
            raise ValueError('source:python-epoch-already-sealed')
        values = {path: value for path, value in bindings.items() if not path.startswith('Lean/')}
        data = {path: value for path, value in values.items() if not path.endswith('.py')}
        if data != self._data:
            self._problem = 'source:python-data-epoch'
        for _, path, value in (*self._bootstrap.values(), *self._enrolled.values()):
            relative = path.relative_to(self.root).as_posix()
            if values.get(relative) != value:
                self._problem = 'source:python-import-epoch'
        self._bindings = MappingProxyType(values)
        self._namespace = namespace
        self._local_modules = tuple(local_modules)
        self._module_names = {
            path.removeprefix(_ARCHIVE).removesuffix('.py').replace('/', '.').removesuffix('.__init__')
            for path in values if path.endswith('.py')}
        for module in self._local_modules:
            if not any(row[0] is module for row in self._enrolled.values()):
                self._problem = 'source:unbound-python-context'
        for name, (module, _, _) in self._bootstrap.items():
            self._capture_bindings(module, module.__dict__, name)
        try:
            for (name, _), module in self._import_records.items():
                if name in self._module_names:
                    self._require_import(name, module)
            self._check_references()
        except ValueError as error:
            # Keep the invalid incarnation even if registry names disappear
            # later. Registry edits cannot enroll an already held producer.
            self._problem = str(error)

    def check(self, root):
        if Path(root).resolve() != self.root:
            raise ValueError('source:python-epoch-root')
        if self._bindings is None:
            raise ValueError('source:python-epoch-unsealed')
        if self._problem is not None:
            raise ValueError(self._problem)
        for module, _, _ in self._module_contexts.values():
            self._import_package(module.__dict__)
        for (name, _), module in tuple(self._import_records.items()):
            if name in self._module_names:
                self._require_import(name, module)
        for relative, expected in self._bindings.items():
            path = self.root/relative
            if not path.is_file() or _digest(path) != expected:
                raise ValueError('source:python-source-epoch')
        contexts = {**self._bootstrap, **self._enrolled}
        for name, (module, _, _) in contexts.items():
            if sys.modules.get(name) is not module:
                raise ValueError('source:python-module-replaced')
        for namespace, name, expected in self._references:
            if namespace.get(name) is not expected:
                raise ValueError('source:python-reference-replaced')
        for relative in self._bindings:
            if not relative.endswith('.py'):
                continue
            local = relative.removeprefix(_ARCHIVE).removesuffix('.py')
            name = local.replace('/', '.').removesuffix('.__init__')
            module = sys.modules.get(name)
            if module is None:
                continue
            record = contexts.get(name)
            if record is None:
                raise ValueError('source:unbound-python-context')
            if (module is not record[0]
                    or Path(module.__file__).resolve() != record[1]):
                raise ValueError('source:python-module-replaced')
        self._check_references()

    def protect(self, function):
        @wraps(function)
        def receiving_call(root=self.root, *, session=None):
            self.check(root)
            try:
                with self.imports():
                    return function(root, session=session)
            finally:
                self.check(root)
        return receiving_call
