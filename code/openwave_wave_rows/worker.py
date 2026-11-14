"""Isolated execution controller; all scientific modules keep their exact bytes."""
from pathlib import Path
import contextlib
import hashlib
import importlib.util
import importlib.machinery
import io
import json
import os
import runpy
import sys
import sysconfig

CONTROL = Path(__file__).resolve().parent
WORKSPACE = CONTROL / 'workspace'
MODE = sys.argv[1]
CONFIG = json.loads((CONTROL / 'runtime.json').read_text())
MANIFEST = json.loads((CONTROL / 'SNAPSHOT_MANIFEST.json').read_text())
STDLIB = Path(sysconfig.get_path('stdlib')).resolve()
SITES = [Path(p).resolve() for p in CONFIG['site_directories']]
RUNTIME_DISTRIBUTIONS = {name: name for name in ('numpy', 'scipy', 'sympy', 'mpmath', 'networkx', 'jsonschema', 'attrs', 'referencing')}
RUNTIME_DISTRIBUTIONS.update({'attr': 'attrs', 'yaml': 'PyYAML',
    'jsonschema_specifications': 'jsonschema-specifications', 'rpds': 'rpds-py',
    'typing_extensions': 'typing-extensions'})
STDLIB_ARCHIVES = {Path(p).resolve() for p in sys.path if p.endswith('.zip')}
for directory in SITES:
    if str(directory) not in sys.path:
        sys.path.append(str(directory))
sys.path[:0] = [str(CONTROL), str(WORKSPACE / 'reverse-engineering-reality')]
sys.dont_write_bytecode = True
READS, WRITES, BLOCKED = set(), set(), []
OPTIONAL_IMPORTS_BLOCKED = set()
PLAN = json.loads((CONTROL / 'PLANS.json').read_text())[MODE]
OUTPUTS = PLAN['outputs']


def is_runtime(path):
    return path.is_relative_to(STDLIB) or path in STDLIB_ARCHIVES or any(path.is_relative_to(p) for p in SITES)


def admitted_package_path(path):
    return any(path.is_relative_to(site / package) or path == site / (package + '.py')
               for site in SITES for package in RUNTIME_DISTRIBUTIONS)


class NumericRuntimeOnly:
    """Keep optional installed packages out of numerical-library imports."""
    def find_spec(self, fullname, path=None, target=None):
        spec = importlib.machinery.PathFinder.find_spec(fullname, path)
        if spec is not None and spec.origin and spec.origin not in ('built-in', 'frozen'):
            origin = Path(spec.origin).resolve()
            if any(origin.is_relative_to(site) for site in SITES):
                if not admitted_package_path(origin):
                    OPTIONAL_IMPORTS_BLOCKED.add(fullname)
                    raise ModuleNotFoundError('Optional external runtime excluded: ' + fullname, name=fullname)
        return None


sys.meta_path.insert(0, NumericRuntimeOnly())


def guard(event, args):
    if event == 'open' and isinstance(args[0], (str, bytes)):
        path = Path(os.fsdecode(args[0])).resolve()
        writing = isinstance(args[1], str) and any(c in args[1] for c in 'wax+')
        writing = writing or isinstance(args[2], int) and bool(args[2] & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC))
        if writing:
            permitted = path == CONTROL / (MODE + '.json') or path in [WORKSPACE / p for p in OUTPUTS]
            if not permitted:
                BLOCKED.append('write:' + str(path))
                raise PermissionError('write outside declared replay outputs: ' + str(path))
            WRITES.add(str(path.relative_to(CONTROL)))
        elif path.is_relative_to(WORKSPACE):
            relative = str(path.relative_to(WORKSPACE))
            aliases = CONFIG.get('runtime_aliases', {})
            if relative not in MANIFEST['files_sha256'] and relative not in aliases:
                raise PermissionError('unmanifested project read: ' + relative)
            READS.add(relative)
        elif path.is_relative_to(CONTROL) and path.name in ('worker.py', 'factors.py'):
            pass
        elif is_runtime(path) or str(path) in ('/dev/null', '/dev/urandom'):
            pass
        else:
            BLOCKED.append('read:' + str(path))
            raise PermissionError('external project/data read denied: ' + str(path))
    if event in ('subprocess.Popen', 'os.system', 'socket.connect', 'socket.bind'):
        raise PermissionError('scientific replay cannot launch external processes')


sys.addaudithook(guard)
try:
    Path(CONFIG['forbidden_probe']).read_bytes()
except PermissionError:
    pass
else:
    raise AssertionError('external project read probe succeeded')


def module_receipt():
    local, external = [], {}
    for name, module in list(sys.modules.items()):
        value = getattr(module, '__file__', None)
        if not value:
            continue
        path = Path(value).resolve()
        if path.is_relative_to(WORKSPACE):
            relative = str(path.relative_to(WORKSPACE))
            if relative not in MANIFEST['files_sha256'] and relative not in CONFIG.get('runtime_aliases', {}):
                raise AssertionError('unmanifested imported module: ' + relative)
            local.append(relative)
        elif path.is_relative_to(CONTROL):
            if path.name not in ('worker.py', 'factors.py'):
                raise AssertionError('unexpected controller import')
        elif any(path.is_relative_to(site) for site in SITES):
            if not admitted_package_path(path):
                raise AssertionError('undeclared external Python dependency: ' + name)
        elif not is_runtime(path):
            raise AssertionError('external project module: ' + name + ' at ' + str(path))
    for name, distribution in RUNTIME_DISTRIBUTIONS.items():
        version = getattr(sys.modules.get(name), '__version__', None)
        if version is not None and version != CONFIG['dependency_versions'][distribution]:
            raise AssertionError('unexpected numerical runtime version: ' + name)
        external[distribution] = CONFIG['dependency_versions'][distribution]
    return {'site_startup_disabled': bool(sys.flags.no_site),
            'environment_pythonpath_ignored': bool(sys.flags.isolated),
            'external_project_read_probe_rejected': 'read:' + CONFIG['forbidden_probe'] in BLOCKED,
            'all_project_modules_inside_snapshot': True,
            'snapshot_module_files': sorted(set(local)),
            'snapshot_read_files': sorted(READS), 'declared_output_writes': sorted(WRITES),
            'blocked_read_write_attempts': [value.replace(CONFIG['forbidden_probe'], '<outside-workspace-probe>')
                .replace(sys.prefix, '<python-prefix>').replace(str(Path.home()), '<home>') for value in BLOCKED],
            'python_version': sys.version,
            'optional_external_imports_excluded': sorted(OPTIONAL_IMPORTS_BLOCKED),
            'numerical_library_versions': external}


captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    script = WORKSPACE / PLAN['script']
    sys.path.insert(0, str(script.parent))
    runpy.run_path(str(script), run_name='__main__')
runtime = module_receipt()
assert runtime['external_project_read_probe_rejected']
outputs = {relative: hashlib.sha256((WORKSPACE / relative).read_bytes()).hexdigest() for relative in OUTPUTS}
result = {'schema': 'oph.wave_rows.isolated_execution.v1', 'mode': MODE,
          'regenerated_outputs': outputs, 'runtime': runtime,
          'execution_log': captured.getvalue()}
(CONTROL / (MODE + '.json')).write_text(json.dumps(result, indent=2) + '\n')
print(MODE + ': isolated scientific replay passed', flush=True)
