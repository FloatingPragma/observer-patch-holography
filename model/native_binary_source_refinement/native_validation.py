"""Check source roots in receiver-owned Lean sessions with current import custody."""
from __future__ import annotations

import hashlib
from contextlib import contextmanager
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

PREFIX = 'EventAlgebra.NativeBinarySource'
ROOTS = (PREFIX+'Endpoints', 'EventAlgebra.NativeL2CountMoment',
         PREFIX+'Witnesses', PREFIX+'Replay', PREFIX+'GrammarCheck')
_LIBRARY_SNAPSHOTS = {}
_LIBRARY_MODULES = {}
_SESSIONS = {}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _resolve_import(module, roots):
    relative = Path(*module.split('.')).with_suffix('.olean')
    return next(((root/relative).resolve() for root in roots
                 if (root/relative).is_file()), None)


def _check_library_members(artifacts, resolutions, roots):
    for row in resolutions:
        current = _resolve_import(row['module'], roots)
        if current is None or str(current) != row['artifact_path']:
            raise ValueError('kernel:library-import-resolution-changed')
    for row in artifacts:
        path = Path(row['path'])
        if not row.get('exists', True):
            if path.exists():
                raise ValueError('kernel:library-artifact-appeared')
            continue
        if not path.is_file():
            raise ValueError('kernel:library-artifact-missing')
        if digest(path) != row['sha256']:
            raise ValueError('kernel:library-artifact-changed')


class NativeValidationSession:
    """A session never accepts receipt booleans, old logs or prebuilt local objects."""

    def __init__(self, root, evidence=None):
        self.root = Path(root).resolve()
        self.lean = self.root/'Lean'
        parent = Path(evidence) if evidence else self.root/'temp/native-binary-source'
        parent.mkdir(parents=True, exist_ok=True)
        self.output = Path(tempfile.mkdtemp(prefix='session-', dir=parent))
        self.units, self.events = {}, []
        lake = shutil.which('lake')
        if lake is None:
            raise ValueError('kernel:lean-required')
        command = [lake, 'env', sys.executable, '-B', '-c',
                   'import json,os,subprocess; print(json.dumps([os.environ.get("LEAN_PATH",""),subprocess.check_output(["lean","--print-prefix"],text=True).strip()]))']
        result = subprocess.run(command, cwd=self.lean, capture_output=True, text=True,
                                check=True, encoding='utf-8')
        paths, prefix = json.loads(result.stdout)
        suffix = '.exe' if os.name == 'nt' else ''
        self.compiler = Path(prefix)/'bin'/('lean'+suffix)
        project_output = (self.lean/'.lake/build/lib/lean').resolve()
        self.libraries = [Path(p).resolve() for p in paths.split(os.pathsep)
                          if p and Path(p).resolve() != project_output]
        self.libraries = list(dict.fromkeys([*self.libraries,
                                            (Path(prefix)/'lib/lean').resolve()]))

    def settings(self):
        files = [self.compiler, self.lean/'lean-toolchain',
                 self.lean/'lake-manifest.json', self.lean/'lakefile.lean']
        manifest = json.loads((self.lean/'lake-manifest.json').read_bytes())
        for package in manifest['packages']:
            base = self.lean/'.lake/packages'/package['name']
            for name in ('lean-toolchain', 'lake-manifest.json', 'lakefile.lean', 'lakefile.toml'):
                if (base/name).is_file():
                    files.append(base/name)
        identity = {'compiler': str(self.compiler),
                    'files': [(str(p.resolve()), digest(p)) for p in files],
                    'libraries': list(map(str, self.libraries))}
        return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()

    def library_snapshot(self, imports, runtime):
        roots = self.libraries
        root_key = tuple(map(str, roots))
        seeds = tuple(sorted(imports))
        key = (seeds, runtime, root_key)
        previous = _LIBRARY_SNAPSHOTS.get(key)
        if previous:
            checked = getattr(self, '_checked_libraries', None)
            if checked is not None and key in checked:
                return previous['identity']
            _check_library_members(previous['artifacts'], previous['resolutions'], roots)
            if checked is not None:
                checked.add(key)
            return previous['identity']
        pending, visited, artifacts, resolutions = list(seeds), set(), [], []
        while pending:
            module = pending.pop()
            if module in visited:
                continue
            visited.add(module)
            module_key = (runtime, root_key, module)
            cached = _LIBRARY_MODULES.get(module_key)
            if cached:
                _check_library_members(cached['artifacts'], [cached['resolution']], roots)
                artifacts.extend(cached['artifacts'])
                resolutions.append(cached['resolution'])
                pending.extend(cached['imports'])
                continue
            artifact = _resolve_import(module, roots)
            if artifact is None:
                raise ValueError(f'kernel:missing-import-artifact:{module}')
            metadata = artifact.with_suffix('.ilean')
            if not metadata.is_file():
                raise ValueError(f'kernel:missing-import-closure:{module}')
            direct = [row[0] for row in json.loads(metadata.read_bytes())['directImports']]
            members = [artifact, metadata, Path(str(artifact)+'.private'), Path(str(artifact)+'.server')]
            rows = [{'path': str(p), 'exists': p.is_file(),
                     'sha256': digest(p) if p.is_file() else None} for p in members]
            resolution = {'module': module, 'artifact_path': str(artifact)}
            _LIBRARY_MODULES[module_key] = {'artifacts': rows, 'imports': direct,
                                           'resolution': resolution}
            artifacts.extend(rows)
            resolutions.append(resolution)
            pending.extend(direct)
        snapshot = {'seeds': seeds, 'roots': root_key,
                    'artifacts': sorted(artifacts, key=lambda row: row['path']),
                    'resolutions': sorted(resolutions, key=lambda row: row['module'])}
        identity = hashlib.sha256(json.dumps(snapshot, sort_keys=True).encode()).hexdigest()
        (self.output/(identity+'-libraries.json')).write_text(json.dumps(snapshot, indent=2)+'\n', encoding='utf-8')
        _LIBRARY_SNAPSHOTS[key] = dict(snapshot, identity=identity)
        checked = getattr(self, '_checked_libraries', None)
        if checked is not None:
            checked.add(key)
        return identity

    def unit(self, name, stack=()):
        if name in stack:
            raise ValueError('kernel:local-import-cycle')
        source = self.lean/Path(*name.split('.')).with_suffix('.lean')
        text = source.read_text(encoding='utf-8')
        imports = [name for line in text.splitlines() if line.startswith('import ')
                   for name in line[7:].split('--', 1)[0].split()]
        local = [m for m in imports if (self.lean/Path(*m.split('.')).with_suffix('.lean')).is_file()]
        vendor = [m for m in imports if m not in local]
        parents = [self.unit(m, (*stack, name)) for m in local]
        runtime = self.settings()
        identity = {'module': name, 'source': digest(source), 'runtime': runtime,
                    'libraries': self.library_snapshot(vendor, runtime),
                    'parents': [(p['identity'], p['object_sha256']) for p in parents]}
        key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        previous = self.units.get(key)
        if previous:
            for row in previous['artifacts']:
                path = Path(row['path'])
                if not row.get('exists', True):
                    if path.exists():
                        raise ValueError('kernel:private-artifact-changed')
                elif not path.is_file() or digest(path) != row['sha256']:
                    raise ValueError('kernel:private-artifact-changed')
            self.events.append({'module': name, 'event': 'validated-reuse', 'identity': key})
            return previous
        folder = self.output/key
        target = folder/Path(*name.split('.')).with_suffix('.olean')
        target.parent.mkdir(parents=True, exist_ok=True)
        parent_roots = list(dict.fromkeys(p for parent in parents for p in parent['import_roots']))
        # Lean resolves the shared EventAlgebra prefix from one import root.
        # Put the complete checked local closure in that receiver-owned root.
        for parent_root in parent_roots:
            for parent_object in sorted(parent_root.rglob('*.olean*')):
                copy = folder/parent_object.relative_to(parent_root)
                copy.parent.mkdir(parents=True, exist_ok=True)
                if copy.is_file() and copy.read_bytes() != parent_object.read_bytes():
                    raise ValueError('kernel:conflicting-local-closure')
                shutil.copyfile(parent_object, copy)
        env = dict(os.environ, LEAN_PATH=os.pathsep.join(map(str, [folder, *self.libraries])))
        command = [str(self.compiler), '-R', str(self.lean), '-o', str(target), str(source)]
        tick = time.perf_counter()
        completed = subprocess.run(command, cwd=self.lean, env=env, capture_output=True)
        stdout, stderr = folder/'stdout.log', folder/'stderr.log'
        stdout.write_bytes(completed.stdout)
        stderr.write_bytes(completed.stderr)
        self.events.append({'module': name, 'event': 'native-check', 'identity': key,
                            'command': command, 'exit_code': completed.returncode,
                            'seconds': time.perf_counter()-tick,
                            'stdout_sha256': digest(stdout), 'stderr_sha256': digest(stderr)})
        self.persist()
        if digest(source) != identity['source']:
            raise ValueError('kernel:source-changed-during-check')
        if completed.returncode:
            raise ValueError(f'kernel:current-source-rejected:{name}')
        output = completed.stdout.decode('utf-8')
        for values in re.findall(r'depends on axioms: \[([^]]*)\]', output, re.S):
            if {v.strip() for v in values.split(',') if v.strip()} - {'propext', 'Classical.choice', 'Quot.sound'}:
                raise ValueError('kernel:unapproved-proof-dependency')
        for theorem in re.findall(r'^(?:theorem|lemma)\s+([A-Za-z0-9_]+)', text, re.M):
            if f'.{theorem}\'' not in output:
                raise ValueError(f'kernel:missing-axiom-audit:{name}.{theorem}')
        artifacts = [{'path': str(p), 'sha256': digest(p)} for p in
                     [stdout, stderr, *sorted(folder.rglob('*.olean*'))]]
        for obj in sorted(folder.rglob('*.olean')):
            for suffix in ('.private', '.server'):
                companion = Path(str(obj)+suffix)
                if not companion.is_file():
                    artifacts.append({'path': str(companion), 'exists': False, 'sha256': None})
        record = {'identity': key, 'object_sha256': digest(target), 'object': target,
                  'artifacts': artifacts, 'import_roots': [folder]}
        self.units[key] = record
        return record

    def persist(self):
        (self.output/'events.json').write_text(json.dumps(self.events, indent=2)+'\n', encoding='utf-8')

    @contextmanager
    def validation_epoch(self):
        self._checked_libraries = set()
        before = self.settings()
        try:
            yield
            if self.settings() != before:
                raise ValueError('kernel:runtime-changed-during-check')
        finally:
            self._checked_libraries = None

    def validate(self):
        with self.validation_epoch():
            for name in ROOTS:
                self.unit(name)
        self.persist()

    def execute(self, module, args):
        if module not in (PREFIX+'Replay', PREFIX+'GrammarCheck'):
            raise ValueError('kernel:unknown-runtime-service')
        with self.validation_epoch():
            record = self.unit(module)
        namespace = 'Replay' if module.endswith('Replay') else 'GrammarCheck'
        wrapper = self.output/(namespace+'Main.lean')
        wrapper.write_text(f'import {module}\ndef main := EventAlgebra.NativeBinarySource.{namespace}.main\n', encoding='utf-8')
        inputs = [{'path': str(Path(p).resolve()), 'sha256': digest(p)} for p in [self.compiler, wrapper, *args]]
        env = dict(os.environ, LEAN_PATH=os.pathsep.join(map(str, [*record['import_roots'], *self.libraries])))
        command = [str(self.compiler), '--run', str(wrapper), *map(str, args)]
        completed = subprocess.run(command, cwd=self.lean, env=env, capture_output=True)
        output = self.output/'runtime'/str(len(self.events))
        output.mkdir(parents=True)
        (output/'stdout.log').write_bytes(completed.stdout)
        (output/'stderr.log').write_bytes(completed.stderr)
        self.events.append({'module': module, 'event': 'execute', 'command': command,
                            'inputs': inputs, 'proof_identity': record['identity'],
                            'stdout_path': str(output/'stdout.log'),
                            'stderr_path': str(output/'stderr.log'),
                            'exit_code': completed.returncode,
                            'stdout_sha256': hashlib.sha256(completed.stdout).hexdigest(),
                            'stderr_sha256': hashlib.sha256(completed.stderr).hexdigest()})
        self.persist()
        if any(digest(row['path']) != row['sha256'] for row in inputs):
            raise ValueError('kernel:runtime-input-changed')
        if completed.returncode:
            raise ValueError(f'kernel:execution-rejected:{module}')
        return completed.stdout


def current_session(root, evidence=None):
    key = (str(Path(root).resolve()), str(Path(evidence).resolve()) if evidence else None)
    if key not in _SESSIONS:
        _SESSIONS[key] = NativeValidationSession(root, evidence)
    return _SESSIONS[key]
