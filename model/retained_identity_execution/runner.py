"""Execute selected original cases in a fresh authenticated source checkout."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from . import live


class RetainedExecutionError(RuntimeError):
    pass


def validate_reports(nodeids, value):
    """A process exit or a cached pass cannot replace any selected case."""
    expected = set(nodeids)
    if len(expected) != len(nodeids) or not expected:
        raise RetainedExecutionError('retained case selection is empty or duplicated')
    by_case = {name: {} for name in nodeids}
    for report in value.get('reports', []):
        name, when = report.get('nodeid'), report.get('when')
        if name not in expected or when not in {'setup', 'call', 'teardown'}:
            raise RetainedExecutionError('unexpected original case or report phase')
        if when in by_case[name] or report.get('outcome') not in {'passed', 'failed', 'skipped'}:
            raise RetainedExecutionError('duplicate or malformed original case outcome')
        by_case[name][when] = report
    missing = [name for name, phases in by_case.items() if not phases]
    if missing:
        raise RetainedExecutionError('original cases were not executed: ' + ', '.join(missing))
    any_failure = False
    for name, phases in by_case.items():
        if set(phases) != {'setup', 'call', 'teardown'}:
            if not any(r['outcome'] != 'passed' for r in phases.values()):
                raise RetainedExecutionError('original case has incomplete successful phases: ' + name)
        any_failure |= any(r['outcome'] != 'passed' or r.get('wasxfail') for r in phases.values())
    if value.get('exitstatus') != 0 and not any_failure:
        raise RetainedExecutionError('original pytest process failed outside its case outcomes')
    return by_case


def case_failure(phases):
    """Retain real failures; an unexpected skip is never turned into a pass."""
    failures = [r for r in phases.values() if r['outcome'] != 'passed' or r.get('wasxfail')]
    if failures:
        return '\n'.join(r.get('longrepr') or ('Original case outcome: ' + r['outcome'])
                         for r in failures)
    if set(phases) != {'setup', 'call', 'teardown'}:
        return 'Original case did not complete every pytest phase.'
    return None


def validate_imported_sources(value, *, view, root):
    """An observed executable change outside the frozen input map is fatal too."""
    view, root = Path(view).resolve(), Path(root).resolve()
    imported = value.get('imported_repository_python_sha256')
    if not isinstance(imported, dict) or not imported:
        raise RetainedExecutionError('original execution did not report its imported source files')
    for name, expected in imported.items():
        before, after = (view / name).resolve(), (root / name).resolve()
        if (not before.is_relative_to(view) or not after.is_relative_to(root)
                or before.suffix != '.py' or not before.is_file() or not after.is_file()
                or live.sha(before) != expected or live.sha(after) != expected):
            raise RetainedExecutionError('unreviewed imported executable difference: ' + name)
    return len(imported)


class Execution:
    def __init__(self, nodeids, *, root=live.ROOT, log_directory=None):
        self.root = Path(root).resolve()
        self.approved = live.authority(self.root)
        self.nodeids = tuple(nodeids)
        allowed = {n for rows in self.approved['nodeids_by_file'].values() for n in rows}
        if (not self.nodeids or len(set(self.nodeids)) != len(self.nodeids)
                or not set(self.nodeids) <= allowed):
            raise RetainedExecutionError('unreviewed, missing, or duplicated retained test identity')
        self.log_directory = Path(log_directory) if log_directory else Path(tempfile.mkdtemp(prefix='oph-retained-case-log-'))
        self.log_directory.mkdir(parents=True, exist_ok=True)
        self.reports = None
        self.error = None
        self.audit = None

    def _command(self, argv, *, cwd, log_name):
        environment = dict(os.environ)
        environment.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
                           OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
        environment.pop('PYTHONPATH', None)
        environment.pop('PYTEST_ADDOPTS', None)
        result = subprocess.run(argv, cwd=cwd, env=environment, capture_output=True, text=True)
        (self.log_directory / log_name).write_text(result.stdout + result.stderr)
        return result

    def run(self):
        if self.reports is not None:
            return self.reports
        if self.error is not None:
            raise RetainedExecutionError(self.error)
        try:
            self._run()
        except Exception as error:
            self.error = str(error)
            raise RetainedExecutionError(self.error) from error
        return self.reports

    def _run(self):
        with tempfile.TemporaryDirectory(prefix='oph-retained-identity-') as folder:
            base = Path(folder).resolve()
            view = base / 'source'
            result = self._command(['git', 'clone', '--quiet', '--shared', '--no-checkout',
                                    str(self.root), str(view)], cwd=self.root, log_name='clone.log')
            if result.returncode:
                raise RetainedExecutionError('cannot create the private historical source view: ' + result.stderr)
            result = self._command(['git', 'checkout', '--quiet', '--detach', self.approved['revision']],
                                   cwd=view, log_name='checkout.log')
            if result.returncode:
                raise RetainedExecutionError('full retained Git history is required: ' + result.stderr)
            revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=view, text=True).strip()
            if revision != self.approved['revision']:
                raise RetainedExecutionError('historical checkout has the wrong revision')
            live.validate_view(view, self.approved, root=self.root)
            self.audit = live.validate_current(view, root=self.root)
            if any(n.startswith('model/particle_atlas_horizontal/') for n in self.nodeids):
                bootstrap = ('import runpy,sys;sys.path.insert(0,' + repr(str(view)) + ');'
                             'sys.argv=["model.particle_atlas_horizontal.build","--write"];'
                             'runpy.run_module("model.particle_atlas_horizontal.build",run_name="__main__")')
                built = self._command([sys.executable, '-I', '-B', '-W', 'error', '-c', bootstrap],
                                      cwd=view, log_name='historical-atlas-build.log')
                if built.returncode:
                    raise RetainedExecutionError('unchanged historical atlas producer failed: ' + built.stdout + built.stderr)
            reporter = base / 'retained_case_reporter.py'
            reporter_source = Path(__file__).with_name('reporter.py')
            shutil.copyfile(reporter_source, reporter)
            if live.sha(reporter) != live.sha(reporter_source):
                raise RetainedExecutionError('outcome observer source changed during copying')
            journal = self.log_directory / 'original-case-reports.json'
            journal.unlink(missing_ok=True)
            arguments = ['-q', *self.nodeids]
            bootstrap = ('import sys;sys.path[:0]=' + repr([str(view), str(base)]) + ';'
                         'import pytest;from retained_case_reporter import Reporter;'
                         'raise SystemExit(pytest.main(' + repr(arguments) + ',plugins=[Reporter(' + repr(str(journal)) + ',' + repr(str(view)) + ')]))')
            result = self._command([sys.executable, '-I', '-B', '-W', 'error', '-c', bootstrap],
                                   cwd=view, log_name='original-pytest.log')
            if not journal.is_file():
                raise RetainedExecutionError('original pytest did not produce its complete case journal:\n'
                                             + result.stdout + result.stderr)
            value = json.loads(journal.read_text())
            if value.get('exitstatus') != result.returncode:
                raise RetainedExecutionError('original pytest exit status differs from its journal')
            checked_reports = validate_reports(self.nodeids, value)
            imported_count = validate_imported_sources(value, view=view, root=self.root)
            live.validate_view(view, self.approved, root=self.root)
            self.audit = live.validate_current(view, root=self.root)
            record = {'revision': revision, 'selected_nodeids': list(self.nodeids),
                      'actual_original_tests_executed': True, 'outcomes_reused_from_receipts': False,
                      'test_and_source_bytes_changed': False, 'live_dependency_validation': self.audit,
                      'current_dependencies_rechecked_after_actual_original_execution': True,
                      'observed_imported_python_files_unchanged': imported_count,
                      'case_journal_sha256': live.sha(journal),
                      'case_properties': {'execution_scope': 'retained_scientific_source',
                                          'source_commit': revision,
                                          'current_physics_verified_by_historical_pass': False}}
            (self.log_directory / 'execution.json').write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
            self.reports = checked_reports

