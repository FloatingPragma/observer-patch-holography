"""Execute two unchanged custody suites on current science and retained release metadata."""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import tempfile

import pytest

from . import model
from model.retained_identity_execution.runner import (RetainedExecutionError,
    case_failure, validate_imported_sources, validate_reports)

NAME = 'oph-current-science-retained-release'
TEST_FILES = frozenset({
    'model/completion/test_field_current_records_custody.py',
    'model/completion/test_field_current_records_navigation.py',
})
RETAINED_SHA = 'a68eefc1d92a87bd30ee2fbc11d2a38d41bb9062b69e83877f91accbf02f971c'


class Execution:
    def __init__(self, nodeids, root, log_directory):
        self.root = Path(root).resolve()
        self.approved = model.authority(root)
        self.nodeids = tuple(nodeids)
        allowed = set(self.approved['original_test_nodeids'])
        model.require(self.nodeids and len(set(self.nodeids)) == len(self.nodeids)
                      and set(self.nodeids) <= allowed, 'unreviewed current-custody test identity')
        self.log_directory = Path(log_directory)
        self.log_directory.mkdir(parents=True, exist_ok=True)
        self.reports = self.error = None
        self.source_snapshot_sha256 = None

    def run(self):
        if self.error is not None:
            raise RetainedExecutionError(self.error)
        if self.reports is not None:
            return self.reports
        try:
            self._run()
        except Exception as error:
            self.error = str(error)
            raise RetainedExecutionError(self.error) from error
        return self.reports

    def _run(self):
        release = model.validate_actual_release(self.root)
        with model.scientific_view(self.root) as (view, snapshot):
            with tempfile.TemporaryDirectory(prefix='oph-current-release-reporter-') as directory:
                observer = Path(directory)/'retained_case_reporter.py'
                original = self.root/'model/retained_identity_execution/reporter.py'
                shutil.copyfile(original, observer)
                model.require(model.sha(observer) == model.sha(original), 'current-case observer copy changed')
                journal = self.log_directory/'original-case-reports.json'
                journal.unlink(missing_ok=True)
                argv = ['-q', '-o', 'cache_dir='+str(Path(directory)/'pytest-cache'), *self.nodeids]
                bootstrap = ('import sys;sys.path[:0]='+repr([str(view),directory])+';'
                             'import pytest;from retained_case_reporter import Reporter;'
                             'raise SystemExit(pytest.main('+repr(argv)+',plugins=[Reporter('
                             +repr(str(journal))+','+repr(str(view))+')]))')
                environment = dict(os.environ)
                environment.pop('PYTHONPATH',None)
                environment.pop('PYTEST_ADDOPTS',None)
                environment['PYTHONDONTWRITEBYTECODE']='1'
                result = subprocess.run([sys.executable,'-I','-B','-W','error','-c',bootstrap],
                                        cwd=view, env=environment, capture_output=True, text=True)
                (self.log_directory/'original-pytest.log').write_text(result.stdout+result.stderr)
                model.require(journal.is_file(), 'current-custody pytest produced no case journal')
                value = json.loads(journal.read_text())
                model.require(value.get('exitstatus') == result.returncode, 'current-custody process/journal disagreement')
                reports = validate_reports(self.nodeids, value)
                imported = validate_imported_sources(value, view=view, root=self.root)
        model.require(model.validate_actual_release(self.root) == release,
                      'actual publication changed during original current-custody tests')
        record = {'scope': 'current scientific source with exactly one retained documentary manifest',
                  'selected_nodeids': list(self.nodeids), 'actual_original_tests_executed': True,
                  'outcomes_reused_from_receipts': False,
                  'current_publication_validation': release, 'current_source_view': snapshot,
                  'observed_imported_python_files_unchanged': imported,
                  'case_journal_sha256': model.sha(journal)}
        (self.log_directory/'execution.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
        self.source_snapshot_sha256 = snapshot['current_source_snapshot_sha256']
        self.reports = reports


class CurrentReleaseItem(pytest.Item):
    def __init__(self, *, original, execution, **kwargs):
        super().__init__(**kwargs)
        self.original_location = original.location
        self.execution = execution
        self.user_properties = list(getattr(original,'user_properties',[])) + [
            ('execution_scope','current_science_retained_documentary_manifest'),
            ('retained_manifest_sha256',RETAINED_SHA),
            ('current_manifest_sha256',execution.approved['current_manifest_sha256']),
            ('original_test_sha256',execution.approved['original_test_sha256'][self.nodeid.split('::')[0]]),
            ('current_physics_replaced_by_historical_source',False),
            ('source_view_and_case_journal',str(execution.log_directory)),
        ]

    def runtest(self):
        reports = self.execution.run()
        self.user_properties.append(('current_source_snapshot_sha256',self.execution.source_snapshot_sha256))
        failure = case_failure(reports[self.nodeid])
        if failure:
            raise RetainedExecutionError(failure)

    def reportinfo(self):
        return self.path, self.original_location[1], self.original_location[2]

    def repr_failure(self, excinfo, style=None):
        if isinstance(excinfo.value, RetainedExecutionError):
            return ('Current scientific source / retained documentary manifest execution\n'
                    +str(excinfo.value)+'\nFull journal: '+str(self.execution.log_directory))
        return super().repr_failure(excinfo, style=style)


class Plugin:
    def __init__(self):
        self.execution = None

    def pytest_report_header(self, config):
        return ('Documentary successor: two enumerated custody suites retain all current science; '
                'only their sealed release manifest is supplied from its archive. Actual publication checked separately.')

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_collection_modifyitems(self, config, items):
        yield
        root = Path(config.rootpath).resolve()
        if root != model.ROOT.resolve():
            return
        selected = [item for item in items if item.nodeid.split('::')[0] in TEST_FILES]
        if not selected or model.sha(root/model.MANIFEST) == RETAINED_SHA:
            return
        approved = model.authority(root)
        model.require(set(approved['original_test_sha256']) == TEST_FILES, 'current-custody test file scope changed')
        before = [item.nodeid for item in items]
        directory = tempfile.mkdtemp(prefix='execution-',dir=config.cache.makedir('documentary_release'))
        self.execution = Execution([item.nodeid for item in selected],root,directory)
        ids = {id(item) for item in selected}
        items[:] = [CurrentReleaseItem.from_parent(item.parent,name=item.name,nodeid=item.nodeid,
                    original=item,execution=self.execution) if id(item) in ids else item for item in items]
        model.require([item.nodeid for item in items] == before, 'documentary routing changed cases, order or count')
        reporter = config.pluginmanager.getplugin('terminalreporter')
        if reporter is not None:
            reporter.write_line(f'{len(selected)} original current-custody cases: current scientific files, '
                                'one archived manifest, and separate strict actual publication validation.')

    def pytest_terminal_summary(self, terminalreporter):
        if self.execution is not None:
            status = 'completed' if self.execution.reports is not None else 'not completed'
            terminalreporter.write_sep('-',f'Current-source documentary execution {status}: {len(self.execution.nodeids)} cases')
            terminalreporter.write_line('No scientific source frozen to history. Journal: '+str(self.execution.log_directory))


def install(config):
    if not config.pluginmanager.hasplugin(NAME):
        config.pluginmanager.register(Plugin(),NAME)
