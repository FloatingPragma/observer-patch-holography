"""Preserve every selected node while executing its original frozen identity."""
from pathlib import Path
import pytest

from . import live
from .runner import Execution, RetainedExecutionError, case_failure

PLUGIN_NAME = 'oph-retained-scientific-identities'


class RetainedItem(pytest.Item):
    def __init__(self, *, original, execution, **kwargs):
        super().__init__(**kwargs)
        self.original_location = original.location
        self.execution = execution
        self.user_properties = list(getattr(original, 'user_properties', [])) + [
            ('execution_scope', 'retained_scientific_source'),
            ('source_commit', execution.approved['revision']),
            ('original_test_sha256', execution.approved['test_files'][self.nodeid.split('::')[0]]),
            ('current_physics_verified_by_historical_pass', False),
        ]

    def runtest(self):
        reports = self.execution.run()
        failure = case_failure(reports[self.nodeid])
        if failure:
            raise RetainedExecutionError(failure)

    def repr_failure(self, excinfo, style=None):
        if isinstance(excinfo.value, RetainedExecutionError):
            return ('Actual retained-source execution at ' + self.execution.approved['revision']
                    + '\n' + str(excinfo.value) + '\nFull case journal: '
                    + str(self.execution.log_directory))
        return super().repr_failure(excinfo, style=style)

    def reportinfo(self):
        return self.path, self.original_location[1], self.original_location[2]


class Plugin:
    def __init__(self):
        self.execution = None

    def pytest_report_header(self, config):
        return ('Frozen-identity policy: explicitly enumerated cases execute unchanged at '
                'bc63b3dd8b7a352337b2080d025cbd46934bc36f; unlisted tests execute live.')

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_collection_modifyitems(self, config, items):
        # The unchanged shard/filter hooks finish first. No node is selected,
        # dropped, or assigned to a different shard by this wrapper.
        yield
        root = Path(config.rootpath).resolve()
        if root != live.ROOT.resolve():
            return
        approved = live.authority(root)
        retained = [item for item in items if item.nodeid.split('::')[0] in approved['test_files']]
        if not retained:
            return
        before = [item.nodeid for item in items]
        logs = Path(config.cache.makedir('retained_identity'))
        import tempfile
        directory = tempfile.mkdtemp(prefix='execution-', dir=logs)
        self.execution = Execution([item.nodeid for item in retained], root=root, log_directory=directory)
        selected = {id(item) for item in retained}
        items[:] = [RetainedItem.from_parent(item.parent, name=item.name, nodeid=item.nodeid,
                                            original=item, execution=self.execution)
                    if id(item) in selected else item for item in items]
        if before != [item.nodeid for item in items]:
            raise pytest.UsageError('retained routing changed test identities, counts or order')
        reporter = config.pluginmanager.getplugin('terminalreporter')
        if reporter is not None:
            reporter.write_line(f'{len(retained)} frozen-identity cases retain their node IDs and original bodies; '
                                'actual bc63 subprocess execution and separate current dependency checks required.')

    def pytest_terminal_summary(self, terminalreporter):
        if self.execution is not None:
            status = 'completed' if self.execution.reports is not None else 'not completed'
            terminalreporter.write_sep('-', f'Retained execution {status}: {len(self.execution.nodeids)} original '
                                      f'cases at {self.execution.approved["revision"]}')
            terminalreporter.write_line('Historical outcomes do not verify current physics. '
                                       'Current source equivalence is checked separately; journal: '
                                       + str(self.execution.log_directory))


def install(config):
    if not config.pluginmanager.hasplugin(PLUGIN_NAME):
        config.pluginmanager.register(Plugin(), PLUGIN_NAME)


def pytest_configure(config):
    install(config)
