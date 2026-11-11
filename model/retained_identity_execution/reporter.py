"""Observe original pytest outcomes without altering collection or execution."""
import json
import hashlib
from pathlib import Path
import sys


class Reporter:
    def __init__(self, output, source_root=None):
        self.output = Path(output)
        self.source_root = Path(source_root).resolve() if source_root else None
        self.reports = []

    def pytest_runtest_logreport(self, report):
        self.reports.append({'nodeid': report.nodeid, 'when': report.when,
                             'outcome': report.outcome, 'duration': report.duration,
                             'longrepr': report.longreprtext if report.failed or report.skipped else '',
                             'wasxfail': str(getattr(report, 'wasxfail', ''))})

    def pytest_sessionfinish(self, session, exitstatus):
        imported = {}
        if self.source_root is not None:
            for module in tuple(sys.modules.values()):
                filename = getattr(module, '__file__', None)
                if not filename:
                    continue
                path = Path(filename).resolve()
                if path.suffix == '.py' and path.is_relative_to(self.source_root) and path.is_file():
                    imported[path.relative_to(self.source_root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.output.write_text(json.dumps({'exitstatus': int(exitstatus),
                                           'reports': self.reports,
                                           'imported_repository_python_sha256': imported}, indent=2, sort_keys=True) + '\n')
