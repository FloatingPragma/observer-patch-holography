"""Isolated standard-library worker; no mutable workspace access."""
from pathlib import Path
import hashlib
import json
import os
import runpy
import sys
import sysconfig

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'workspace'
MODE = sys.argv[1]
sys.dont_write_bytecode = True
snapshot = json.loads((HERE / 'SNAPSHOT_MANIFEST.json').read_text())
row = 'openwave-tuning/row_completion/ROW_01/'
entries = {
    'candidate': [row + 'run.py', row + 'vacuum_control.py'],
    'auditor_root': ['openwave-tuning/row_completion/reviews/ROW_01_ROOT/audit.py'],
    'auditor_dyn17': ['openwave-tuning/row_completion/reviews/ROW_01_DYN17/audit.py'],
    'codex': ['openwave-tuning/next/CHAIN3/row1_full_criterion_review/review.py'],
}
outputs = {
    'candidate': [row + n for n in ('GEOMETRY.json', 'EXAMPLES.json', 'RESULT.json', 'VACUUM.json')],
    'auditor_root': ['openwave-tuning/row_completion/reviews/ROW_01_ROOT/RESULT.json'],
    'auditor_dyn17': ['openwave-tuning/row_completion/reviews/ROW_01_DYN17/RESULT.json'],
    'codex': ['openwave-tuning/next/CHAIN3/row1_full_criterion_review/' + n for n in ('REVIEW.json', 'READOUTS.json')],
}
for name in outputs[MODE]:
    (ROOT / name).unlink()

stdlib = Path(sysconfig.get_path('stdlib')).resolve()
reads = set()


def guard(event, args):
    if event == 'open' and isinstance(args[0], (str, bytes, os.PathLike)):
        p = Path(os.fsdecode(args[0])).resolve()
        mode = args[1]
        flags = args[2]
        writing = (isinstance(mode, str) and any(x in mode for x in 'wax+')) or (
            isinstance(flags, int) and bool(flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT)))
        if writing:
            if p != HERE / 'REPLAY.json' and not p.is_relative_to(ROOT):
                raise PermissionError('Write outside isolated workspace')
        elif p.is_relative_to(ROOT):
            reads.add(str(p.relative_to(ROOT)))
        elif p != HERE / 'worker.py' and not p.is_relative_to(stdlib):
            raise PermissionError('Read outside pinned workspace and standard library: ' + str(p))
    if event in ('socket.connect', 'subprocess.Popen', 'os.system'):
        raise PermissionError('External execution and network are outside this replay')


sys.addaudithook(guard)
denied = False
try:
    (HERE.parent / 'outside_mutable_workspace.py').read_text()
except PermissionError:
    denied = True
assert denied
for entry in entries[MODE]:
    sys.path.insert(0, str((ROOT / entry).parent))
    runpy.run_path(str(ROOT / entry), run_name='__main__')
for name, expected in snapshot['files_sha256'].items():
    actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    if actual != expected:
        raise AssertionError('Regenerated output or pinned source differs: ' + name)
project_imports = {}
for name, module in tuple(sys.modules.items()):
    value = getattr(module, '__file__', None)
    if value and Path(value).resolve().is_relative_to(ROOT):
        p = Path(value).resolve()
        project_imports[name] = str(p.relative_to(ROOT))
        assert project_imports[name] in snapshot['files_sha256']
receipt = {'mode': MODE, 'all_source_and_regenerated_output_bytes_match': True,
           'outside_workspace_read_rejected': True, 'stdlib_only': True,
           'outputs_deleted_and_regenerated': outputs[MODE],
           'project_imports': project_imports, 'pinned_file_count': len(snapshot['files_sha256']),
           'project_files_read': sorted(reads)}
(HERE / 'REPLAY.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
