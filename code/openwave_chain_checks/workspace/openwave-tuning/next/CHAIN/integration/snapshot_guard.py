"""Reject mutable-project reads in the historical execution workspace."""
from pathlib import Path
import os
import sys
import json

ACTIVE=False
READS=set()
BLOCKED=[]
ROOT=None
LIVE=None


def install(directory):
    global ACTIVE,ROOT,LIVE
    here=Path(directory).resolve();ROOT=here.parents[3]
    metadata=json.loads((here/'HISTORICAL_SOURCE.json').read_text())
    LIVE=Path(metadata['original_live_root']).resolve()
    if ROOT==LIVE:raise RuntimeError('historical snapshot must be separate from live source')
    sys.dont_write_bytecode=True
    def audit(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes)):return
        path=Path(os.fsdecode(args[0])).resolve()
        if path.is_relative_to(LIVE) and not path.is_relative_to(ROOT):
            BLOCKED.append(str(path))
            raise PermissionError('mutable project access denied: '+str(path))
        if path.is_relative_to(ROOT):READS.add(str(path.relative_to(ROOT)))
    sys.addaudithook(audit)
    ACTIVE=True
    try:
        (LIVE/'openwave-push/issues/N3D_cubic_routing/n3dlib.py').read_bytes()
    except PermissionError:
        pass
    else:
        raise RuntimeError('mutable source read guard failed')


def receipt():
    if not ACTIVE:raise RuntimeError('execute the historical source with run.py')
    imported=[]
    for module in list(sys.modules.values()):
        value=getattr(module,'__file__',None)
        if value:
            path=Path(value).resolve()
            if path.is_relative_to(LIVE) and not path.is_relative_to(ROOT):
                raise RuntimeError('mutable project module imported: '+str(path))
            if path.is_relative_to(ROOT):imported.append(str(path.relative_to(ROOT)))
    return {'historical_workspace':True,'all_successful_project_reads_inside_snapshot':True,
            'mutable_source_read_probe_rejected':True,'blocked_mutable_reads':BLOCKED,
            'snapshot_files_read':len(READS),'snapshot_module_files':sorted(set(imported)),
            'python_version':sys.version,'numerical_library_versions':{name:getattr(sys.modules.get(name),'__version__',None) for name in ('numpy','scipy','sympy','mpmath')},
            'site_startup_disabled':bool(sys.flags.no_site),'environment_pythonpath_ignored':bool(sys.flags.isolated)}
