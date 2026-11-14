"""Conservative imported-source and declared evidence closure for the walk."""
import hashlib
import json
from pathlib import Path
import sys

import common_walk as common
import compile_walk as compiled
from model.native_clock import configuration

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RER = ROOT/'reverse-engineering-reality'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    common.graph()
    summary = compiled.compiler.source_summary()
    receipt = json.loads((HERE/'COMMON_WALK_RECEIPT.json').read_text())
    for name, pin in receipt['source_sha256'].items():
        assert sha(ROOT/name) == pin, name
    assert sha(HERE/'COMMON_WALK_PROGRAM.json') == receipt['program_sha256']
    paths = {HERE/name for name in ('common_walk.py', 'compile_walk.py',
        'COMMON_WALK_PROGRAM.json', 'COMMON_WALK_RECEIPT.json', 'COMMON_SOURCE_API.md',
        'pin_source_closure.py')}
    for module in list(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.is_file() and (path.is_relative_to(RER) or path.is_relative_to(HERE)):
                paths.add(path)
    paths.update(RER/name for name in compiled.pointer.SOURCE_PATHS)
    paths.update(RER/name for name in configuration.SOURCE_PATHS)
    paths.update(RER/name for name in summary['source_sha256'])
    paths.add(common.spatial.slow_geometry.SPEC)
    pins = {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}
    return {'schema': 'oph.common.local12.source_custody.v1',
            'model_id': common.MODEL_ID,
            'checks': {'primary_receipt_and_program_pins_match': True,
                       'imported_repository_source_closure_pinned': True,
                       'declared_pointer_reset_and_control_evidence_pinned': True},
            'files_sha256': pins, 'file_count': len(pins),
            'scope': 'Conservative imported Python-source closure plus explicit native geometry, processor, pointer and reset evidence; unused imports convey no scientific transfer.'}


if __name__ == '__main__':
    result = run()
    (HERE/'COMMON_SOURCE_CUSTODY.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': result['checks'], 'file_count': result['file_count']}))
