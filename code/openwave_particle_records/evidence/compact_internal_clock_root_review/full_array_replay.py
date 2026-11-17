"""Run unchanged source into a disposable output copy and compare every array."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'compact_internal_clock'
spec = importlib.util.spec_from_file_location('sealed_arc_source', SOURCE / 'run_serialized.py')
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
original = json.loads((SOURCE / 'RECEIPT.json').read_text())
start = time.perf_counter()
with tempfile.TemporaryDirectory(prefix='oph-arc-array-review-') as temporary:
    output = Path(temporary) / 'packet'
    shutil.copytree(SOURCE, output)
    # Only the output location changes. ROOT, PROVIDER, native code and all
    # scientific operations remain those of the sealed implementation.
    module.HERE = output
    scientific, seconds = module.experiment(record=True)
    assert scientific == original['scientific']
    hashes = {}
    for relative, expected in original['artifacts'].items():
        actual = hashlib.sha256((output / relative).read_bytes()).hexdigest()
        assert actual == expected, relative
        hashes[relative] = actual
print(json.dumps({'all_arrays_byte_identical': True, 'source_seconds': seconds,
                  'elapsed_seconds': time.perf_counter() - start,
                  'artifacts': hashes}, indent=2))
