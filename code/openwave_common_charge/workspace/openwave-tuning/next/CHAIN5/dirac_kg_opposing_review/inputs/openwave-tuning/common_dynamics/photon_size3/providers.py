"""Immutable ancestry for the size-three member of the quadratic source."""
from pathlib import Path
import hashlib
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = HERE.parent/'photon_audit/MANIFEST.json'
PARENT_SHA = 'ed6d2ccd6cba8021b148d640626118c38e390e443f8b2db50dbf94bcc6ee6fde'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    assert sha(PARENT) == PARENT_SHA
    data = json.loads(PARENT.read_text())
    result = {str(PARENT.relative_to(ROOT)):PARENT_SHA}
    for group in ('files', 'source_pins'):
        for name, digest in data[group].items():
            assert sha(ROOT/name) == digest, name
            assert name not in result or result[name] == digest
            result[name] = digest
    return result


def native():
    sys.path.insert(0, str(HERE.parent/'vector_field'))
    import native_field
    return native_field
