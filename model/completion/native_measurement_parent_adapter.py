"""Live custody for the complete 66-word measurement receipt.

A consumer reads model/native_energy_measurement/results.json only when the
receipt's recorded source map equals the measurement producer's current one.
This adapter performs no measurement, spectral or clock calculation.
"""
import hashlib
import json
from pathlib import Path

from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
PARENT = 'model/native_energy_measurement/results.json'
SOURCE_FILES = ('model/completion/native_measurement_parent_adapter.py',)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def authenticate(path, retained_sha256, live_sources, root=ROOT):
    """Return the digest of the current parent receipt after checking its sources.

    ``retained_sha256`` is accepted for the consumers' call signature; a
    receipt is identified by its current bytes, not by a recorded digest.
    """
    root = Path(root).resolve()
    if Path(path).resolve() != root / PARENT:
        raise ValueError('unreviewed complete measurement parent')
    data = (root / PARENT).read_bytes()
    current = strict_json.loads(data.decode())
    if canonical(bound_pins(current.get('source_bindings'))) != canonical(bound_pins(live_sources)):
        raise ValueError('complete measurement current source map differs')
    for relative, expected in bound_pins(live_sources).items():
        source = (root / relative).resolve()
        if (not source.is_relative_to(root) or not source.is_file()
                or hashlib.sha256(source.read_bytes()).hexdigest() != expected):
            raise ValueError('complete measurement source bytes changed: ' + relative)
    return hashlib.sha256(data).hexdigest()
