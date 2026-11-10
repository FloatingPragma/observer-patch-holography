"""Typed live custody for the retained complete 66-word measurement replay.

The original receipt is immutable. The current receipt was produced by the
unchanged measurement producer; its scientific JSON is exactly identical.
Only its source hashes and the four reviewed inherited source paths differ.
This adapter performs no measurement, spectral or clock calculation.
"""
import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
PARENT = 'model/native_energy_measurement/results.json'
ARCHIVE = 'model/completion/historical/r2063_live_inputs/' + PARENT
RETAINED_SHA256 = '6c4f66f29e34eac29bac5d7ebe1c857606888eede03113761556c41eb1b0ffeb'
REPLAYED_SHA256 = '5666787d48cdb4e96112cc5ec97b438767499e06076c120208597cdcd90f4cc6'
SOURCE_FILES = (
    'model/completion/native_measurement_parent_adapter.py',
    'model/completion/test_native_measurement_parent_adapter.py',
)
REVIEWED_APPEND = {
    'Lean/SourceLaw/RecordRealization.lean': '26dfd52ccd75e29d3b3f518719ba370dfc502ba350c5b8b41de6723a659b7864',
    'model/native_record_response/SHORT_RECORD_RECONSTRUCTION.md': '194f10986e81a096be290ba615b1d0f5b72f99d39e4f3c50bf7adc6b81b1c98d',
    'model/native_record_response/short_records.py': '3e9ecafdc322d386afb83532113ea4b55466b5be6655a4eff325692101cc161e',
    'model/native_record_response/test_short_records.py': '3dde9821de4b6065c8bd65f7f67e4eb94210f9396fd30d6ff0181287e9a40231',
}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def compare_payloads(retained, current, live_sources):
    """Retain all typed science and admit exactly the reviewed input union."""
    if (type(retained) is not dict or type(current) is not dict
            or set(retained) != set(current)):
        raise ValueError('complete measurement field inventory changed')
    old = retained['source_bindings']
    new = current['source_bindings']
    if set(new) != set(old) | set(REVIEWED_APPEND):
        raise ValueError('complete measurement input scope changed')
    if any(new.get(path) != expected for path, expected in REVIEWED_APPEND.items()):
        raise ValueError('reviewed complete measurement input changed')
    if canonical(new) != canonical(live_sources):
        raise ValueError('complete measurement current source map differs')
    left = {key: value for key, value in retained.items() if key != 'source_bindings'}
    right = {key: value for key, value in current.items() if key != 'source_bindings'}
    if canonical(left) != canonical(right):
        raise ValueError('complete measurement scientific payload changed')


def authenticate(path, retained_sha256, live_sources, root=ROOT):
    root = Path(root).resolve()
    if Path(path).resolve() != root / PARENT or retained_sha256 != RETAINED_SHA256:
        raise ValueError('unreviewed complete measurement parent')
    old_bytes = (root / ARCHIVE).read_bytes()
    current_bytes = (root / PARENT).read_bytes()
    if hashlib.sha256(old_bytes).hexdigest() != RETAINED_SHA256:
        raise ValueError('retained complete measurement bytes changed')
    if hashlib.sha256(current_bytes).hexdigest() != REPLAYED_SHA256:
        raise ValueError('current complete measurement differs from the actual replay')
    retained = strict_json.loads(old_bytes.decode())
    current = strict_json.loads(current_bytes.decode())
    compare_payloads(retained, current, live_sources)
    for relative, expected in live_sources.items():
        source = (root / relative).resolve()
        if (not source.is_relative_to(root) or not source.is_file()
                or hashlib.sha256(source.read_bytes()).hexdigest() != expected):
            raise ValueError('complete measurement source bytes changed: ' + relative)
    return REPLAYED_SHA256
