"""Exact 114-factor matter input retained by the official construction."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SNAPSHOT_SHA256 = '456eae4b8b95801bc2213255d8df8448e383682be0c924f584cf70d6edfc13cc'

@lru_cache(None)
def snapshot():
    raw = (HERE / 'source_snapshot.json').read_bytes()
    if hashlib.sha256(raw).hexdigest() != SNAPSHOT_SHA256:
        raise ValueError('The specified native matter snapshot has changed')
    data = json.loads(raw)
    assert len(data['terms']) == 114 and data['native_CAR_modes'] == 192
    assert len(data['mode_charge_native_thirds']) == len(data['mode_vertex']) == 192
    return data
