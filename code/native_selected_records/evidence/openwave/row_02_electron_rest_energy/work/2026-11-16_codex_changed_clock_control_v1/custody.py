"""Validate the prospectively frozen code and source-input bytes, without writes."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
OPENWAVE=HERE.parents[2]
def check():
    record=json.loads((HERE/'PREREG.json').read_text())
    for name,digest in record['local_files'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Frozen local file changed: '+name)
    for name,digest in record['source_inputs'].items():
        if hashlib.sha256((OPENWAVE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Frozen source input changed: '+name)
    return record
