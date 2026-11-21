"""Read a published acquisition or explicitly create one at an unused path."""
from __future__ import annotations

import hashlib
import io
from pathlib import Path

import numpy as np


def snapshot(path, digest, *, arrays=None):
    path = Path(path)
    if arrays is None:
        blob = path.read_bytes()
        if hashlib.sha256(blob).hexdigest() != digest:
            raise ValueError(f"changed published acquisition: {path.name}")
    else:
        buf = io.BytesIO()
        np.savez_compressed(buf, **arrays)
        blob = buf.getvalue()
        with path.open("xb") as output:
            output.write(blob)
    return blob


def finite_float_record(value, expected, tolerance):
    return (value.dtype == expected.dtype and value.shape == expected.shape
            and np.isfinite(value).all() and np.max(abs(value-expected)) <= tolerance)
