"""Small helpers shared by the modules: plain JSON values and timing."""
from __future__ import annotations

import math
import time

import numpy as np


def clean(obj):
    """Convert NumPy scalars, arrays and tuples into plain JSON values; reject non-finite floats."""
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [clean(v) for v in obj.tolist()]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        x = float(obj)
        if not math.isfinite(x):
            raise ValueError('non-finite value in a result')
        return x
    if obj is None or isinstance(obj, str):
        return obj
    raise TypeError(f'unexpected value of type {type(obj).__name__}')


class Clock:
    """Wall time of the sections of a run (reported, never compared)."""

    def __init__(self):
        self.start = time.monotonic()
        self.sections = {}

    def mark(self, name):
        now = time.monotonic()
        self.sections[name] = round(now - self.start - sum(self.sections.values()), 2)
