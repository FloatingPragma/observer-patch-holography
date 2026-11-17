#!/usr/bin/env python3
"""Exact per-vertex volume, capacity and memory matrix columns."""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sy

try:
    from .geometry import GeometryConfig, build_geometry_extension
except ImportError:
    from .geometry import GeometryConfig, build_geometry_extension

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def _exact_snapshot(value):
    z = complex(value)
    return sy.Rational(z.real)+sy.I*sy.Rational(z.imag)


@lru_cache(maxsize=16)
def auxiliary_columns(config=GeometryConfig()):
    """Complete auxiliary columns, exact for the frozen numeric snapshots."""
    extension = build_geometry_extension(np.zeros((1, 1)), config=config)
    matrices = []
    for matrix in (extension.auxiliary_hamiltonian, extension.lapse):
        csc = matrix.tocsc()
        columns = []
        for column in range(csc.shape[1]):
            lo, hi = csc.indptr[column:column+2]
            columns.append(tuple((int(row), _exact_snapshot(value))
                for row, value in zip(csc.indices[lo:hi], csc.data[lo:hi]) if value != 0))
        matrices.append(tuple(columns))
    return tuple(matrices)

