#!/usr/bin/env python3
"""Prepare a family reference and ordered moments from complete native records.

The native memory and ray chart are source-bound. Quarter-turn measurement,
fresh auxiliary preparation, and ordered family feedback are declared logical
services. Every success and failure is retained. No empirical mass is read.
"""
from __future__ import annotations

import argparse
import ast
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys


import numpy as np
import sympy as sp

from .._source import PIN, ROOT, verify_sources
HERE = Path(__file__).resolve().parent
I3 = sp.eye(3)
SQRT5 = sp.sqrt(5)
PHASE_PAIRS = ((2, 3), (2, 10), (3, 10))
SOURCE_FILES = (
    "model/common_source/slow_geometry.py",
    "model/common_source/native_family_record_memory.py",
    "model/common_source/native_family_protocol.py",
    "model/common_source/native_family_protocol_results.npz",
    "model/common_source/native_repair_mass.py",
    "model/common_source/native_repair_mass_readout.py",
    "model/common_source/native_gauge_object.py",
    "model/observer/native.py",
    "model/native_modular_reference/README.md",
    "model/native_modular_reference/construction.py",
    "model/native_modular_reference/results.json",
    "model/native_modular_reference/preparation.py",
    "model/native_modular_reference/preparation_results.json",
    "model/native_modular_reference/full_flow.py",
    "model/native_modular_reference/full_flow_results.json",
)


def exact_scalar(x):
    return sp.expand(sp.radsimp(sp.simplify(x)))


def exact_matrix(x):
    return sp.ImmutableMatrix(sp.Matrix(x).applyfunc(exact_scalar))


def is_zero(x):
    return all(exact_scalar(v) == 0 for v in x)


def matrix_json(x):
    return [[str(exact_scalar(v)) for v in row] for row in sp.Matrix(x).tolist()]


def numeric(x):
    return np.asarray(sp.Matrix(x).evalf(30), dtype=np.complex128)


def exact_real_float(x):
    """Simplify a proved-real exact scalar before crossing the numeric API."""
    value = exact_scalar(x)
    if exact_scalar(sp.im(value)) != 0:
        raise ArithmeticError("A Hermitian moment has a nonzero imaginary part")
    return float(exact_scalar(sp.re(value)))


def _pin(repo):
    return verify_sources(repo, SOURCE_FILES)


@lru_cache(maxsize=4)
def _exact_data(repo_string):
    repo = Path(repo_string)
    hashes = _pin(repo)
    path = repo / "model/common_source/slow_geometry.py"
    tree = ast.parse(path.read_text())
    fn = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == "raw_frame")
    ret = next(x.value for x in fn.body if isinstance(x, ast.Return))
    if not all(isinstance(x, (ast.Tuple, ast.Name, ast.Load)) for x in ast.walk(ret)):
        raise ValueError("Source frame does not have the pinned literal grammar")
    phi = (1 + SQRT5) / 2
    names = {"zero": sp.Integer(0), "one": sp.Integer(1), "minus": sp.Integer(-1),
             "phi": phi, "negative_phi": -phi}
    frame = sp.Matrix(eval(compile(ast.Expression(ret), str(path), "eval"),
                           {"__builtins__": {}}, names))
    if frame.shape != (12, 3):
        raise ValueError("Native chart dimension differs")
    projectors = tuple(exact_matrix(frame[j, :].T * frame[j, :] /
                                    (frame[j, :] * frame[j, :].T)[0])
                       for j in range(12))
    assert all(is_zero(P * P - P) and sp.trace(P) == 1 for P in projectors)
    assert is_zero(sum(projectors, sp.zeros(3)) - 4 * I3)
    seams = tuple((i, j) for i in range(12) for j in range(i+1, 12)
                  if exact_scalar(sum((frame[i, k] - frame[j, k])**2
                                      for k in range(3))) == 4)
    assert len(seams) == 30
    endpoint_counts = [sum(int(j in edge) for edge in seams) for j in range(12)]
    assert endpoint_counts == [5] * 12

    with np.load(repo / "model/common_source/native_family_protocol_results.npz",
                 allow_pickle=False) as packet:
        raw_labels = packet["memory_raw_labels"].copy()
        raw_weights = packet["memory_raw_weights"].copy()
        labels = packet["memory_labels"].copy()
        recorded_T = packet["memory_transition"].copy()
        source_pulses = packet["memory_pulses"].copy()
    expected_labels = np.array([(cc, cr, ec, er) for cc in (0, 1) for cr in range(4)
                                for ec in (0, 1) for er in range(8)], dtype=np.int64)
    assert np.array_equal(raw_labels, expected_labels)
    assert raw_weights.shape == (128, 2) and labels.shape == (1536, 6)
    weights = tuple(sp.Rational(int(n), int(d)) for n, d in raw_weights)
    assert all(x >= 0 for x in weights) and sum(weights) == 1
    T = sp.zeros(12)
    destinations = np.empty((128, 12), dtype=np.int64)
    for r, row in enumerate(raw_labels):
        for k in range(12):
            label = labels[12*r+k]
            assert np.array_equal(label[:4], row) and label[4] == k
            j = int(label[5])
            destinations[r, k] = j
            T[j, k] += weights[r]
    assert np.max(np.abs(numeric(T).real - recorded_T)) < 2e-16
    assert all(sum(T[:, k]) == 1 for k in range(12))
    null = (T-sp.eye(12)).nullspace()
    assert len(null) == 1
    pi = exact_matrix(null[0] / sum(null[0]))
    assert T*pi == pi and sum(pi) == 1 and all(x >= 0 for x in pi)
    support = tuple(k for k in range(12) if pi[k] > 0)
    assert support == (2, 3, 9, 10)
    assert min(T[2, k] for k in range(12)) == sp.Rational(11, 12)
    assert min(T[2, k] for k in support) == sp.Rational(67, 72)
    p = sum(w for row, w in zip(raw_labels, weights) if tuple(row[:2]) == (0, 1))
    q = sum(w for row, w in zip(raw_labels, weights) if row[2] == 0 and row[3] in (4, 5))
    pq = sum(w for row, w in zip(raw_labels, weights)
             if tuple(row[:2]) == (0, 1) and row[2] == 0 and row[3] in (4, 5))
    assert (p, q, pq) == (sp.Rational(1, 12), sp.Rational(1, 6), sp.Rational(1, 72))
    pulse_error = max(float(np.max(np.abs(
        -source_pulses[a, :, :3, 3:].imag / np.sin(h) -
        np.asarray([numeric(P) for P in projectors]))))
        for a, h in enumerate((1/8, 1/4)))
    assert pulse_error < 2e-15

    Y = exact_matrix(sum((pi[k]*projectors[k] for k in range(12)), sp.zeros(3)))
    W = exact_matrix(T * sp.diag(*pi))
    jordan = sp.zeros(3)
    pinched = sp.zeros(3)
    commutator = sp.zeros(3)
    pinched_reverse = sp.zeros(3)
    for j in range(12):
        for k in range(12):
            P, Q = projectors[j], projectors[k]
            if W[j, k]:
                jordan += W[j, k] * (P*Q + Q*P)
                pinched += W[j, k] * P*Q*P
                commutator += W[j, k] * (P*Q - Q*P)
            if W[k, j]:
                pinched_reverse += W[k, j] * P*Q*P
    D = exact_matrix(2*pinched-jordan)
    J = exact_matrix(-sp.I*commutator)
    assert is_zero(pinched-pinched_reverse)
    assert is_zero(D-D.T) and sp.trace(D) == 0
    assert is_zero(J-J.H) and is_zero(J+J.conjugate()) and sp.trace(J) == 0
    C = tuple(exact_matrix(sp.I*(projectors[a]*projectors[b]-projectors[b]*projectors[a]))
              for a, b in PHASE_PAIRS)
    F = tuple(exact_matrix(projectors[k]-I3/3) for k in (0, 1, 4, 5, 8)) + C
    detY = exact_scalar(Y.det())
    assert detY.is_positive is True
    return dict(repo=repo, source_hashes=hashes, frame=exact_matrix(frame),
                projectors=projectors, seams=seams, raw_labels=raw_labels,
                weights=weights, destinations=destinations, T=exact_matrix(T), pi=pi,
                support=support, Y=Y, W=W, D=D, J=J, features=F, phase_features=C,
                pulse_replay_error=pulse_error, detY=detY)


def exact_source_packet(repo):
    """Exact SymPy source data. Treat returned matrices as immutable."""
    resolved = Path(repo).resolve()
    _pin(resolved)
    return _exact_data(str(resolved))
