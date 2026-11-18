"""Uniform edgewise cone refinements and the real quartic-wave subaction.

The analytic paper proves the conditional continuum trajectory estimate.
This finite packet checks refinement, Ritz approximation and an autonomous
nonlinear trajectory; numerical comparisons are not interval error bounds.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations, permutations, product
import json
from math import factorial, gcd
from pathlib import Path
import re

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import factorized, spsolve

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "runtime/whitney_real_continuum_receipt.json"
SCOPE = "UNIFORM_CONE_REFINEMENT__CONDITIONAL_REAL_QUARTIC_WAVE_CONTINUUM__FINITE_NUMERIC_CHECKS"
PINS = (
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "code/electromagnetism/whitney_interacting_quantum.py",
    "code/electromagnetism/whitney_real_continuum.py",
    "code/electromagnetism/verify_whitney_real_continuum.py",
    "code/electromagnetism/test_whitney_real_continuum.py",
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode("ascii")


def macro_mesh():
    text = (ROOT/PINS[0]).read_text(encoding="utf-8")
    text = text.split("def portVector", 1)[1].split("theorem portVector_positivePort", 1)[0]
    phi = (1+np.sqrt(5))/2
    values = {"0": 0., "1": 1., "-1": -1., "φ": phi, "-φ": -phi}
    vertices = np.array([[0., 0., 0.]]+[[values[x.strip()] for x in row.split(",")]
                        for row in re.findall(r"!\[([^\[\]]+)\]", text)])
    text = (ROOT/"Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    text = text.split("def orientedFaces", 1)[1].split("def faceEdges", 1)[0]
    cells = [(0,)+tuple(sorted(int(x)+1 for x in face))
             for face in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", text)]
    if vertices.shape != (13, 3) or len(set(cells)) != 20:
        raise ValueError("canonical cone census")
    return vertices, sorted(cells)


def vertex_key(macro, bary):
    divisor = gcd(*map(int, bary))
    return tuple((int(v), int(w)//divisor) for v, w in zip(macro, bary, strict=True) if w)


def mesh(n):
    """Restrict the n-grid Kuhn triangulation to cumulative barycentric space."""
    if type(n) is not int or n < 1:
        raise ValueError("positive integer refinement required")
    vertices, macros = macro_mesh()
    points, cells = {}, []
    steps = np.eye(3, dtype=int)
    for macro in macros:
        for base in product(range(n), repeat=3):
            for order in permutations(range(3)):
                cumulative = np.vstack((base, np.array(base)+np.cumsum(steps[list(order)], axis=0)))
                if not np.all(cumulative[:, 0] >= cumulative[:, 1]) or not np.all(cumulative[:, 1] >= cumulative[:, 2]):
                    continue
                bary = np.column_stack((n-cumulative[:, 0], cumulative[:, 0]-cumulative[:, 1],
                                        cumulative[:, 1]-cumulative[:, 2], cumulative[:, 2]))
                keys = [vertex_key(macro, row) for row in bary]
                for key in keys:
                    points[key] = sum(vertices[v]*w for v, w in key)/sum(w for _, w in key)
                cells.append(tuple(keys))
    keys = sorted(points)
    index = {key: i for i, key in enumerate(keys)}
    tets = sorted(tuple(sorted(index[key] for key in cell)) for cell in cells)
    return np.array([points[key] for key in keys]), np.array(tets, dtype=int), keys


def geometry(vertices, cells):
    x = vertices[cells]
    affine = np.concatenate((np.ones(x.shape[:-1]+(1,)), x), axis=-1)
    grad = np.linalg.inv(affine)[:, 1:, :].transpose(0, 2, 1)
    volume = abs(np.linalg.det(x[:, 1:]-x[:, :1]))/6
    return x, grad, volume


def matrices(vertices, cells):
    x, grad, volume = geometry(vertices, cells)
    rows = np.broadcast_to(cells[:, :, None], (len(cells), 4, 4)).ravel()
    cols = np.broadcast_to(cells[:, None, :], (len(cells), 4, 4)).ravel()
    mass = volume[:, None, None]*(np.ones((4, 4))+np.eye(4))[None]/20
    stiffness = np.einsum("tic,tjc,t->tij", grad, grad, volume)
    shape = (len(vertices), len(vertices))
    return (coo_matrix((mass.ravel(), (rows, cols)), shape=shape).tocsr(),
            coo_matrix((stiffness.ravel(), (rows, cols)), shape=shape).tocsr())


def quadrature(order=4):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes+1)/2, weights/2
    r, s, t = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    wr, ws, wt = np.meshgrid(weights, weights, weights, indexing="ij")
    lam = np.stack((1-r, r*(1-s), r*s*(1-t), r*s*t), axis=-1).reshape(-1, 4)
    weight = (wr*ws*wt*r*r*s).ravel()
    return lam, weight


def cubic_load(values, cells, volume):
    """Exact integral of (sum lambda_i u_i)^3 lambda_j, via h_3."""
    u = values[cells]
    s1, s2, s3 = (np.sum(u**p, axis=1)[:, None] for p in (1, 2, 3))
    h2 = (s1*s1+s2)/2
    h3 = (s1**3+3*s1*s2+2*s3)/6
    local = volume[:, None]*(h3+u*h2+u*u*s1+u**3)/140
    return np.bincount(cells.ravel(), weights=local.ravel(), minlength=len(values))


def energy(u, velocity, mass, stiffness, cells, volume, m2=.5, g=.25):
    # Integral u_h^4 = sum_i u_i integral(u_h^3 lambda_i).
    return float((velocity@mass@velocity+u@stiffness@u+m2*(u@mass@u))/2
                 +g*(u@cubic_load(u, cells, volume))/4)


def ritz_probe(vertices, cells):
    x, grad, volume = geometry(vertices, cells)
    mass, stiffness = matrices(vertices, cells)
    lam, weights = quadrature()
    q = np.einsum("qi,tic->tqc", lam, x)
    ref = np.sum(q*q, axis=2)
    local = np.einsum("tqc,tic,t,q->ti", 2*q, grad, 6*volume, weights)
    local += np.einsum("tq,qi,t,q->ti", ref, lam, 6*volume, weights)
    load = np.bincount(cells.ravel(), weights=local.ravel(), minlength=len(vertices))
    projected = spsolve(stiffness+mass, load)
    result = {}
    for name, coefficients in (("nodal", np.sum(vertices**2, axis=1)), ("ritz", projected)):
        value = coefficients[cells]@lam.T
        gradient = np.einsum("ti,tic->tc", coefficients[cells], grad)[:, None, :]
        l2sq = np.einsum("tq,t,q->", (value-ref)**2, 6*volume, weights)
        h1sq = l2sq+np.einsum("tqc,tqc,t,q->", gradient-2*q, gradient-2*q, 6*volume, weights)
        result[name+"_l2_error"] = float(np.sqrt(l2sq))
        result[name+"_h1_error"] = float(np.sqrt(h1sq))
    result["ritz_weak_residual"] = float(np.max(abs((stiffness+mass)@projected-load)))
    return result


def mesh_record(n):
    vertices, cells, keys = mesh(n)
    x, grad, volume = geometry(vertices, cells)
    diameters = np.max([np.linalg.norm(x[:, i]-x[:, j], axis=1) for i, j in combinations(range(4), 2)], axis=0)
    inradius = 1/np.linalg.norm(grad, axis=2).sum(axis=1)
    faces = Counter(tuple(sorted(face)) for cell in cells for face in combinations(cell, 3))
    return {"n": n, "vertices": len(vertices), "tetrahedra": len(cells),
        "boundary_triangles": sum(count == 1 for count in faces.values()),
        "interior_triangles": sum(count == 2 for count in faces.values()),
        "mesh_sha256": hashlib.sha256(canonical({"keys": keys, "cells": cells.tolist()})).hexdigest(),
        "h_min_times_n": float(diameters.min()*n), "h_max_times_n": float(diameters.max()*n),
        "inradius_min_times_n": float(inradius.min()*n),
        "diameter_inradius_max": float((diameters/inradius).max()),
        "volume": float(volume.sum()), **ritz_probe(vertices, cells)}


def trajectory(n=4, steps=32):
    vertices, cells, _ = mesh(n)
    _, _, volume = geometry(vertices, cells)
    mass, stiffness = matrices(vertices, cells)
    solve = factorized(mass.tocsc())
    u = .2*np.maximum(1-np.sum(vertices**2, axis=1), 0)**4
    velocity = np.zeros(len(vertices))
    initial = energy(u, velocity, mass, stiffness, cells, volume)
    def rhs(state):
        field, speed = state[:len(vertices)], state[len(vertices):]
        force = stiffness@field+.5*(mass@field)+.25*cubic_load(field, cells, volume)
        return np.concatenate((speed, -solve(force)))
    state = np.concatenate((u, velocity)); step = .125/steps
    for _ in range(steps):
        k1 = rhs(state); k2 = rhs(state+step*k1/2)
        k3 = rhs(state+step*k2/2); k4 = rhs(state+step*k3)
        state += step*(k1+2*k2+2*k3+k4)/6
    return {"n": n, "steps": steps, "time_interval": ["0", "1/8"],
        "initialization": "nodal compact real pulse; this run is not a continuum-error certificate",
        "final_state": state.tolist(), "initial_energy": initial,
        "final_energy": energy(state[:len(vertices)], state[len(vertices):], mass, stiffness, cells, volume)}


def build():
    return {"schema": "oph.whitney_real_continuum.v1", "scope": SCOPE,
        "source_pins": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in PINS},
        "analytic_scope": {"uniform_refinement_proved_in_paper": True,
            "conditional_real_sector_trajectory_bound": True, "full_charged_complex_trajectory_bound": False,
            "continuum_reference_existence_assumed": True, "numerical_trajectory_error_certified": False,
            "physical_source_or_clock_selected": False, "formalized_in_lean": False},
        "uniform_bounds": {"h_min_times_n_lower": "1/2", "h_max_times_n_upper": "6",
            "inradius_times_n_lower": "1/12", "shape_ratio_upper": "72"},
        "mesh_checks": [mesh_record(n) for n in (1, 2, 4, 8)],
        "trajectories": [trajectory(4, steps) for steps in (32, 64)]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    packet = build(); args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({"receipt": str(args.output), "bytes": args.output.stat().st_size}))
