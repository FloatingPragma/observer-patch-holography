"""Independent finite verification of the real-matter continuum packet.

No producer is imported. Ordered lattice chains reconstruct the mesh;
exact simplex moments assemble the Ritz load and cubic force. Adaptive
DOP853 independently checks the fixed-mesh RK4 endpoints. Analytic
continuum convergence and uniform mesh quality remain paper proofs.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from math import comb, factorial, gcd, isfinite
from pathlib import Path
import re

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import factorized, spsolve
from scipy.special import roots_jacobi

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent/"runtime/whitney_real_continuum_receipt.json"
SCOPE = "UNIFORM_CONE_REFINEMENT__CONDITIONAL_REAL_QUARTIC_WAVE_CONTINUUM__FINITE_NUMERIC_CHECKS"
PINS = {
    "Lean/Screen/SeamCurrentEdge30Moment.lean", "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "code/electromagnetism/whitney_interacting_quantum.py",
    "code/electromagnetism/whitney_real_continuum.py",
    "code/electromagnetism/verify_whitney_real_continuum.py", "code/electromagnetism/test_whitney_real_continuum.py",
}
ANALYTIC = {"uniform_refinement_proved_in_paper": True,
    "conditional_real_sector_trajectory_bound": True, "full_charged_complex_trajectory_bound": False,
    "continuum_reference_existence_assumed": True, "numerical_trajectory_error_certified": False,
    "physical_source_or_clock_selected": False, "formalized_in_lean": False}
BOUNDS = {"h_min_times_n_lower": "1/2", "h_max_times_n_upper": "6",
          "inradius_times_n_lower": "1/12", "shape_ratio_upper": "72"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode("ascii")


def exact(actual, expected, message):
    require(canonical(actual) == canonical(expected), message)


def finite(value):
    require(type(value) in (float, int) and isfinite(value), "finite numeric diagnostic")
    return float(value)


def load(path=OUTPUT):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result
    def number(value):
        result = float(value)
        require(isfinite(result), "nonfinite JSON number")
        return result
    def constant(_):
        raise ValueError("nonfinite JSON constant")
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs,
                      parse_float=number, parse_constant=constant)


def source_geometry():
    text = (ROOT/"Lean/Screen/SeamCurrentEdge30Moment.lean").read_text(encoding="utf-8")
    text = text.split("def portVector", 1)[1].split("theorem portVector_positivePort", 1)[0]
    phi = (1+np.sqrt(5))/2
    values = {"0": 0., "1": 1., "-1": -1., "φ": phi, "-φ": -phi}
    vertices = np.array([[0., 0., 0.]]+[[values[x.strip()] for x in row.split(",")]
                         for row in re.findall(r"!\[([^\[\]]+)\]", text)])
    text = (ROOT/"Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    text = text.split("def orientedFaces", 1)[1].split("def faceEdges", 1)[0]
    faces = [tuple(sorted(int(i)+1 for i in row)) for row in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", text)]
    require(vertices.shape == (13, 3) and len(set(faces)) == 20, "source geometry census")
    return vertices, [(0,)+face for face in sorted(faces)]


def reference_chains(n):
    cells = []
    for a in range(n+1):
        for b in range(a+1):
            for c in range(b+1):
                base = np.array([a, b, c])
                neighbors = []
                for increment in product((0, 1), repeat=3):
                    if not any(increment):
                        continue
                    point = base+increment
                    if n >= point[0] >= point[1] >= point[2] >= 0:
                        neighbors.append(point)
                for triple in combinations(neighbors, 3):
                    ordered = sorted(triple, key=lambda point: int(point.sum()))
                    if all(np.all(left <= right) for left, right in zip([base]+ordered[:-1], ordered)):
                        cell = np.vstack((base, ordered))
                        require(abs(round(np.linalg.det(cell[1:]-cell[0]))) == 1, "unimodular lattice chain")
                        cells.append(cell)
    require(len(cells) == n**3, "reference cell count")
    return cells


def mesh(n):
    original, macros = source_geometry()
    points, cells = {}, []
    chains = reference_chains(n)
    for macro in macros:
        for chain in chains:
            tet = []
            for a, b, c in chain:
                weights = (n-a, a-b, b-c, c)
                common = gcd(*map(int, weights))
                key = tuple((v, int(w)//common) for v, w in zip(macro, weights) if w)
                points[key] = sum(original[v]*w for v, w in key)/sum(w for _, w in key)
                tet.append(key)
            cells.append(tet)
    keys = sorted(points); lookup = {key: i for i, key in enumerate(keys)}
    cells = np.array(sorted(tuple(sorted(lookup[key] for key in tet)) for tet in cells), dtype=int)
    return np.array([points[key] for key in keys]), cells, keys


def assemble(vertices, cells):
    x = vertices[cells]
    edge = (x[:, 1:]-x[:, :1]).transpose(0, 2, 1)
    tail = np.linalg.inv(edge)
    grad = np.concatenate((-tail.sum(axis=1)[:, None, :], tail), axis=1)
    volume = abs(np.linalg.det(edge))/6
    rows, cols, masses, stiffness = [], [], [], []
    for i in range(4):
        for j in range(4):
            rows.extend(cells[:, i]); cols.extend(cells[:, j])
            masses.extend(volume*(1+int(i == j))/20)
            stiffness.extend(volume*np.einsum("tc,tc->t", grad[:, i], grad[:, j]))
    shape = (len(vertices), len(vertices))
    mass = coo_matrix((masses, (rows, cols)), shape=shape).tocsr()
    stiff = coo_matrix((stiffness, (rows, cols)), shape=shape).tocsr()
    return x, grad, volume, mass, stiff


def moment(indices):
    counts = Counter(indices)
    return Fraction(6*np.prod([factorial(count) for count in counts.values()]), factorial(len(indices)+3))


def ritz_check(vertices, cells, data):
    x, grad, volume, mass, stiffness = data
    cubic = np.array([float(moment(indices)) for indices in product(range(4), repeat=3)]).reshape(4, 4, 4)
    gram = np.einsum("tic,tjc->tij", x, x)
    local = 2*volume[:, None]*np.einsum("tc,tic->ti", x.mean(axis=1), grad)
    local += volume[:, None]*np.einsum("ijk,tjk->ti", cubic, gram)
    load = np.bincount(cells.ravel(), weights=local.ravel(), minlength=len(vertices))
    projected = spsolve(mass+stiffness, load)
    nodes, weights = [], []
    for power in (2, 1, 0):
        node, weight = roots_jacobi(4, 0, power)
        nodes.append((node+1)/2); weights.append(weight/2**(power+1))
    r, s, t = np.meshgrid(*nodes, indexing="ij")
    wr, ws, wt = np.meshgrid(*weights, indexing="ij")
    lam = np.stack((1-r, r*(1-s), r*s*(1-t), r*s*t), axis=-1).reshape(-1, 4)
    weight = (wr*ws*wt).ravel()
    q = np.einsum("qi,tic->tqc", lam, x); reference = (q*q).sum(axis=2)
    result = {}
    for name, coefficient in (("nodal", (vertices*vertices).sum(axis=1)), ("ritz", projected)):
        error = coefficient[cells]@lam.T-reference
        derivative = np.einsum("ti,tic->tc", coefficient[cells], grad)[:, None, :]-2*q
        l2 = np.einsum("tq,t,q->", error*error, 6*volume, weight)
        h1 = l2+np.einsum("tqc,tqc,t,q->", derivative, derivative, 6*volume, weight)
        result[name+"_l2_error"] = float(np.sqrt(l2)); result[name+"_h1_error"] = float(np.sqrt(h1))
    result["ritz_weak_residual"] = float(np.max(abs((mass+stiffness)@projected-load)))
    return result


def verify(packet):
    require(type(packet) is dict and set(packet) == {"schema", "scope", "source_pins", "analytic_scope", "uniform_bounds", "mesh_checks", "trajectories"}, "packet schema")
    exact(packet["schema"], "oph.whitney_real_continuum.v1", "schema version")
    exact(packet["scope"], SCOPE, "scientific scope")
    exact(packet["analytic_scope"], ANALYTIC, "analytic versus finite evidence scope")
    exact(packet["uniform_bounds"], BOUNDS, "exact uniform constants")
    require(type(packet["source_pins"]) is dict and set(packet["source_pins"]) == PINS, "source pin census")
    for relative in PINS:
        exact(packet["source_pins"][relative], hashlib.sha256((ROOT/relative).read_bytes()).hexdigest(), "source pin: "+relative)
    require(type(packet["mesh_checks"]) is list and len(packet["mesh_checks"]) == 4, "mesh census")
    record_keys = {"n", "vertices", "tetrahedra", "boundary_triangles", "interior_triangles", "mesh_sha256",
        "h_min_times_n", "h_max_times_n", "inradius_min_times_n", "diameter_inradius_max", "volume",
        "nodal_l2_error", "nodal_h1_error", "ritz_l2_error", "ritz_h1_error", "ritz_weak_residual"}
    for row, n in zip(packet["mesh_checks"], (1, 2, 4, 8), strict=True):
        require(type(row) is dict and set(row) == record_keys, "mesh record schema")
        for key, expected in {"n": n, "vertices": 13+42*(n-1)+50*comb(n-1, 2)+20*comb(n-1, 3),
                              "tetrahedra": 20*n**3, "boundary_triangles": 20*n*n,
                              "interior_triangles": 40*n**3-10*n*n}.items():
            exact(row[key], expected, "exact mesh count: "+key)
    require(type(packet["trajectories"]) is list and len(packet["trajectories"]) == 2, "trajectory census")
    diagnostic, stored = [], {}
    for row in packet["mesh_checks"]:
        n = row["n"]; vertices, cells, keys = mesh(n)
        exact(row["mesh_sha256"], hashlib.sha256(canonical({"keys": keys, "cells": cells.tolist()})).hexdigest(), "independent mesh hash")
        faces = Counter(tuple(sorted(face)) for cell in cells for face in combinations(cell, 3))
        require(all(value in (1, 2) for value in faces.values()), "conforming triangle multiplicity")
        require(sum(value == 1 for value in faces.values()) == row["boundary_triangles"]
                and sum(value == 2 for value in faces.values()) == row["interior_triangles"], "actual face census")
        _, macros = source_geometry()
        outer_faces = [set(macro[1:]) for macro in macros]
        for face, count in faces.items():
            if count == 1:
                support = {vertex for index in face for vertex, _ in keys[index]}
                require(any(support <= outer for outer in outer_faces), "no exposed internal macro-face")
        data = assemble(vertices, cells); x, grad, volume, _, _ = data
        h = np.max([np.linalg.norm(x[:, i]-x[:, j], axis=1) for i, j in combinations(range(4), 2)], axis=0)
        radius = 1/np.linalg.norm(grad, axis=2).sum(axis=1)
        actual = {"h_min_times_n": float(n*h.min()), "h_max_times_n": float(n*h.max()),
                  "inradius_min_times_n": float(n*radius.min()), "diameter_inradius_max": float((h/radius).max()),
                  "volume": float(volume.sum()), **ritz_check(vertices, cells, data)}
        require(actual["h_min_times_n"] > .5 and actual["h_max_times_n"] < 6
                and actual["inradius_min_times_n"] > 1/12 and actual["diameter_inradius_max"] < 72, "finite geometric bounds")
        for key, value in actual.items():
            require(abs(finite(row[key])-value) <= 3e-10*(1+abs(value)), "independent geometric/Ritz diagnostic: "+key)
        require(actual["ritz_weak_residual"] < 1e-10 and actual["ritz_h1_error"] <= actual["nodal_h1_error"]+1e-11, "Ritz projection identity")
        diagnostic.append(actual)
        if n == 4:
            stored = {"vertices": vertices, "cells": cells, "data": data}
    vertices, cells = stored["vertices"], stored["cells"]
    _, _, volume, mass, stiffness = stored["data"]
    fourth = np.array([float(moment(indices)) for indices in product(range(4), repeat=4)]).reshape((4,)*4)
    def cubic(values):
        u = values[cells]
        local = np.einsum("ijkl,tj,tk,tl,t->ti", fourth, u, u, u, volume, optimize=True)
        return np.bincount(cells.ravel(), weights=local.ravel(), minlength=len(vertices))
    solve = factorized(mass.tocsc()); size = len(vertices)
    initial = np.r_[.2*np.maximum(1-(vertices*vertices).sum(axis=1), 0)**4, np.zeros(size)]
    def rhs(_, state):
        u, speed = state[:size], state[size:]
        return np.r_[speed, -solve(stiffness@u+.5*(mass@u)+.25*cubic(u))]
    def energy(state):
        u, speed = state[:size], state[size:]
        return float((speed@mass@speed+u@stiffness@u+.5*(u@mass@u))/2+.25*(u@cubic(u))/4)
    independent = solve_ivp(rhs, (0, .125), initial, method="DOP853", rtol=2e-12, atol=2e-13)
    require(independent.success, "independent wave integration")
    trajectory_errors = []
    for row, steps in zip(packet["trajectories"], (32, 64), strict=True):
        require(type(row) is dict and set(row) == {"n", "steps", "time_interval", "initialization", "final_state", "initial_energy", "final_energy"}, "trajectory record schema")
        for key, value in {"n": 4, "steps": steps, "time_interval": ["0", "1/8"],
                           "initialization": "nodal compact real pulse; this run is not a continuum-error certificate"}.items():
            exact(row[key], value, "trajectory declared inputs")
        require(type(row["final_state"]) is list and len(row["final_state"]) == 2*size, "state dimension")
        state = np.array([finite(value) for value in row["final_state"]])
        error = float(np.max(abs(state-independent.y[:, -1])))
        require(error < 2e-7, "independent nonlinear wave endpoint")
        require(abs(finite(row["initial_energy"])-energy(initial)) < 1e-12, "initial nonlinear energy")
        require(abs(finite(row["final_energy"])-energy(state)) < 1e-12, "final nonlinear energy")
        require(abs(energy(state)-energy(initial)) < 2e-9, "numerical energy drift")
        trajectory_errors.append(error)
    require(trajectory_errors[1] < trajectory_errors[0]/8, "RK4 temporal refinement control")
    return {"accepted": True, "scope": SCOPE, "refinement_parameters": [1, 2, 4, 8],
        "tetrahedra": [20, 160, 1280, 10240], "uniform_shape_bound": "72",
        "conditional_real_sector_trajectory_bound": True, "full_charged_complex_trajectory_bound": False,
        "numerical_trajectory_error_certified": False, "physical_source_or_clock_selected": False,
        "formalized_in_lean": False, "numeric_diagnostics": {"meshes": diagnostic, "endpoint_errors": trajectory_errors}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), sort_keys=True))
