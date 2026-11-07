#!/usr/bin/env python3
"""A complete finite chiral reference/propagator/measure patch.

The patch is local in gauge-configuration space. Spacetime locality of its
measure current, fixed-reference gauge invariance, the quantum gauge law,
reflection positivity and an interacting continuum are not inferred.
"""
from __future__ import annotations

import argparse
from collections import Counter
import importlib
from itertools import product
import json
from pathlib import Path
import platform
import sys


import numpy as np
import scipy


SITES = tuple(product(range(2), repeat=4))
POSITIONS = {point: i for i, point in enumerate(SITES)}
DELTA = 1/2500
TARGET = (.4, .6)
TORON = np.array([.17, -.11, .07, .03])


def small(value, name, bound=2e-10):
    value = float(value)
    if not np.isfinite(value) or value < 0 or value > bound:
        raise ArithmeticError(f"{name}: {value} > {bound}")
    return value


def adjoint(a):
    return a.conj().T


def complex_record(z):
    return {"real": float(complex(z).real), "imag": float(complex(z).imag)}


def inverse_sqrt(a):
    a = (a+adjoint(a))/2
    eig, vec = np.linalg.eigh(a)
    if eig[0] <= 0:
        raise ValueError("A positive matrix is required")
    root = (vec*eig**(-.5))@adjoint(vec)
    return root, eig, vec


def canonical_transport(p, p0):
    """Canonical unitary mapping ran(P0) to ran(P), when ||P-P0||<1."""
    n = len(p)
    identity = np.eye(n)
    c = p@p0+(identity-p)@(identity-p0)
    metric = adjoint(c)@c
    inv, values, vectors = inverse_sqrt(metric)
    if values[0] < 1e-5:
        raise ValueError("This coordinate patch approaches its singular boundary")
    q = c@inv
    small(np.linalg.norm(adjoint(q)@q-identity, "fro"), "Canonical transport unitarity")
    small(np.linalg.norm(p@q-q@p0, "fro"), "Canonical projector intertwining")
    return q, {"c": c, "inverse_sqrt": inv, "metric_eigenvalues": values,
               "metric_eigenvectors": vectors}


def canonical_transport_variation(p0, dp, aux):
    """Analytic Frechet derivative, including the degenerate metric spectrum."""
    c, inv = aux["c"], aux["inverse_sqrt"]
    values, vec = aux["metric_eigenvalues"], aux["metric_eigenvectors"]
    dc = dp@(2*p0-np.eye(len(p0)))
    da = adjoint(dc)@c+adjoint(c)@dc
    roots = np.sqrt(values)
    divided = -1/(roots[:, None]*roots[None, :]*(roots[:, None]+roots[None, :]))
    dinv = vec@(divided*(adjoint(vec)@da@vec))@adjoint(vec)
    return dc@inv+c@dinv


def projector_basis(p):
    values, vectors = np.linalg.eigh((p+adjoint(p))/2)
    basis = vectors[:, values > .5]
    small(np.linalg.norm(basis@adjoint(basis)-p, "fro"), "Complete reference range")
    return basis


def chiral_frame(d, p, p0, v0, vbar):
    q, aux = canonical_transport(p, p0)
    v = q@v0
    m = adjoint(vbar)@d@v
    singular = np.linalg.svd(m, compute_uv=False)
    if singular[-1] <= 1e-5:
        raise ValueError("This Gaussian response patch requires an invertible Weyl matrix")
    green = v@np.linalg.solve(m, adjoint(vbar))
    pbar = vbar@adjoint(vbar)
    errors = {
        "frame_orthonormality": small(np.linalg.norm(adjoint(v)@v-np.eye(v.shape[1]), "fro"), "Chiral frame"),
        "complete_chiral_range": small(np.linalg.norm(v@adjoint(v)-p, "fro"), "Complete chiral range"),
        "left_inverse_identity": small(np.linalg.norm(d@green-pbar, "fro"), "D G = Pbar"),
        "right_inverse_identity": small(np.linalg.norm(green@d-p, "fro"), "G D = P"),
        "left_chiral_support": small(np.linalg.norm(p@green-green, "fro"), "P G = G"),
        "right_chiral_support": small(np.linalg.norm(green@pbar-green, "fro"), "G Pbar = G"),
    }
    return {"Q": q, "V": v, "M": m, "G": green, "aux": aux,
            "errors": errors, "minimum_Weyl_singular_value": float(singular[-1])}


def control_fields():
    a = np.empty((16, 4)); b = np.empty_like(a); h = np.empty_like(a)
    for i in range(16):
        for mu in range(4):
            a[i, mu] = (((11*i+7*mu+3) % 17)-8)/8
            b[i, mu] = (((13*i+5*mu+1) % 19)-9)/9
            h[i, mu] = (((7*i+11*mu+2) % 23)-11)/11
    return a, b, h


FIELDS = control_fields()


def link_control(charge, parameters=None):
    theta = np.broadcast_to(TORON, (16, 4)).copy()
    if parameters is not None:
        ss, tt = parameters
        theta += DELTA*(FIELDS[2]+ss*FIELDS[0]+tt*FIELDS[1])
    links = np.exp(1j*charge*theta)
    # A fermion spin structure, not an extra gauge field. One temporal
    # boundary crossing is antiperiodic; the plaquette signs cancel.
    for i, x in enumerate(SITES):
        if x[3] == 1:
            links[i, 3] *= -1
    return links


def wilson_independent(links, gamma, *, derivative=False):
    n = 64
    out = np.zeros((n, n), complex) if derivative else 4*np.eye(n, dtype=complex)
    for i, x in enumerate(SITES):
        span = slice(4*i, 4*i+4)
        for mu in range(4):
            y = tuple(1-z if j == mu else z for j, z in enumerate(x))
            iy = POSITIONS[y]
            other = slice(4*iy, 4*iy+4)
            out[span, other] += -.5*(np.eye(4)-gamma[mu])*links[i, mu]
            out[other, span] += -.5*(np.eye(4)+gamma[mu])*links[i, mu].conjugate()
    return out


def overlap_data(source, charge, parameters=None, *, links_override=None):
    links = link_control(charge, parameters) if links_override is None else links_override
    w, g5 = source.wilson_operator(links[:, :, None, None], side=2)
    gamma, _ = source.gamma_matrices()
    small(np.linalg.norm(w-wilson_independent(links, gamma), "fro"), "Independent complete Wilson stencil")
    d, hat, gap = source.overlap(w, g5)
    h = g5@(w-np.eye(64))
    eig, vec = np.linalg.eigh(h)
    p = (np.eye(64)-hat)/2
    pbar = (np.eye(64)+g5)/2
    small(abs(np.trace(p).real-32), "Full chiral rank")
    gw = small(np.linalg.norm(g5@d+d@g5-d@g5@d, "fro"), "GW")
    plaquette = source.plaquette_defect(links[:, :, None, None], side=2)
    if plaquette >= 1/40:
        raise ValueError("Declared background left the source admissible plaquette domain")
    return {"links": links, "W": w, "gamma5": g5, "D": d, "P": p,
            "Pbar": pbar, "H": h, "eigenvalues": eig, "eigenvectors": vec,
            "gap": gap, "GW": gw, "plaquette": plaquette}


def overlap_variation(data, dlinks, gamma):
    dw = wilson_independent(dlinks, gamma, derivative=True)
    dh = data["gamma5"]@dw
    small(np.linalg.norm(dh-adjoint(dh), "fro"), "Hermitian sign-kernel variation")
    eig, vec = data["eigenvalues"], data["eigenvectors"]
    signs = np.where(eig > 0, 1., -1.)
    divided = np.zeros((64, 64))
    opposite = signs[:, None] != signs[None, :]
    denominator = eig[:, None]-eig[None, :]
    numerator = signs[:, None]-signs[None, :]
    divided[opposite] = numerator[opposite]/denominator[opposite]
    ds = vec@(divided*(adjoint(vec)@dh@vec))@adjoint(vec)
    return data["gamma5"]@ds, ds/2
