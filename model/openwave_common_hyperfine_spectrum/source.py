"""One positive native action with seven declared calibration inputs."""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import sympy as sy
from model.openwave_spin_gradient_nuclei import source as native
from model.openwave_native_composites.spectrum import one_body
from model.openwave_calibrated_hadrons.source import angular4, integral

HERE = Path(__file__).resolve().parent
INPUTS = json.loads((HERE / 'inputs.json').read_text())


def parameters():
    return {k: F(v) for k, v in INPUTS['declared_rational_parameters'].items()}


def weights(p):
    return {**native.weights(p), 'Pair4': p['pair_triplet'] / 4}


def species_spin4(modes, basis, species):
    index = {m: i for i, m in enumerate(modes)}
    A = np.zeros((len(modes), len(modes)), dtype=np.int64)
    for z, x, spin, f, c in modes:
        if z == species and spin == 0:
            A[index[z, x, 0, f, c], index[z, x, 1, f, c]] = 1
    R = integral(one_body(A, basis))
    z = np.array([sum((1 - 2 * s) for i, (a, x, s, f, c) in enumerate(modes)
                      if a == species and mask >> i & 1) for mask in basis], dtype=np.int64)
    return (sp.diags(z * z + 2 * z, dtype=np.int64) + 4 * R.T @ R).tocsr()


@lru_cache(None)
def components(kind):
    d = native.components(kind)
    modes, basis = d['modes'], d['basis']
    nq = np.array([sum(1 for i, m in enumerate(modes) if m[0] == 0 and mask >> i & 1)
                   for mask in basis], dtype=np.int64)
    na = np.array([sum(1 for i, m in enumerate(modes) if m[0] == 1 and mask >> i & 1)
                   for mask in basis], dtype=np.int64)
    S4 = angular4(modes, basis, (0, 1), True)
    pair = (sp.diags(6 * nq * na, dtype=np.int64) + S4
            - species_spin4(modes, basis, 0) - species_spin4(modes, basis, 1)).tocsr()
    pair.eliminate_zeros()
    if kind == 'baryon':
        assert pair.nnz == 0
    else:
        assert (pair - S4).nnz == 0
        assert (pair @ pair - 8 * pair).nnz == 0
    reduced = (d['Binteger'].T @ pair @ d['Binteger']).toarray()
    hpair = reduced / np.sqrt(np.outer(d['Gram'], d['Gram']))
    assert np.max(abs(pair @ d['B'] - d['B'] @ hpair)) < 1e-12
    return {**d, 'parts': {**d['parts'], 'Pair4': pair},
            'integer_reduced_parts': {**d['integer_reduced_parts'], 'Pair4': reduced},
            'hparts': {**d['hparts'], 'Pair4': hpair}}


def matrix(kind, p=None):
    p = parameters() if p is None else p
    return sum(float(v) * components(kind)['hparts'][k] for k, v in weights(p).items())


@lru_cache(None)
def sector(kind):
    d, p = components(kind), parameters()
    H = sum(float(v) * d['parts'][k] for k, v in weights(p).items())
    h = matrix(kind, p)
    assert np.max(abs(H @ d['B'] - d['B'] @ h)) < 3e-10
    return {**d, 'H': H, 'h': h, 'parameters': p}


@lru_cache(None)
def exact_matrix(kind):
    d, p = components(kind), parameters()
    A = sy.zeros(len(d['Gram']))
    for k, v in weights(p).items():
        A += sy.Rational(v.numerator, v.denominator) * sy.Matrix(d['integer_reduced_parts'][k].tolist())
    g = d['Gram']
    return sy.Matrix(len(g), len(g), lambda i, j: A[i, j] / sy.sqrt(int(g[i]) * int(g[j])))
