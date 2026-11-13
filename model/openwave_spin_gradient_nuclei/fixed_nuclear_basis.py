"""Exact native projector columns for the complete nuclear proof coordinates.

The indices only choose coordinates on the existing Casimir/parity range.
Every replay rebuilds that range and proves that the selected columns span it;
no floating rank-revealing factorization chooses a proof basis.
"""
from functools import lru_cache
import json
from math import gcd
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import sympy as sy
from model.openwave_nuclear_distance import nuclear
from model.openwave_six_quark import source as native

WITNESSES = json.loads((Path(__file__).with_name('spectral_witnesses.json')).read_text())


@lru_cache(None)
def projector(S, I, parity):
    assert type(S) is int and S in range(4)
    assert type(I) is int and I in range(4)
    assert type(parity) is int and parity in (-1, 1)
    d = nuclear.carrier()
    eye = sp.eye(len(d['keep']), dtype=np.int64, format='csr')
    P, den = d['P'].copy(), d['den']
    for op, j in ((d['Sg'], S), (d['Ig'], I)):
        for eigenvalue in (0, 2, 6, 12):
            if eigenvalue == j * (j + 1):
                continue
            P = native.exact_product(op - eigenvalue * eye, P)
            den *= j * (j + 1) - eigenvalue
            common = gcd(abs(den), *map(int, P.data))
            P.data //= common
            den //= common
    P = native.exact_product(eye + parity * d['parity'], P)
    den *= 2
    common = gcd(abs(den), *map(int, P.data))
    P.data //= common
    den //= common
    if den < 0:
        P, den = -P, -den
    assert den > 0 and (P - P.T).nnz == 0
    assert (native.exact_product(P, P) - den * P).nnz == 0
    trace = sum(map(int, P.diagonal()))
    assert trace % den == 0
    rank = trace // den
    assert 0 <= rank <= P.shape[0]
    assert (native.exact_product(d['Sg'], P) - S * (S + 1) * P).nnz == 0
    assert (native.exact_product(d['Ig'], P) - I * (I + 1) * P).nnz == 0
    assert (native.exact_product(d['parity'], P) - parity * P).nnz == 0
    assert (native.exact_product(d['P'], P) - d['den'] * P).nnz == 0
    return P, den, rank


def channel_from_record(S, I, parity, record):
    P, den, rank = projector(S, I, parity)
    if not rank:
        assert record is None
        return {'rank': 0}
    assert isinstance(record, dict)
    indices = record['native_projector_column_indices']
    assert type(indices) is list and len(indices) == rank
    assert all(type(i) is int and 0 <= i < P.shape[1] for i in indices)
    assert len(set(indices)) == rank
    B = P[:, indices].tocsr()
    assert (native.exact_product(P, B) - den * B).nnz == 0
    Gram = sy.Matrix(native.exact_product(B.T, B).toarray().tolist())
    _, D = Gram.LDLdecomposition(hermitian=True)
    assert all(x > 0 for x in D.diagonal())
    # rank(P/den) is its exact trace. This many independent range columns
    # therefore span every vector in the native spin/isospin/parity channel.
    return {'rank': rank, 'projector_denominator': den, 'Binteger': B,
            'Gram': Gram, 'spin': S, 'isospin': I, 'parity': parity}


@lru_cache(None)
def channel(S, I, parity):
    return channel_from_record(S, I, parity, WITNESSES.get(f'nuclear_{S}_{I}_{parity}'))
