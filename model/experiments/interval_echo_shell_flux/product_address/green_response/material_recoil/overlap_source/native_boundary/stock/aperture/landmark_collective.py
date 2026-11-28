"""Certify the unrestricted range response at the eight-source circuit launches."""

from fractions import Fraction as Q
import math

import numpy as np
from sympy import Matrix, Rational, ZZ, eye
from sympy.polys.matrices import DomainMatrix
from scipy.sparse.csgraph import laplacian

from . import closed_passage as passage
from . import collective_source as source


def eight_source_witness(witness):
    board = source.tidal.lab.energy.pair.stock.observations.federation(1)
    layout = source.packing(board, witness["mode"], 8)
    assert [list(p) for p, _ in layout] == witness["source_owner_paths"]
    ends = [s[-1] for _, s in layout]
    blank, slots = passage.contacts(board, [e for _, s in layout for e in s])
    L = DomainMatrix.from_list(laplacian(blank)[:-1, :-1].tolist(), ZZ).to_sparse()
    B = np.zeros((board.n, len(ends)), int)
    for j, (u, v) in enumerate(ends):
        B[u, j], B[v, j] = 1, -1
    rhs = DomainMatrix.from_list(B[:-1].tolist(), ZZ).to_sparse()
    X, den = L.solve_den(rhs)
    assert L * X == rhs * den
    Z = np.zeros(B.shape, object)
    Z[:-1] = X.to_list()
    divisor = math.gcd(den, *Z.ravel())
    Z, den = Z // divisor, den // divisor
    impedance = Matrix([list(Z[u] - Z[v]) for u, v in ends]) / den
    assert impedance == impedance.T
    prices = [source.chances(board.n, [name] * 8) for name in ("original", "source")]
    p, q = [Rational(row[0].numerator, row[0].denominator) for row in prices]
    assert all(len(set(row)) == 1 for row in prices) and p < q
    identity = eye(len(ends))
    update = (
        (q - p) * (identity + p * impedance).inv() * (identity + q * impedance).inv()
    )
    assert (
        update
        == q * (identity + q * impedance).inv() - p * (identity + p * impedance).inv()
    )
    fields = {}
    for x in [witness["start"], witness["probe"]]:
        star = sorted(set(slots[x]) - {x})
        assert all(x not in edge for _, seams in layout for edge in seams)
        assert list(sum((Z[j] - Z[x] for j in star))) == [0] * len(ends)
        grad = Matrix([list(Z[j] - Z[x]) for j in star]) / den
        response = -2 * (grad * update * grad.T).trace()
        assert response < 0
        fields[str(x)] = dict(
            raw_range_change=str(Q(int(response.p), int(response.q))),
            independent_bearing_change=0,
        )
    return dict(
        stations=board.n,
        contacts=ends,
        launches=fields,
        original_probability=str(p),
        repaired_probability=str(q),
    )
