"""Exact additive internal generators on the ternary port and eleven-state memory."""

import argparse
from functools import reduce
import json
from math import gcd, lcm
from pathlib import Path

import numpy as np
import sympy as sp

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate


def integer_terms(kappa):
    write = law._write_pattern()
    swap, cost = law._load_ops()
    reverse, magnitude = law._label_ops()
    twice_g = np.kron(
        np.eye(10, dtype=int) - reverse + 2 * kappa * magnitude, np.eye(9, dtype=int)
    )
    complement = 2 * np.eye(90, dtype=int) - write @ write.T
    assert np.array_equal(complement @ twice_g, twice_g @ complement)
    blank = 2 * (np.eye(9, dtype=int) - swap + 2 * kappa * cost)
    return blank, complement @ twice_g, write


def unit(n, i):
    a = np.zeros((n, n), dtype=np.int64)
    a.flat[i] = 1
    return a


def equations(kappa):
    blank, marked, write = integer_terms(kappa)
    free, echo = [], []
    for i in range(9 + 100):
        port = unit(3, i) if i < 9 else np.zeros((3, 3), dtype=int)
        head = unit(10, i - 9) if i >= 9 else np.zeros((10, 10), dtype=int)
        load = np.kron(port, np.eye(3, dtype=int)) + np.kron(np.eye(3, dtype=int), port)
        written = np.kron(np.eye(10, dtype=int), load) + np.kron(
            head, np.eye(9, dtype=int)
        )
        free.append(
            np.r_[
                (load @ blank - blank @ load).ravel(),
                (written @ marked - marked @ written).ravel(),
            ]
        )
        echo.append(
            np.r_[
                (written @ write - write @ load).ravel(),
                (write.T @ written - load @ write.T).ravel(),
            ]
        )

    def reduced(cols):
        mat = np.asarray(cols).T
        return np.unique(mat[np.any(mat != 0, axis=1)], axis=0)

    return reduced(free), reduced(echo)


def exact_basis(rows):
    mat = sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(rows)).to_sparse()
    raw = mat.nullspace().to_Matrix().tolist()
    out = []
    for vec in raw:
        den = reduce(lcm, (int(sp.denom(z)) for z in vec), 1)
        nums = [int(z * den) for z in vec]
        div = reduce(gcd, nums, 0)
        out.append(np.array([v // div for v in nums], dtype=np.int64))
    for v in out:
        assert not np.any(rows @ v)
    return out


def matrix_rank(rows):
    if not len(rows):
        return 0
    return (
        sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(np.stack(rows)))
        .to_sparse()
        .rank()
    )


def split(v):
    return v[:9].reshape(3, 3), v[9:].reshape(10, 10)


def characterize(rows):
    basis = exact_basis(rows)
    brackets = []
    for i, v in enumerate(basis):
        a, b = split(v)
        for w in basis[:i]:
            c, d = split(w)
            ab = np.r_[(a @ c - c @ a).ravel(), (b @ d - d @ b).ravel()]
            if ab.any():
                brackets.append(ab)
    neutral_labels = [i for i, (_, d) in enumerate(law.LABELS) if d == 0]
    outside = np.ones((10, 10), dtype=bool)
    outside[np.ix_(neutral_labels, neutral_labels)] = False
    return dict(
        dimension=len(basis),
        port_dimension=matrix_rank([v[:9] for v in basis]),
        derived_dimension=matrix_rank(brackets),
        brackets_confined_to_zero_contrast_coin=all(
            not v[:9].any() and not split(v)[1][outside].any() for v in brackets
        ),
    )


def projector(coins, contrasts):
    head = np.array([coins[c] * contrasts.get(d, 0) for c, d in law.LABELS])
    return np.outer(head, head)


def native_generators():
    basis = []
    correction = projector([1, 1], {2: 1}) + projector([1, 1], {-2: 1})
    for i in range(3):
        a = unit(3, 4 * i)
        # Twice the generator avoids halves in the symmetric coin projectors.
        b = (1, -2, 1)[i] * correction
        basis.append(np.r_[2 * a.ravel(), b.ravel()])
    for contrasts in ({0: 1}, {2: 1, -2: 1}, {2: 1, -2: -1}):
        b = projector([1, -1], contrasts)
        basis.append(np.r_[np.zeros(9, dtype=int), b.ravel()])
    return basis


def run():
    result = {}
    for k in (int(candidate.KAPPA), 0):
        free, echo = equations(k)
        both = np.vstack([free, echo])
        result[str(k)] = {
            "free": characterize(free),
            "free_and_echo": characterize(both),
        }
        if k == candidate.KAPPA:
            proposed = native_generators()
            assert all(not np.any(both @ v) for v in proposed)
            assert matrix_rank(proposed) == result[str(k)]["free_and_echo"]["dimension"]
            q = proposed[2] - proposed[0]
            a, b = split(q)
            assert np.array_equal(a, np.diag([-2, 0, 2])) and not b.any()
            write = law._write_pattern()
            # All purely memory generators above annihilate the populated write range.
            result[str(k)]["dark_write_norm_squared"] = [
                int(np.sum((np.kron(split(v)[1], np.eye(9, dtype=int)) @ write) ** 2))
                for v in proposed[3:]
            ]
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    result = run()
    with parser.parse_args().output.open("x") as f:
        f.write(json.dumps(result, indent=2) + "\n")
