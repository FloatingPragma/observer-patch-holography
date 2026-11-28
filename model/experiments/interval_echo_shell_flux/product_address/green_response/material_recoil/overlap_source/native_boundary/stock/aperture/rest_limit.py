"""Rational spectral enclosures for the free rest preparation at all long mean ages."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np
from scipy.linalg import eigh

from . import isolated_medium as lab

HERE = lab.HERE
SCALE = 10**14
RADIUS = Q(1, 10**9)


def spectral_disks(h):
    """Enclose eigenvalue clusters of the integer generator A=2H."""
    a = (2 * h).astype(np.int64)
    assert (2 * h != a).nnz == 0
    e, v = eigh(a.toarray().astype(float))
    centers = [int(round(x * SCALE)) for x in e]
    m = np.rint(v * SCALE).astype(np.int64)
    n = len(e)
    gram = m.astype(object).T @ m.astype(object)
    square = SCALE**2
    gram_error = max(
        abs(gram - square * np.eye(n, dtype=np.int64).astype(object)).sum(axis=1)
    )
    assert Q(int(gram_error), square) < Q(3, 4)  # inverse norm below two
    assert max(np.asarray(abs(a).sum(axis=1)).ravel()) * SCALE < 2**62
    raw = (a @ m).astype(object) * SCALE - m.astype(object) * np.array(
        centers, dtype=object
    )
    error = max(max(abs(raw).sum(axis=i)) for i in (0, 1))
    assert 2 * Q(int(error), square) < RADIUS
    return centers, m, gram


def eigen_witness(h):
    centers, m, gram = spectral_disks(h)
    n = len(centers)
    bank = []
    for i, center in enumerate(centers):
        distance = min(
            Q(abs(center - other), SCALE) for j, other in enumerate(centers) if j != i
        )
        if distance <= 2 * RADIUS:
            continue
        # Eigenvectors outside this singleton are at least distance-RADIUS
        # from its candidate center. The normalized column residual is < RADIUS.
        eta = RADIUS / (distance - RADIUS)
        col = m[:, i]
        norm2 = int(gram[i, i])
        weight = Q(sum(map(int, col)) ** 2, n * norm2)
        bank.append((col, norm2, weight, eta, distance - 2 * RADIUS))
    return bank


def limits(obs, bank, n):
    # Retain four times the unnormalized exterior operator exactly.
    a = np.rint(4 * obs).astype(np.int64)
    np.testing.assert_array_equal(a, 4 * obs)
    norm = Q(int(max(abs(a).sum(axis=1))), 4)
    # The blank uniform wave is a zero mode of the blank Laplacian.
    blank = Q(int(a.sum()), 4 * n)
    granularity = 10**13
    center = error = missing = rate = Q(0)
    for col, norm2, p, eta, gap in bank:
        # Keep the integer products in range before the exact scalar reduction.
        assert max(abs(a).sum(axis=1)) * SCALE < 2**63 - 1
        moment = Q(
            sum(int(i) * int(j) for i, j in zip(col, a @ col, strict=True)), 4 * norm2
        )
        center += Q(math.floor(p * moment * granularity), granularity)
        error += Q(1 + math.ceil(norm * eta * (1 + 2 * p) * granularity), granularity)
        missing += Q(math.floor((p - eta) * granularity), granularity)
        rate += Q(math.ceil(8 * norm / gap), 1)
    error += norm * (1 - missing)
    low, high = [(center - blank + sign * error) / (n - 1) for sign in (-1, 1)]
    rate /= n - 1
    margin = max(0, low, -high)
    # Enclose the infinite-age interval and coefficient outwards to readable decimals.
    grid = 10**7
    return dict(
        limit_lower=math.floor(low * grid) / grid,
        limit_upper=math.ceil(high * grid) / grid,
        finite_window_error_numerator=math.ceil(rate),
        sign_guaranteed_from_age=math.floor(rate / margin) + 1 if margin else None,
    )


def run():
    source = HERE / "mediated_release.json"
    tape = json.loads(source.read_bytes())
    board = lab.lab.cold.federation(0)
    marks, _ = lab.lab.cold.read_map(board, tape["spectral_tape"])
    shells = lab.lab.cold.distances(board, marks[0])
    h, meters = lab.instruments(board, shells, marks)
    hb, blank = lab.instruments(board, shells, [])
    assert np.all(hb @ np.ones(board.n) == 0)
    for current, reference in zip(meters, blank, strict=True):
        np.testing.assert_array_equal(current[2], reference[2])
    bank = eigen_witness(h)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (source, HERE / "REST_LIMIT.md", HERE / "isolated_medium.json")
        },
        ports=board.n,
        singletons=len(bank),
        centers_radius=str(RADIUS),
        cuts=[limits(o, bank, board.n) for _, _, o in meters],
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "rest_limit.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), result)
    print(result["singletons"], result["cuts"])
