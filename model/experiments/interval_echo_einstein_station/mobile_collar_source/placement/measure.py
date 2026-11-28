"""Census of the original sharp two-head source at fixed physical reader orientations."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy import sparse

from .. import measure as source

HERE = Path(__file__).resolve().parent
TIMES = (1, 2, 4)
ORDER = 256
SCALE = 32
SEED = 1827096


def routes(board, all_neighbours):
    for a, b in sorted([*board.external, *((y, x) for x, y in board.external)]):
        left = [p for p in board.neighbours[a] if p // 12 == a // 12]
        right = [p for p in board.neighbours[b] if p // 12 == b // 12]
        assert len(left) == len(right) == 5
        for home, remote in product(
            left if all_neighbours else left[:1], right if all_neighbours else right[:1]
        ):
            yield (home, a, b, remote)


def initial_index(route, charge):
    a, b = source.starting_coins(route, charge)
    return 2 * a + b


def sign_class(i):
    return "alike" if i in (0, 3) else "unlike"


def geometry_counts(board):
    counts = Counter()
    for route in routes(board, True):
        home, a, b, remote = route
        # Any walk with a distant flip in its first three hops must
        # reach an end of the remote armed seam in its first two.
        for target in (b, remote):
            assert target not in board.neighbours[home]
            common = set(board.neighbours[home]) & set(board.neighbours[target])
            assert common == ({a} if target == b else set())
        assert source.starting_coins(route, -1) != source.starting_coins(route, 1)
        counts[sign_class(initial_index(route, -1))] += 1
    return dict(counts)


def native_bounds():
    result = {}
    for charge in (-1, 1):
        blocks = source.local_parts(charge, 8)
        for forward in (False, True):
            col = int((not forward) ^ (charge == 1))
            old, new = (1, 2) if forward else (2, 1)
            hop = blocks[new, old][:, col]
            np.testing.assert_array_equal(hop, np.eye(2)[:, 1 - col] * (-4))
        # An armed site has four/five blank neighbours and one live head.
        # Its diagonal correction is 0 or -4, with at most one off-diagonal 4.
        bounds = []
        for d, armed in product((5, 6), (0, 1)):
            diag = 8.5 * (d - armed) + 8 * (2 * charge + charge**2) * (6 - d) - 40
            for corr in (-4, 0) if armed else (0,):
                bounds.append(abs(diag + corr) + 0.5 * (d - armed) + 4 * armed)
        result[str(charge)] = max(bounds)
    assert max(result.values()) <= SCALE
    return result


@lru_cache(maxsize=None)
def coefficients(time):
    # Jacobi-Anger weights. Each Bessel power series at integer z=16t
    # is truncated after K with an exact rational partial sum and remainder.
    z, k = SCALE * time // 2, 224
    weights = []
    for n in range(ORDER + 1):
        top, bottom = 1, 1
        for j in range(k, 0, -1):
            d = j * (n + j)
            top, bottom = bottom * d - z * z * top, bottom * d
        exact = Fraction(z**n * top, factorial(n) * bottom)
        tail = Fraction(z ** (n + 2 * k + 2), factorial(k + 1) * factorial(n + k + 1))
        assert tail < Fraction(1, 10**23)
        weight = float(exact) * (1 if n == 0 else 2)
        assert abs(
            Fraction.from_float(weight) - exact * (1 if n == 0 else 2)
        ) < Fraction(1, 2**51)
        weights.append(weight)
    tail = (
        2
        * Fraction(z ** (ORDER + 1), factorial(ORDER + 1))
        / (1 - Fraction(z, ORDER + 2))
    )
    assert tail < Fraction(1, 10**20)
    return np.array(weights)


def scaled_h(board, route, charge):
    h = source.hamiltonian(board, route, charge, 8, center=False)
    h = (h - 40 * sparse.eye(h.shape[0], format="csr")) * (1 / SCALE)
    assert (h != h.T).nnz == 0 and abs(h).sum(axis=1).max() <= 1
    assert np.diff(h.indptr).max() <= 7
    return h


def propagate_scaled(h, route, init):
    v = np.zeros(h.shape[0])
    v[4 * route[0] + init] = 1
    current = v
    weights = np.stack([coefficients(t) for t in TIMES])
    real = weights[:, 0, None] * v
    imag = np.zeros_like(real)
    for n in range(1, ORDER + 1):
        nxt = h @ current
        if n > 1:
            nxt = 2 * nxt - v
        v, current = current, nxt
        s = 1 if n % 4 in (0, 3) else -1
        accumulator = imag if n % 2 else real
        accumulator += s * weights[:, n, None] * current
    rows = (real**2 + imag**2).reshape(len(TIMES), -1, 4).sum(axis=1)
    np.testing.assert_allclose(rows.sum(axis=1), 1, atol=2e-11, rtol=0)
    return rows


def propagate(board, route, charge):
    return propagate_scaled(
        scaled_h(board, route, charge), route, initial_index(route, charge)
    )


def covariance(p):
    return p[..., 3] - (p[..., 2] + p[..., 3]) * (p[..., 1] + p[..., 3])


def statistics(records, radius, sample_n, population):
    by_time = {}
    for i, t in enumerate(TIMES):

        def group(rows):
            all_p = np.array([row["p"][i] for row in rows])
            cov = covariance(all_p)
            return dict(
                n=len(rows),
                positive=int(sum(cov > radius)),
                negative=int(sum(cov < -radius)),
                unresolved=int(sum(abs(cov) <= radius)),
                mean=float(cov.mean()),
                pooled=float(covariance(all_p.mean(axis=0))),
                minimum=float(min(cov)),
                maximum=float(max(cov)),
            )

        groups = {
            key: group([row for row in records if row["class"] == key])
            for key in ("alike", "unlike")
        }
        groups["all"] = group(records)
        # Simultaneous six sampled means (two charges and three clocks) on
        # the held board. Serfling for individual binary covariances in [-1/4,1/4].
        halfwidth = (
            0.0
            if population == sample_n
            else float(
                0.5
                * np.sqrt(
                    (1 - (sample_n - 1) / population)
                    * np.log(2 * 6 / 0.01)
                    / (2 * sample_n)
                )
            )
        )
        groups["all"]["placement_mean_radius_99"] = halfwidth + radius
        by_time[str(t)] = groups
    return by_time


def uncertainty():
    # IEEE binary64 round-to-nearest. Each symmetric sparse recurrence uses
    # <=7 real terms per row, ||abs(A)||2 <= 1 and exact binary scaling.
    # If its inputs have norm <=2, one local recurrence error is < 2e-14.
    # Propagation by U_j(A) amplifies that by no more than j+1.
    recur = ORDER * (ORDER + 1) / 2 * 2e-14
    assert recur < 1
    weight_norm = max(sum(abs(coefficients(t))) for t in TIMES)
    # Coefficient and summation errors for the two real accumulators.
    state = (
        weight_norm * recur + 2 * (ORDER + 1) * (2**-50 + 2**-50 * weight_norm) + 1e-20
    )
    # Squaring and summing <=3840 positions (both real channels); true
    # individual, marginal and joint projector probabilities obey this bound.
    probability = state * (2 + state) + 1e-11
    return dict(
        state=state,
        projection=probability,
        covariance=3 * probability + probability**2,
        coefficient_l1=float(weight_norm),
        scaled_matrix_norm_upper=1,
    )


def controls():
    # With zero contrast all written hopping blocks and coin rotations vanish.
    # Uniform unit material at k=8 keeps both written heads scalar for all time.
    for charge in (-1, 1):
        np.testing.assert_array_equal(source.local_parts(charge, 0), 0)
        np.testing.assert_array_equal(
            source.local_parts(charge, 8)[0, 0], 8 * np.eye(2)
        )
    return dict(changed_contrast=0, uniform_load=0, product_of_measured_marginals=0)


def acquire():
    uncertainty_record = uncertainty()
    data = {}
    for level in range(3):
        b = source.federation(level)
        all_routes = list(routes(b, False))
        n = len(all_routes)
        selected = (
            list(range(n))
            if level < 2
            else sorted(
                map(
                    int,
                    np.random.Generator(np.random.PCG64(SEED)).choice(
                        n, 96, replace=False
                    ),
                )
            )
        )
        home = all_routes.index(source.region(b))
        rows = {str(q): [] for q in (-1, 1)}
        for index in sorted(set([*selected, home])):
            route = all_routes[index]
            for q in (-1, 1):
                coin = initial_index(route, q)
                rows[str(q)].append(
                    dict(
                        index=index,
                        route=route,
                        initial=coin,
                        **{"class": sign_class(coin)},
                        p=propagate(b, route, q).tolist(),
                    )
                )
        report = {
            q: statistics(
                [r for r in rr if r["index"] in selected],
                uncertainty_record["covariance"],
                len(selected),
                n,
            )
            for q, rr in rows.items()
        }
        data[str(b.cells)] = dict(
            all_arm_placements=geometry_counts(b),
            all_minimum_arm_placements=dict(
                Counter(sign_class(initial_index(route, -1)) for route in all_routes)
            ),
            directed_bridges=n,
            sample=selected,
            home=home,
            summary=report,
            records=rows,
        )
    return dict(
        protocol_sha256=sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        times=TIMES,
        row_bounds=native_bounds(),
        error=uncertainty_record,
        controls=controls(),
        boards=data,
    )


def compare(old, new):
    if isinstance(old, dict):
        assert old.keys() == new.keys()
        for k in old:
            compare(old[k], new[k])
    elif isinstance(old, list):
        assert len(old) == len(new)
        for a, b in zip(old, new, strict=True):
            compare(a, b)
    elif isinstance(old, float):
        assert abs(old - new) < 2e-11
    else:
        assert old == new


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = acquire()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n")
    else:
        compare(json.loads(path.read_text()), data)
    for size, board in data["boards"].items():
        print(
            size,
            board["all_arm_placements"],
            {q: row["4"]["all"] for q, row in board["summary"].items()},
        )
