"""Two one-stop physical record routes form a finite deposit on the guide chart."""

import argparse
import json
from itertools import combinations

import numpy as np
from scipy.fft import idstn

from . import patch, source
from .. import measure as apparatus
from ..closed_region.measure import check

HERE = source.HERE
ALPHAS = (1 / 64, 1 / 512, 1 / 4096)


def triangle(n, level):
    _, external = patch.board(level)
    length = 2**level
    face = (external // 12) // (length * length)
    sides = {}
    for (a, b), (left, right) in zip(external, face, strict=True):
        if (left == 0) == (right == 0):
            continue
        inside, outside = (a, right) if left == 0 else (b, left)
        sides.setdefault(int(outside), []).append(int(inside))
    assert len(sides) == 3 and all(len(side) == length for side in sides.values())
    origins = []
    for axes, side in zip(combinations(range(3), 2), sorted(sides), strict=True):
        for rank, _ in enumerate(sorted(sides[side])):
            point = np.full(3, n // 2 - 1 - round(2 * length / 3), int)
            point[axes[0]] += 2 * rank + 1
            point[axes[1]] += 2 * (length - rank) - 1
            assert min(point) >= 0 and max(point) < n
            origins.append(tuple(int(x) for x in point))
    assert len(set(origins)) == 3 * length
    return origins


def propagate(live, steps, rng):
    """Step all active cohorts once, retaining the boundary loss."""
    following = np.zeros_like(live)
    for j, step in enumerate(steps[:-2]):
        take = rng.binomial(live, 1 / (8 - j))
        live -= take
        axis = np.flatnonzero(step)[0]
        a, b = np.moveaxis(following, axis, 0), np.moveaxis(take, axis, 0)
        if step[axis] > 0:
            a[1:] += b[:-1]
        else:
            a[:-1] += b[1:]
    return following + live


def sample_stops(origin, n, alpha, steps, rng, num):
    live = np.zeros((n,) * 3, dtype=np.int64)
    live[origin] = num
    totals = np.zeros_like(live)
    for age in range(4096):
        if not live.any():
            break
        stopped = rng.binomial(live, alpha)
        totals += stopped
        live -= stopped
        live = propagate(live, steps, rng)
    assert not live.any()
    return np.r_[totals.ravel(), num - totals.sum()]


def spectral(n, q):
    modes = np.arange(1, n + 1) * np.pi / (n + 1)
    axis = 2 * q * (1 - np.cos(modes))
    return modes, axis[:, None, None] + axis[None, :, None] + axis[None, None, :]


def green(origin, n, modes, eigen, alpha):
    basis = [np.sqrt(2 / (n + 1)) * np.sin((x + 1) * modes) for x in origin]
    freq = basis[0][:, None, None] * basis[1][None, :, None] * basis[2][None, None, :]
    return idstn(freq / (alpha + (1 - alpha) * eigen), type=1, norm="ortho")


def whole_space_bound(alpha, q):
    """Fourier bound for a common stopping address of independent 3D couriers."""
    return min(1.0, np.pi**2 / 64 * (alpha / ((1 - alpha) * q)) ** 1.5)


def displacement_controls(n, q, modes, eigen):
    center = n // 2 - 1
    records = []
    for label, offsets in (("same origin", (0, 0)), ("separated origins", (-8, 8))):
        origins = [(center + x, center, center) for x in offsets]
        alpha = ALPHAS[1]
        stock = sum((alpha * green(s, n, modes, eigen, alpha)) ** 2 for s in origins)
        full = float(stock.sum())
        records.append(
            dict(
                placement=label,
                origins=origins,
                stop_probability=alpha,
                total_collision_per_write=full / 2,
                fraction_within_central_eight=stock[
                    center - 3 : center + 5,
                    center - 3 : center + 5,
                    center - 3 : center + 5,
                ].sum()
                / full,
            )
        )
    return records


def pair_sample(origins, expected, steps, alpha, gain):
    rng = np.random.default_rng(271128)
    n = 10
    trials = 2**26
    launches = np.zeros((n,) * 3, dtype=np.int64)
    rows = []
    for origin, deposit in zip(origins, expected, strict=True):
        left = sample_stops(origin, n, alpha, steps, rng, trials)
        right = sample_stops(origin, n, alpha, steps, rng, trials)
        original_right = right.copy()
        same = np.zeros(n**3, np.int64)
        for i, number in enumerate(left[:-1]):
            if not number:
                continue
            team = rng.multivariate_hypergeometric(right, number)
            same[i] = team[i]
            right -= team
        lost_both = int(right[-1])
        lost_one = int(left[-1] + original_right[-1] - 2 * lost_both)
        distinct = int(trials - lost_one - lost_both - same.sum())
        assert distinct >= 0
        p = float((deposit**2).sum())
        width = float(
            np.sqrt(2 * p * np.log(1200) / trials) + 2 * np.log(1200) / (3 * trials)
        )
        assert abs(same.sum() / trials - p) < width
        # Use the same energy ticket in both banks. A difference in the
        # nonzero energy branch occurs on 1/4 of joint destinations.
        extra_spectrum = rng.binomial(same, 1 / 4)
        odd = rng.binomial(extra_spectrum, 1 / 2)
        fired = rng.binomial(odd, 9 / 512) + rng.binomial(extra_spectrum - odd, 8 / 512)
        launches += fired.reshape((n,) * 3)
        rows.append(
            dict(
                origin=origin,
                coincident=int(same.sum()),
                distinct=distinct,
                one_lost=lost_one,
                both_lost=lost_both,
                expected_coincident=p * trials,
                probability_half_width=width,
                excess_tickets=int(fired.sum()),
            )
        )
    visits = np.zeros_like(launches)
    live = launches.copy()
    while live.any():
        visits += live
        live = propagate(live, steps, rng)
    forecasts = apparatus.green.poisson(
        sum(p * p for p in expected) * gain / 512, (1 / 8,) * 3
    )[1:-1, 1:-1, 1:-1]
    panels = []
    for probe in (1, 2):
        face = (6, 5, 5, 0)
        setting = apparatus.receiver.Setting("matched", face, probe, (2, 2, 2))
        arms = apparatus.receiver.count_panel(rng, visits, setting)
        score = sum(sum(a["right"]) - sum(a["starts"]) / 2 for a in arms)
        reward = apparatus.receiver.signed_readout(setting)
        panels.append(
            dict(
                probe=probe,
                arms=arms,
                signed_endpoints=score,
                expected_signed_endpoints=float((forecasts * reward).sum() * trials),
            )
        )
    return dict(
        trials_per_origin=trials,
        origins=rows,
        unmatched_launches=int(launches.sum()),
        measured_appointments=int(visits.sum()),
        panels=panels,
    )


def run():
    steps, conductance = apparatus.green.measured_steps(None)
    np.testing.assert_array_equal(conductance, [1 / 8] * 3)
    gain = float(source.run_case(8)[1]["sharp"]["gain"])
    out = []
    for n, levels in ((10, (0,)), (34, (0, 1, 2)), (82, (0, 1, 2, 3))):
        modes, eigen = spectral(n, conductance[0])
        diag = apparatus.diagonal_visits(n, conductance[0])
        for level in levels:
            points = triangle(n, level)
            caches = [green(point, n, modes, eigen, 0) / diag for point in points]
            first_pass = float(sum((p * p).sum() for p in caches))
            del caches
            rows = []
            for alpha in ALPHAS:
                deposits = [alpha * green(p, n, modes, eigen, alpha) for p in points]
                # An independent source on every written boundary edge. Collisions
                # are counted only between the two records of that specific write.
                stock = sum(p * p for p in deposits)
                assert min(p.min() for p in deposits) > -1e-15
                assert max(p.sum() for p in deposits) <= 1 + 1e-12
                collision = float(stock.sum()) / len(points)
                ceiling = whole_space_bound(alpha, conductance[0])
                assert max((p * p).sum() for p in deposits) <= ceiling
                field = apparatus.green.poisson(gain * stock / 512, conductance)
                ks = sorted({1, n // 2, *(k for k in (2, 4, 8, 16, 32) if k <= n // 2)})
                rows.append(
                    dict(
                        stop_probability=alpha,
                        deposited_per_route=float(np.mean([p.sum() for p in deposits])),
                        collision_per_write=collision,
                        whole_space_upper_bound=ceiling,
                        extra_energy_per_write=gain * collision,
                        shells=[
                            apparatus.shell(
                                k, gain * stock / 512, field, conductance[0]
                            )
                            for k in ks
                        ],
                    )
                )
                if n == 10 and alpha == ALPHAS[0]:
                    sampled = pair_sample(points, deposits, steps, alpha, gain)
            out.append(
                dict(
                    n=n,
                    patch_side=2**level,
                    origins=points,
                    separate_first_passage_queries_per_write=first_pass / len(points),
                    panels=rows,
                )
            )
    return dict(
        charts=out,
        short_guide_acquisition=sampled,
        equal_stock_displacement=displacement_controls(n, conductance[0], modes, eigen),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = json.loads(json.dumps(run()))
    path = HERE / "transport.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), result)
