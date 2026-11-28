"""Infer the initial unmonitored material tide from finite sharp-contact reads."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json

import numpy as np
from scipy.linalg import expm
from scipy.stats import beta

from . import tidal_release as tidal

flight = tidal.flight
lab = tidal.lab
HERE = tidal.HERE


def spectral_contact():
    h = lab.energy.pair.stock.kernels()[2]
    ind = np.arange(99)
    material = ind[:, None] % 9 != ind[None, :] % 9
    h2 = (h**2 * material).sum(0)
    full = lab.gates()[0]
    p = (full * material).sum(0)
    # Probability p/4; the ordinary transfer coefficient is 1/4.
    # Disagreement in two independent copies is 2(p/4)(1-p/4).
    np.testing.assert_array_equal(h2, p + 72 * p * (4 - p))
    assert set(p) == {0, 2, 4}
    eig = (0, 8, 9, 32, 33)
    poly = {}
    identity = np.eye(99, dtype=np.int64)
    for e in eig:
        num = identity.copy()
        den = 1
        for f in eig:
            if f != e:
                num = num @ (h - 4 * f * identity)
                den *= 4 * (e - f)
        np.testing.assert_array_equal(h @ num, 4 * e * num)
        poly[e] = num.astype(object) * Q(1, den)
    assert np.array_equal(sum(poly.values()), identity)
    # A real symmetric local law has a finite cosine spectrum. Taylor's
    # fourth-order remainder bounds the finite-dwell curvature estimate.
    coeff = defaultdict(lambda: np.full(99, Q(0), object))
    for e, f in combinations(eig, 2):
        coeff[f - e] += 2 * (poly[e] * poly[f] * material).sum(0)
    curv = -sum(w * k * k for k, w in coeff.items()) / 2
    np.testing.assert_array_equal(16 * curv, h2)
    dwell = Q(1, 128)
    error = sum(abs(w) * k**4 for k, w in coeff.items()) * dwell**2 / 6
    held = (abs(expm(-1j * float(dwell) * h / 4)) ** 2 * material).sum(0)
    rng = np.random.default_rng(514799)
    N = 2**34
    counts = rng.binomial(N, held.clip(0, 1))
    # Independent period-pi shots classify the three possible native rates.
    # The count of discordant pairs is measured separately.
    fulls = rng.binomial(64, p / 4)
    pairs = rng.binomial(64, 2 * (p / 4) * (1 - p / 4))
    classes = np.where(fulls == 0, 0, np.where(fulls == 64, 4, 2))
    assert np.array_equal(classes, p)
    assert np.all(pairs[p == 2] > 0) and np.all(pairs[p != 2] == 0)
    eps = 0.01 / (2 * 99)
    lo = np.zeros(99)
    hi = np.ones(99)
    pos = counts > 0
    lo[pos] = beta.ppf(eps, counts[pos], N - counts[pos] + 1)
    pos = counts < N
    hi[pos] = beta.ppf(1 - eps, counts[pos] + 1, N - counts[pos])
    band = np.array([lo, hi]) * 4 / float(dwell) ** 2
    band[0] -= np.asarray(error, float)
    band[1] += np.asarray(error, float)
    exact = h2 / 4
    assert np.all((band[0] <= exact) & (band[1] >= exact))
    assert np.all(
        abs(4 * held / float(dwell) ** 2 - exact) <= np.asarray(error, float) + 1e-12
    )
    return h2, dict(
        full_probability_quarters=p.tolist(),
        twice_tested_pairs=64,
        full_swaps=fulls.tolist(),
        discordant_pairs=pairs.tolist(),
        classification_error_upper=float(2 * 99 / 2**64),
        ordinary_clock_coefficient="1/4",
        dwell=str(dwell),
        held_shots_per_column=N,
        held_swaps=counts.tolist(),
        finite_dwell_error=list(map(str, error)),
        held_normalized_coefficient_99=band.T.tolist(),
    )


@lru_cache(maxsize=None)
def expected_rate(pop, labels, t, launch, cats, h2tuple):
    value = Q(0)
    for word, p in tidal.draws(pop, cats.count(2)):
        bg = iter(word)
        ls = [next(bg) if c == 2 else launch[c] for c in cats]
        if len(ls) == 3 and ls[2] != t:
            continue
        a, b = ls[:2]
        if (a == t) == (b == t):
            continue
        val = sum(h2tuple[9 * m + 3 * a + b] for m in labels)
        # h²/4 is the rate divided by the measured blank single rate.
        value += p * Q(val, 4 * len(labels)) * (int(a == t) - int(b == t))
    return value


def closure(dist, t, board, marks, h2):
    A, B = marks[-1]
    h2 = tuple(map(int, h2))
    lanc = ((0, 0), (t, 0), (0, t), (t, t))
    result = []
    for i, j in board.seams:
        if {i, j}.isdisjoint((A, B)):
            continue

        def kind(x):
            return 0 if x == A else 1 if x == B else 2

        cats = (kind(i), kind(j))
        pops = defaultdict(Q)
        for (pop, heads), p in dist.items():
            pops[pop, heads[marks.index((i, j))] if (i, j) in marks else (0,)] += p

        def rate(context):
            return [
                sum(
                    p * expected_rate(pop, labels, t, x, context, h2)
                    for (pop, labels), p in pops.items()
                )
                for x in lanc
            ]

        j0 = rate(cats)
        by_kind = {}
        for k in (0, 1, 2):
            if k < 2 and k in cats:
                continue
            f = rate((*cats, k))
            by_kind[k] = f[3] - f[2] - f[1] + f[0]
        result.append((int(i), int(j), j0[2] - j0[1], by_kind))
    return result


def observe(h, ends, cost, edges):
    A, B = ends
    base, use = tidal.chart(h, ends, cost)
    # Only rows meeting an input and their incident neighbors are needed.
    sites = {A, B} | {v for i, j, _, _ in edges for v in (i, j)}
    d = {v: base[v] + cost * use[v] for v in sites}
    length = d[A][B]
    axis = (d[A] ** 2 - d[B] ** 2) / (2 * length**2)
    proj = Q(0)
    joint = Q(0)
    for i, j, s, co in edges:
        proj += s * (axis[j] - axis[i])
        diff = d[j] - d[i]
        other = diff.sum() - sum(diff[k] for k in {A, B, i, j})
        joint += co[2] * other
        for k, index in [(0, A), (1, B)]:
            if k in co:
                joint += co[k] * diff[index]
    return dict(
        initial_distance=str(length),
        normalized_chord_acceleration=str(2 * proj),
        normalized_pair_acceleration=str(2 * joint / length),
    )


def run():
    h2, meter = spectral_contact()
    old = json.loads((HERE / "tidal_release.json").read_bytes())["schedules"]["forward"]
    _, words, cells, *_ = lab.energy.pair.setup()
    board, _, marks, *_ = lab.energy.returned(words)
    _, h, _ = flight.records.thermal_radar.chart(board, marks)
    u = np.full(len(words), Q(1, len(words)), object)
    z = lab.energy.least_change(cells) / len(words)
    w = json.loads((HERE / "material_release.json").read_bytes())["matched_mixtures"]
    arms = dict(
        original=flight.records.ensemble(words, u),
        source=flight.records.ensemble(words, z),
        single=flight.records.ensemble(words, u, single=1),
        ordinary=flight.records.ensemble(words, u, extra=1),
    )
    arms.update(positive=flight.mix(arms, w[0]), reference=flight.mix(arms, w[1]))
    for name, e, s in [
        ("reference_plus_single", 0, 1),
        ("reference_plus_double", 1, 0),
    ]:
        arms[name] = flight.mix(
            dict(
                source=flight.records.ensemble(words, z, extra=e, single=s),
                ordinary=flight.records.ensemble(words, u, extra=e + 1, single=s),
            ),
            w[1],
        )
    out = {}
    for mode in ("all", "clean"):
        panel = {}
        for t in (1, 2):
            result = {}
            for name, dist in arms.items():
                _, arm = tidal.conditional(dist, t) if mode == "clean" else (Q(1), dist)
                if not arm:
                    continue
                edges = closure(arm, t, board, marks, h2)
                entry = old[mode][str(t)][name]
                costs = entry["chart_cost"]
                surveys = {
                    k: observe(h, marks[-1], Q(c), edges) for k, c in costs.items()
                }
                if name in ("source", "positive"):
                    ref = "original" if name == "source" else "reference"
                    rc = old[mode][str(t)][ref]["chart_cost"]["fastest_relay"]
                    surveys["fixed_reference_fastest"] = observe(
                        h, marks[-1], Q(rc), edges
                    )
                result[name] = dict(acceptance=entry["acceptance"], charts=surveys)
            panel[str(t)] = result
        out[mode] = panel
    # A modified law retains native diagonal and head-only amplitudes; it
    # forbids material exchange on the third contact. Evaluate at a common
    # blank reference geometry since that contact can no longer deliver loads.
    for t in (1, 2):
        vals = []
        for tag in ("original", "source"):
            edges = [
                e for e in closure(arms[tag], t, board, marks, h2) if e[:2] != marks[-1]
            ]
            vals.append(observe(h, marks[-1], Q(int(h[marks[-1]])), edges))
        assert vals[0] == vals[1]
    return dict(
        protocol_sha256=sha256((HERE / "COHERENT_TIDE.md").read_bytes()).hexdigest(),
        monitored_sha256=sha256((HERE / "tidal_release.json").read_bytes()).hexdigest(),
        contact=meter,
        panels=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    raw = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    path = HERE / "coherent_tide.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    print("Checked all 99 local columns and both relative-release instruments.")
