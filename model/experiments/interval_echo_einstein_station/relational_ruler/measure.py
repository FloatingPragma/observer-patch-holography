"""Joint native address survey with the second tag's bearing suppressed."""

import argparse
from collections import deque
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.stats import beta

from model.executable_model import cell
from model.executable_model.refinement import load_refinement
from model.candidate_laws.interval_echo.candidate import KAPPA

HERE = Path(__file__).resolve().parent
SHOTS = 2**32
T = 1 / 64
TERMS = 9


class Routes:
    def __init__(self, car):
        self.internal = cell().neighbours
        self.external = dict(zip(car.external[:, 0], car.external[:, 1], strict=True))
        self.external.update(zip(car.external[:, 1], car.external[:, 0], strict=True))

    @lru_cache(maxsize=70000)
    def neighbours(self, site):
        base = site // 12 * 12
        result = tuple(base + n for n in self.internal[site - base])
        if site in self.external:
            return result + (int(self.external[site]),)
        return result

    def range_ball(self, start, radius):
        seen = {start: 0}
        pending = deque([start])
        while pending:
            site = pending.popleft()
            dist = seen[site]
            if dist < radius:
                for dest in self.neighbours(site):
                    if dest not in seen:
                        seen[dest] = dist + 1
                        pending.append(dest)
        return seen


def census(routes, home, radius):
    allsites = np.array(sorted(k for k, v in home.items() if v <= 2 * radius))
    ids = {int(k): i for i, k in enumerate(allsites)}
    points = np.array([ids[k] for k, v in home.items() if v <= radius])
    irow, icol = [], []
    for i, site in enumerate(allsites):
        for dest in routes.neighbours(int(site)):
            if dest in ids:
                irow.append(i)
                icol.append(ids[dest])
    graph = csr_matrix(
        (np.ones(len(irow)), (irow, icol)), shape=(len(allsites), len(allsites))
    )
    multiplicity = 0
    occupancies = []
    for start in range(0, len(points), 64):
        dist = dijkstra(
            graph,
            directed=False,
            indices=points[start : start + 64],
            unweighted=True,
            limit=radius,
        )
        # A path of ≤radius from an enrolled anchor cannot exit the 2*radius home ball.
        pop = np.count_nonzero(np.isfinite(dist), axis=1) - 1
        assert np.all(
            np.any(dist == radius, axis=1)
        )  # Every intermediate integer occurs.
        occupancies.extend(pop.tolist())
        multiplicity += int(pop.sum())
    return dict(
        radius=radius,
        anchors=len(points),
        labels=len(points) * radius,
        ordered_pairs=multiplicity,
        minimum_pairs_per_anchor=min(occupancies),
        maximum_pairs_per_anchor=max(occupancies),
    )


def bearings(routes, home, radius):
    classes = {}
    for y, dist in home.items():
        if dist != radius:
            continue
        ns = routes.neighbours(y)
        inward = sum(home.get(s) == radius - 1 for s in ns)
        classes.setdefault(len(ns), []).append((inward, y))
    candidates = []
    for degree, items in classes.items():
        worst, lo = min(items)
        best = max(n for n, _ in items)
        hi = min(y for n, y in items if n == best)
        candidates.append((best - worst, -degree, lo, hi, worst, best))
    change, negdegree, lo, hi, worst, best = max(candidates)
    return -negdegree, [(lo, worst), (hi, best)]


def sharp_evolution(routes, anchor, vacancy):
    # Native blank seam reduction on the full registered filled background.
    # H-14 here is the shifted ordered-pair Hamiltonian with radius bound 30.
    coefficients = {(anchor, vacancy): 1.0}
    amplitudes = coefficients.copy()
    prefactor = 1
    for order in range(1, TERMS + 1):
        next_coeff = {}
        for (x, y), v in coefficients.items():
            nx, ny = routes.neighbours(x), routes.neighbours(y)
            adjacent = int(y in nx)
            diagonal = (len(nx) + len(ny) - adjacent) / 2 + 2 * KAPPA * (
                len(ny) - len(nx) + adjacent
            )
            next_coeff[x, y] = next_coeff.get((x, y), 0.0) + (diagonal - 14) * v
            for new in nx:
                pair = (y, x) if new == y else (new, y)
                next_coeff[pair] = next_coeff.get(pair, 0.0) - 0.5 * v
            for new in ny:
                if new != x:
                    pair = (x, new)
                    next_coeff[pair] = next_coeff.get(pair, 0.0) - 0.5 * v
        coefficients = next_coeff
        prefactor *= (-1j * T) / order
        for pair, v in coefficients.items():
            amplitudes[pair] = amplitudes.get(pair, 0.0) + prefactor * v
    epsilon = (30 * T) ** (TERMS + 1) / math.factorial(TERMS + 1)
    # Taylor's integral remainder uses the norm-one unitary along the real time axis.
    return amplitudes, epsilon


def survey(routes, home, radius, level):
    degree, choices = bearings(routes, home, radius)
    results = []
    for i, (y, inward) in enumerate(choices):
        amplitude, eps = sharp_evolution(routes, 0, y)
        probabilities = np.zeros(4)
        local_ranges = routes.range_ball(0, radius + TERMS)
        for (x, dest), a in amplitude.items():
            k = (
                3
                if x != 0
                else (
                    0
                    if local_ranges[dest] < radius
                    else 2
                    if local_ranges[dest] > radius
                    else 1
                )
            )
            probabilities[k] += abs(a) ** 2
        norm_squared = probabilities.sum()
        assert abs(norm_squared - 1) <= 2 * eps + eps**2 + 1e-13
        p = probabilities / norm_squared
        draws = np.random.default_rng([229191, level, i]).multinomial(SHOTS, p)
        tail = 0.01 / (2 * 8)
        hits = draws[0]
        interval = (
            beta.ppf(tail, hits, SHOTS - hits + 1),
            beta.ppf(1 - tail, hits + 1, SHOTS - hits),
        )
        # Density error after normalization is ≤4 eps + O(eps²) for probabilities.
        results.append(
            dict(
                vacancy=y,
                degree=degree,
                closer_neighbours=inward,
                quadratic_inward_probability_coefficient=inward / 4,
                probabilities=p.tolist(),
                counts=draws.tolist(),
                closer_interval_99=list(interval),
                norm_error_bound=eps,
                probability_error_bound=4 * eps + 4 * eps**2 + 1e-13,
            )
        )
    return results


def run(directory):
    prior = json.loads((HERE.parent / "proper_ball_ruler" / "coarse.json").read_text())
    results = []
    for panel in prior["panels"]:
        level = panel["level"]
        car = load_refinement(directory / f"FEDERATION_LEVEL{level}.npz", level=level)
        routes = Routes(car)
        radius = panel["microscopic_ball_radius"]
        home = routes.range_ball(0, max(4 * radius, radius + TERMS))
        counts = [census(routes, home, r) for r in (radius, 2 * radius)]
        assert [c["anchors"] for c in counts] == [
            panel["ball_count"],
            panel["doubled_ball_count"],
        ]
        dimensions = {
            key: math.log2(counts[1][key] / counts[0][key])
            for key in ("labels", "ordered_pairs")
        }
        experiment = survey(routes, home, radius, level)
        results.append(
            dict(
                level=level,
                cells=car.cells,
                host_sha256=car.sha256,
                radius=radius,
                census=counts,
                dimensions=dimensions,
                bearings=experiment,
            )
        )
        routes.neighbours.cache_clear()
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        native_exposure=T,
        shots_per_bearing=SHOTS,
        panels=results,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = run(args.directory)
    with args.output.open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for r in data["panels"]:
        print(
            r["cells"],
            r["dimensions"],
            [
                (v["degree"], v["closer_neighbours"], v["counts"][0] / SHOTS)
                for v in r["bearings"]
            ],
            flush=True,
        )
