"""Physical pair-coincidence screen on the 80-cell monitored sharp-word region."""

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from .. import measure as clock

HERE = Path(__file__).resolve().parent
FACES = tuple(tuple(j for j in range(4) if j != i) for i in range(4))
EDGES = clock.EDGES
PAIRS = tuple(tuple(EDGES.index(e) for e in combinations(face, 2)) for face in FACES)


def faces(edges):
    areas = []
    for pairs in PAIRS:
        a, b, c = edges[list(pairs)] ** 2
        areas.append(np.sqrt(4 * a * b - (a + b - c) ** 2) / 4)
    return np.asarray(areas)


def shape3(u):
    a, b, c = u.T
    return np.sqrt(3 * (a * b + b * c + c * a))


@dataclass
class Screens:
    occupied: np.ndarray
    terminal: np.ndarray
    support: np.ndarray
    rim_areas: np.ndarray
    clock_starts: np.ndarray
    clock_scores: np.ndarray

    def face_probs(self, p):
        # One heralded resident word and one shared monitoring trajectory.
        flux = (p[:, None] * self.occupied).T @ self.terminal
        h = p @ self.occupied
        hit = np.array([flux[i, list(idx)] for i, idx in enumerate(PAIRS)])
        # Two-copy sham: start each trajectory independently with the same apex
        # occupied. On each end draw one of the two residents uniformly, with
        # no tags. Multiply separate single-site rates for distinct face hits.
        singles = flux[:, 6:] / (2 * h[:, None])
        twin = np.array(
            [
                [2 * singles[i, a] * singles[i, b] for a, b in combinations(face, 2)]
                for i, face in enumerate(FACES)
            ]
        )
        return h, hit, twin

    def observables(self, p):
        h, hit, twin = self.face_probs(p)
        intensity = shape3(hit)
        shape = intensity / hit.sum(axis=1)
        sham = shape3(twin) / twin.sum(axis=1)
        coefficients = np.stack([shape, intensity, sham, h**2 * shape3(twin)])
        live = faces(clock.distances(p, self.clock_starts, self.clock_scores))
        return np.r_[coefficients @ self.rim_areas, coefficients @ live]


def instrument(s, c=clock.MEAN_TIMEOUT):
    z = 1 - 1 / c
    end = np.column_stack(
        [*(s.occupied[:, a] & s.occupied[:, b] for a, b in EDGES), *s.occupied.T]
    ).astype(float)
    gram = np.eye(len(s.words)) - z * s.transition
    # Read after at least one half-tick; G in {1,2,...} with mean c.
    term = cho_solve(
        cho_factor(gram, check_finite=False), s.transition @ end / c, check_finite=False
    )
    rim = faces(clock.distances(s.uniform, s.starts, s.scores))
    return Screens(s.occupied, term, end, rim, s.starts, s.scores)


def dyadic_stock(s, im):
    """Concrete ticket approximation to the rank-36 record-preserving direction."""
    start_edges = (*EDGES, (1, 0), (2, 0), (3, 0))
    en = s.starts[:, [clock.RAYS.index(e) for e in start_edges]]
    returns = s.starts - s.scores / 1024
    pairs = np.vstack(
        [(s.occupied[:, i, None] * im.terminal[:, e]).T for i, e in enumerate(PAIRS)]
    )
    rows = np.vstack(
        [
            np.ones(len(s.uniform)),
            s.h.diagonal() / 256,
            en.T,
            returns.T,
            256 * pairs,
            s.cut_count == 2,
        ]
    )
    rhs = np.zeros(len(rows))
    rhs[-1] = 1
    d = np.linalg.lstsq(rows, rhs, rcond=1e-13)[0]
    tickets = np.rint((s.uniform - d / 65536) * 2**50).astype(np.int64)
    tickets[tickets.argmax()] += 2**50 - tickets.sum()
    assert min(tickets) > 0
    return tickets


def interior_reallocation(s):
    """Redirect a fair coin's share between equal-energy pairs that miss every apex."""
    groups = defaultdict(list)
    for i, word in enumerate(s.words):
        if not s.occupied[i].any():
            energy_twice = int(2 * s.h[i, i])
            assert energy_twice == 2 * s.h[i, i]
            groups[energy_twice].append(i)
    # Each of the 276 original sharp words has two equiprobable coin outcomes.
    tickets = np.full(len(s.words), 2)
    for group in groups.values():
        paired = len(group) - len(group) % 2
        tickets[group[:paired:2]] -= 1
        tickets[group[1:paired:2]] += 1
    return dict(
        sorted((energy, len(group)) for energy, group in groups.items())
    ), tickets


def series(s, im):
    u = s.uniform
    grad, _ = clock.differentials(u, s.starts, s.scores)
    moves, ctrl, ctrl_words = clock.survey_directions(s, grad)
    directions = {
        name: d - ctrl * (d @ grad[:, 1]) / (ctrl @ grad[:, 1])
        for name, d in moves.items()
    }
    ticket_data = json.loads((clock.HERE / "counts.json").read_text())["held_lottery"]
    w = np.array(ticket_data["source_word_tickets"]) / ticket_data["total_tickets"]
    directions["held_clock_lottery"] = (u - w) * 1024
    base = im.observables(u)
    slopes = {
        k: im.observables(u + 1e-20j * d).imag / 1e-20 for k, d in directions.items()
    }
    eta = -(directions["zero_to_one"] @ s.entropy_scores) / slopes["zero_to_one"]
    report = {}
    for key, d in directions.items():
        ds = float(d @ s.entropy_scores)
        report[key] = dict(
            image_slopes=slopes[key].tolist(),
            entropy_slope=ds,
            calibrated_residual=(ds + eta * slopes[key]).tolist(),
            energy_slope=float(d @ s.h.diagonal()),
            survey_volume_slope=float(d @ grad[:, 1]),
        )
    h, hit, twin = im.face_probs(u)
    return dict(
        clock=clock.MEAN_TIMEOUT,
        stations=s.rim,
        baseline_launch=h.tolist(),
        absolute_pair_probabilities=hit.tolist(),
        conditional_sham_pair_probabilities=twin.tolist(),
        physical_face_areas=im.rim_areas.tolist(),
        metrics=[
            f"{clock_mode}_{metric}"
            for clock_mode in ("frozen", "live")
            for metric in (
                "conditional_pair_shape",
                "all_start_pair_image",
                "conditional_independent_shape",
                "all_start_independent_image",
            )
        ],
        baseline_images=base.tolist(),
        ordinary_gain=eta.tolist(),
        volume_control=ctrl_words,
        movements=report,
        held_clock_lottery_finite_image_differences=(im.observables(w) - base).tolist(),
    )


def run():
    s = clock.survey(1)
    return series(s, instrument(s))


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
