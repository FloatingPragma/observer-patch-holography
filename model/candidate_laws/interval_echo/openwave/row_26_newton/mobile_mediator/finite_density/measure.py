"""Two clamped neutral sources in a thermal medium with all vacancy numbers."""

from collections import deque
from functools import cache
from fractions import Fraction
from itertools import combinations, product
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp

from model.executable_model import federation
from .. import measure as mediator

HERE = Path(__file__).parent
DISTANCES = (1, 2, 3, 5, 6)
BETA = 1 / 16
MU = 30
SHOTS = 2**22


@cache
def apparatus():
    board = federation(0)
    adj = [set(map(int, row)) for row in board.neighbours]
    paths = {0: (0,)}
    todo = deque([0])
    while todo:
        u = todo.popleft()
        for v in sorted(adj[u]):
            if v not in paths:
                paths[v] = (*paths[u], v)
                todo.append(v)
    rail = max(paths.values(), key=len)[1:10]
    ports = set(rail)
    assert len(ports) == 9
    assert all(
        adj[x] & ports == set(rail[max(0, j - 1) : j] + rail[j + 1 : j + 2])
        for j, x in enumerate(rail)
    )
    pins = {}
    for j in (0, *DISTANCES):
        x = rail[j]
        pins[j] = min(
            (a, b)
            for a in adj[x] - ports
            if adj[a] & ports == {x} and board.degree[a] == 5
            for b in adj[a] - ports
            if not adj[b] & ports
            and board.degree[b] == 5
            and (j == 0 or not ({a, b} | adj[a] | adj[b]) & set(pins[0]))
        )
    return rail, pins, board


@cache
def sectors(kappa, cut_after):
    rail, _, board = apparatus()
    n_sites = len(rail)
    law = mediator.candidate.native_law
    seam = law.seam_term(kappa)
    hop = seam[law.SRC.index((0, 1)), law.SRC.index((1, 0))]
    interaction = mediator.contact(kappa)[0, 0] / 2
    single = (kappa + 0.5) * board.degree[list(rail)] - kappa * (
        6 - board.degree[list(rail)]
    )
    result = []
    for number in range(n_sites + 1):
        words = tuple(
            sum(1 << j for j in sites) for sites in combinations(range(n_sites), number)
        )
        address = {word: j for j, word in enumerate(words)}
        bits = (np.asarray(words)[:, None] >> np.arange(n_sites)) & 1
        h = np.diag(bits @ (single - MU)).astype(float)
        for j in range(n_sites - 1):
            if j == cut_after:
                # Remove the whole native seam, including both single-site diagonal increments.
                h[np.diag_indices_from(h)] -= (kappa + 0.5) * (
                    bits[:, j] + bits[:, j + 1]
                )
                continue
            h[np.diag_indices_from(h)] += interaction * bits[:, j] * bits[:, j + 1]
            for u, word in enumerate(words):
                if bits[u, j] != bits[u, j + 1]:
                    h[u, address[word ^ (3 << j)]] = hop
        result.append((h, bits))
    return result


def readings(kappa, d, a, b, cut_after=None, canonical=None):
    log_parts, log_occupied, right, sector_counts = [], [], [], []
    for n, (h, bits) in enumerate(sectors(kappa, cut_after)):
        if canonical is not None and n != canonical:
            continue
        energies, states = eigh(
            h + np.diag(a * bits[:, 0] + b * bits[:, d]), driver="evd"
        )
        exponent = -BETA * energies
        part = logsumexp(exponent)
        weights = np.exp(exponent - part)
        prob_words = (states * states) @ weights
        log_parts.append(part)
        sector_counts.append(n)
        right.append(float(prob_words @ bits[:, d]))
        if n == 0:
            log_occupied.append(part)
    log_z = float(logsumexp(log_parts))
    probabilities = np.exp(np.array(log_parts) - log_z)
    return dict(
        free_energy=-log_z / BETA,
        right_vacancy=float(probabilities @ right),
        empty=float(np.exp(log_occupied[0] - log_z)) if log_occupied else 0.0,
        mean_vacancies=float(probabilities @ sector_counts),
    )


def panel(kappa, d, cut_after=None, canonical=None):
    v = mediator.contact(kappa)[0] / 2
    bas = readings(kappa, d, 0, 0, cut_after, canonical)
    # The spectator-facing end has both orientations, with the same sharp static energy.
    arms = {
        (a, b): readings(kappa, d, a, b, cut_after, canonical)
        for a, b in product((0, v[0], v[2]), repeat=2)
    }
    signed = []
    for a, b in product((v[0], v[2]), repeat=2):
        f = (
            arms[a, b]["free_energy"]
            - arms[a, 0]["free_energy"]
            - arms[0, b]["free_energy"]
            + bas["free_energy"]
        )
        signed.append(float(f))
    coef = Fraction(-9, 16) * (2 * int(kappa) + 1) ** d * Fraction(1, 4) ** (d + 1)
    out = dict(
        distance=d,
        source_pairs=[list(apparatus()[1][j]) for j in (0, d)],
        signed_interactions=signed,
        quenched_mean=float(np.mean(signed)),
        unperturbed_mean_vacancies=bas["mean_vacancies"],
        right_vacancy=[arms[a, 0]["right_vacancy"] for a in (0, v[0], v[2])],
        high_temperature_coefficient=str(coef)
        if cut_after is None and canonical is None
        else None,
    )
    if canonical is None:
        out["empty"] = [
            [arms[a, b]["empty"] for b in (0, v[0], v[2])] for a in (0, v[0], v[2])
        ]
    return out


def acquire(row, seed):
    rng = np.random.default_rng(seed)
    probs = np.array(row["empty"])
    counts = rng.binomial(SHOTS, probs)
    # Simultaneous Hoeffding bounds for the nine empty-bath meters and three target meters.
    delta = np.sqrt(np.log(2 * 12 / 0.01) / (2 * SHOTS))
    right_counts = rng.binomial(SHOTS, row["right_vacancy"])
    if (
        np.max(abs(counts / SHOTS - probs)) > delta
        or np.max(abs(right_counts / SHOTS - row["right_vacancy"])) > delta
    ):
        raise ValueError("finite thermal count disagrees with its forecast")
    return dict(
        seed=seed,
        per_setting=SHOTS,
        empty_counts=counts.tolist(),
        right_vacancy=right_counts.tolist(),
        simultaneous_99_half_width=float(delta),
    )


def run():
    rail, pins, _ = apparatus()
    out = []
    for kappa in (mediator.candidate.KAPPA, 4):
        rows = [panel(kappa, d) for d in DISTANCES]
        out.append(
            dict(
                kappa=kappa,
                panels=rows,
                acquisition=acquire(rows[0], 915538 + int(kappa)),
            )
        )
    cut_after = 1
    return dict(
        rail=list(rail),
        pair_order=[
            "vacancy-vacancy",
            "vacancy-excess",
            "excess-vacancy",
            "excess-excess",
        ],
        beta=BETA,
        mu=MU,
        settings=out,
        cut_grand=panel(8, 3, cut_after),
        cut_four_vacancies=panel(8, 3, cut_after, 4),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        (HERE / "results.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
    else:
        saved = json.loads((HERE / "results.json").read_text())

        def check(actual, old):
            if isinstance(actual, dict):
                assert actual.keys() == old.keys()
                for key in actual:
                    check(actual[key], old[key])
            elif isinstance(actual, list):
                for a, b in zip(actual, old, strict=True):
                    check(a, b)
            elif isinstance(actual, float):
                np.testing.assert_allclose(actual, old, rtol=0, atol=1e-12)
            else:
                assert actual == old

        check(result, saved)
        print("Replayed both complete finite-density spectra and both count banks.")
