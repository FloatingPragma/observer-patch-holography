"""Modular spectra of regional equilibrium references on an interacting written board."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path

from model.executable_model import rotations
from .. import source
from ..free_probe import measure as flight
from ...closed_region.measure import check

HERE = Path(__file__).parent


def folded():
    board, seam, dist, _ = flight.make_board()
    h = flight.generator(board, seam, 1, 8, True)[0].toarray()
    a, b = seam
    perm = next(p for p in rotations() if 12 + p[a] == b)
    order = np.array(
        [
            4 * site + m
            for i in range(12)
            for site in (i, 12 + perm[i])
            for m in range(4)
        ]
    )
    hf = h[np.ix_(order, order)]
    J = source.LAW._write_pattern()
    tags = [i for i, (_, delta) in enumerate(source.LAW.LABELS) if abs(delta) == 1]
    ids = [
        9 * j + source.LAW.SRC.index(pair) for pair in ((1, 0), (0, 1)) for j in tags
    ]
    cols = J[np.ix_(ids, [source.LAW.SRC.index(w) for w in ((0, 1), (1, 0))])]
    p2 = cols @ cols.T
    refl = np.kron(np.eye(2, dtype=int), np.eye(4, dtype=int)[[1, 0, 3, 2]])
    identity = np.eye(8, dtype=int)
    channels = {}
    for e, twog in ((8, identity + refl), (9, identity - refl)):
        for image, p in ((1, p2), (0, 2 * identity - p2)):
            four = twog @ p
            assert np.array_equal(four.T, four) and np.array_equal(
                four @ four, 4 * four
            )
            channels[(e, image)] = four
    np.testing.assert_array_equal(sum(channels.values()), 4 * identity)
    root = np.zeros((12, 12))
    root[a, a] = 1
    internal = -board.adjacency().toarray()[:12, :12] / 2
    ordinary = internal + np.diag(
        8.5 * np.count_nonzero(internal, axis=1) + 8 * (6 - board.degree[:12])
    )

    def owner_h(e, image):
        return ordinary + e * (np.eye(12) - image * root)

    recombined = sum(np.kron(owner_h(*c), proj / 4) for c, proj in channels.items())
    np.testing.assert_array_equal(4 * recombined, 4 * hf)
    return board, seam, dist, h, order, channels, owner_h


def stocks():
    # 8 coordinates ordered by owner side then (c,d); d=-1 on this input.
    one = np.eye(8)
    first, second = one[4], one[2]
    P = (np.outer(first, first) + np.outer(second, second)) / 2
    F = P - (np.outer(first, second) + np.outer(second, first)) / 4
    E = np.diag([0.25, 0, 0.25, 0, 0.25, 0, 0.25, 0])
    return np.array([P, E, F])


def modular_spectra(thermal, sites, seam, states):
    # The local occupied block equals a principal compression when there is
    # exactly one global load. The empty block sums over external placements.
    ids = [4 * p + m for p in sites for m in range(4)]
    outside = np.setdiff1d(np.arange(24), sites)
    vacuum = sum(
        (thermal[4 * p : 4 * p + 4, 4 * p : 4 * p + 4] for p in outside),
        np.zeros((4, 4)),
    )
    compressed = thermal[np.ix_(ids, ids)]
    np.testing.assert_allclose(np.trace(vacuum) + np.trace(compressed), 1)
    vals, v = eigh(compressed)
    assert vals[0] > 0
    tags = [4 * sites.index(p) + m for p in seam for m in range(4)]
    v = v[tags]
    # Combine degenerate bins before comparing probability laws.
    groups = np.split(
        np.arange(len(vals)), np.flatnonzero(np.diff(-np.log(vals)) < -1e-9) + 1
    )
    probs = np.stack(
        [np.einsum("ki,nkl,li->n", v[:, ix], states, v[:, ix]) for ix in groups], axis=1
    )
    modular = np.array([-np.log(vals[ix]).mean() for ix in groups])
    means = probs @ modular
    # Finite resolution adds uniform independent read noise of full width w.
    # Integrate the piecewise constant difference of the two convolved laws.
    resolution = {}
    for width in (0.001, 0.01, 0.1):
        edges = np.r_[modular - width / 2, modular + width / 2]
        delta = (probs[1] - probs[2]) / width
        jumps = np.r_[delta, -delta]
        ix = np.argsort(edges)
        resolution[str(width)] = float(
            np.dot(np.diff(edges[ix]), abs(np.cumsum(jumps[ix])[:-1])) / 2
        )
    return dict(
        sites=sites,
        level_nats=modular.tolist(),
        probabilities=probs.tolist(),
        modular_means=means.tolist(),
        erased_increment=means[1] - means[0],
        feedback_increment=means[2] - means[0],
        uniform_resolution_total_variation=resolution,
        final_distribution_total_variation=float(abs(probs[1] - probs[2]).sum() / 2),
    )


def run():
    board, seam, dist, h, order, channels, owner_h = folded()
    states = stocks()
    weights = np.array(
        [np.einsum("ijk,kj->i", states, p / 4) for p in channels.values()]
    ).T
    np.testing.assert_array_equal(weights[1], weights[2])
    # Compressions on the entire signed-negative four-word source space.
    # Source order is (0,A),(0,B),(1,A),(1,B). Entries of 2B are integers.
    source_ids = [0, 4, 2, 6]
    native_branch = np.array([0, 1, 1, 0])
    twice_bright = np.outer(native_branch, native_branch)
    for (_, image), four_projector in channels.items():
        np.testing.assert_array_equal(
            four_projector[np.ix_(source_ids, source_ids)],
            twice_bright if image else 2 * np.eye(4) - twice_bright,
        )
    es, v = eigh(h)
    individual = shortest_path(
        board.adjacency(), indices=seam, directed=False, unweighted=True
    ).astype(int)
    banks = []
    for beta in np.log(8 / 7) * np.array([0.5, 1, 2]):
        prob = np.exp(-beta * (es - es[0]))
        rho = (v * (prob / prob.sum())) @ v.T
        cuts = {}
        for name, metric in (
            ("both", dist),
            ("first", individual[0]),
            ("second", individual[1]),
        ):
            for radius in range(max(metric) + 1):
                sites = np.flatnonzero(metric <= radius).tolist()
                if not set(seam) <= set(sites):
                    continue
                cuts[f"{name}_{radius}"] = modular_spectra(rho, sites, seam, states)
        # Independent channel calculation of an occupied reduced density and
        # its log from twelve-dimensional matrices instead of the 96 eigensystem.
        blocks = {}
        for key, p in channels.items():
            e, u = eigh(owner_h(*key))
            blocks[key] = (u * np.exp(-beta * (e - es[0]))) @ u.T
        z = sum(np.trace(channels[c] / 4) * np.trace(m) for c, m in blocks.items())
        np.testing.assert_allclose(
            sum(np.kron(m, channels[c] / 4) for c, m in blocks.items()) / z,
            rho[np.ix_(order, order)],
            atol=4e-15,
        )
        for r in range(max(dist) + 1):
            inside = np.flatnonzero(dist[:12] <= r)
            raw = []
            for c, block in blocks.items():
                mu, u = eigh(block[np.ix_(inside, inside)] / z)
                raw.append(np.dot(-np.log(mu), u[list(inside).index(seam[0])] ** 2))
            np.testing.assert_allclose(
                weights @ raw, cuts[f"both_{r}"]["modular_means"], atol=3e-13
            )
        banks.append(dict(beta=beta, cuts=cuts))
    return dict(
        stocks=["paired", "erased", "phase_feedback"],
        channels=[list(c) for c in channels],
        channel_weights=weights.tolist(),
        temperature_panels=banks,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    res = json.loads(json.dumps(run()))
    path = HERE / "results.json"
    if p.parse_args().write:
        path.write_text(json.dumps(res, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), res)
