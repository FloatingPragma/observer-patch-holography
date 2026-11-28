"""Compare conditional interference and inclusive transport on the same exterior links."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import expm

from . import isolated_medium as lab

HERE = lab.HERE
READS = lab.READS


def graph(board, mask):
    near = board.adjacency().tolil()
    for u, v in mask:
        near[u, v] = near[v, u] = 0
    return near.tocsr()


def direct_generator(board, mask, rho, hb, stiffness):
    h, i, j = lab.lab.ordinary.two_loads(board, stiffness, mask)
    entries = h.tocoo()
    g = np.zeros((board.n, board.n))
    a, b = entries.row, entries.col
    np.add.at(g, (i[a], i[b]), entries.data * rho[j[a], j[b]])
    free = np.diag(hb @ rho)
    g[np.diag_indices(board.n)] -= free.sum() - free
    return g / np.sqrt(np.outer(1 - rho.diagonal(), 1 - rho.diagonal()))


def contact_swap():
    law = lab.lab.cold.emission.LAW
    h = lab.lab.cold.emission.four_h(lab.lab.K)
    swap = np.zeros((9, 9))
    for i, (a, b) in enumerate(law.SRC):
        swap[law.SRC.index((b, a)), i] = 1
    contrasts = np.array([(a - b) ** 2 for a, b in law.SRC])
    np.testing.assert_array_equal(
        h[:9, :9], 2 * (np.eye(9) - swap) + 4 * lab.lab.K * np.diag(contrasts)
    )
    # Blank words remain in the native sector and the contrast phases rephase.
    np.testing.assert_array_equal(h[9:, :9], 0)
    np.testing.assert_allclose(expm(-1j * np.pi * h[:9, :9] / 4), swap, atol=4e-14)
    return swap


def sampled_arm(board, mask, shells, rng, stiffness=lab.lab.K):
    hb = lab.lab.cold.motion(board, mask, 1, stiffness)[0]
    factor = lab.medium.density_factors(hb)[0]["equilibrium"]
    rho = factor @ factor.T
    p = rho.diagonal()
    adj = graph(board, mask)
    src = 4 * stiffness
    incidence = board.degree - np.asarray(adj.sum(axis=1)).ravel()
    s = adj @ p
    edges = np.array(board.seams)
    live = np.array([edge not in mask for edge in board.seams])
    u, v = edges[live].T
    cross = 1 - p[u] - p[v] + rho[u, v]
    J = np.zeros_like(rho)
    J[u, v] = J[v, u] = -cross / (2 * np.sqrt((1 - p[u]) * (1 - p[v])))
    coherent = np.zeros(board.n)
    np.add.at(coherent, u, rho[u, v])
    np.add.at(coherent, v, rho[u, v])
    clock = (
        (board.degree - incidence) / 2
        - src * incidence
        - (src + 0.5) * s / (1 - p)
        + coherent / (2 * (1 - p))
    )
    np.fill_diagonal(J, clock)
    np.testing.assert_allclose(
        J, direct_generator(board, mask, rho, hb, stiffness), atol=5e-14
    )
    launch = np.sqrt((1 - p) / (board.n - 1))
    hpsi = J @ launch
    np.testing.assert_allclose(
        hpsi / launch, -src * (incidence + s / (1 - p)), atol=8e-13
    )
    acc = 2 * (hpsi**2 - launch * (J @ hpsi))
    predicted = np.bincount(shells, weights=acc).cumsum()[1:-1]
    fringe = []
    for x, y in edges:
        q = p[x] + p[y]
        fringe.append(
            rng.multinomial(
                READS, [(q + 2 * rho[x, y]) / 2, (q - 2 * rho[x, y]) / 2, 1 - q]
            ).tolist()
        )
    swaps = []
    for cut in lab.medium.cuts(board, shells):
        x, y = sorted(zip(*cut, strict=True))[0]
        qx, qy = 1 - p[x], 1 - p[y]
        t = 1 - p[x] - p[y] + rho[x, y]
        # Controlled native swap on one successful injection arm, parked on the other.
        probs = [(qx + qy + 2 * t) / 4, (qx + qy - 2 * t) / 4, 1 - (qx + qy) / 2]
        swaps.append(
            dict(
                edge=[int(x), int(y)],
                interference=rng.multinomial(READS, probs).tolist(),
                transferred=int(rng.binomial(READS, qy)),
                wrong_port=0,
            )
        )
    return dict(
        sites=rng.multinomial(READS, p).tolist(), fringes=fringe, swaps=swaps
    ), predicted


def interval_product(x, y):
    terms = np.stack([a * b for a in x for b in y])
    return np.stack([terms.min(axis=0), terms.max(axis=0)])


def difference(a, b):
    return np.asarray(a) - np.asarray(b)[::-1]


def decode(data, tape):
    board = lab.lab.cold.federation(0)
    marks, _ = lab.lab.cold.read_map(board, tape["spectral_tape"])
    shells = lab.lab.cold.distances(board, marks[0])
    idx = {edge: k for k, edge in enumerate(board.seams)}
    gain = lab.ordinary_gain(tape)
    eps = np.sqrt(lab.lab.RISK / (2 * READS))
    result = {}
    for tag, mask in (("blank", []), ("written", marks)):
        arm = data["arms"][tag]
        counts = np.asarray(arm["sites"])
        assert min(counts) >= 0 and counts.sum() == READS
        p = counts / READS
        adj = graph(board, mask)
        s = adj @ p
        density = np.stack([np.maximum(0, p - eps), np.minimum(1, p + eps)])
        neighbors = np.stack([np.maximum(0, s - eps), np.minimum(1, s + eps)])
        conditional = neighbors / (1 - density)
        shots = np.asarray(arm["fringes"])
        assert shots.min() >= 0 and np.all(shots.sum(axis=1) == READS)
        coherence = (shots[:, 0] - shots[:, 1]) / (2 * READS)
        correlations = np.stack([coherence - eps, coherence + eps])
        result[tag] = []
        for swap, (u, v) in zip(
            arm["swaps"], lab.medium.cuts(board, shells), strict=True
        ):
            edges = [idx[tuple(sorted((x, y)))] for x, y in zip(u, v, strict=True)]
            div = np.bincount(u, minlength=board.n) - np.bincount(v, minlength=board.n)
            # Correlated residual retains its leading term as one observable.
            residual = div @ p + np.array([-1, 1]) * np.ptp(div) * eps
            for a, b, sign in ((u, v, 1), (v, u, -1)):
                product = interval_product(
                    difference(density[:, b], correlations[:, edges]), conditional[:, a]
                )
                residual += sign * product.sum(axis=1)[::sign]
            residual = interval_product(residual, gain) / (board.n - 1)
            bond = (
                1 - density[:, u][::-1] - density[:, v][::-1] + correlations[:, edges]
            )
            gradient = difference(conditional[:, u], conditional[:, v])
            coherent = interval_product(
                interval_product(bond, gradient).sum(axis=1), gain
            ) / (board.n - 1)
            x, y = swap["edge"]
            assert [x, y] == list(sorted(zip(u, v, strict=True))[0]) and swap[
                "wrong_port"
            ] == 0
            good, bad, fail = swap["interference"]
            assert good + bad + fail == READS and min(good, bad, fail) >= 0
            t = (good - bad) / READS
            vacuum = np.sqrt((1 - density[:, x]) * (1 - density[:, y]))
            # Higher p means a smaller denominator. The observed t is positive.
            vis = (t + np.array([-2, 2]) * eps) / vacuum
            arrival = swap["transferred"] / READS + np.array([-1, 1]) * eps
            assert arrival[1] >= 1 - density[1, y] and arrival[0] <= 1 - density[0, y]
            result[tag].append(
                dict(
                    residual=residual.tolist(),
                    coherent=coherent.tolist(),
                    conditional_swap_visibility=vis.tolist(),
                    wrong_port_upper=min(
                        1.0, -np.expm1(-lab.lab.RISK / max(1, swap["transferred"]))
                    ),
                )
            )
    diff = []
    for w, b in zip(result["written"], result["blank"], strict=True):
        diff.append(
            {
                k: difference(w[k], b[k]).tolist()
                for k in ("residual", "coherent", "conditional_swap_visibility")
            }
        )
    return dict(arms=result, contrast=diff)


def run():
    path = HERE / "mediated_release.json"
    tape = json.loads(path.read_bytes())
    board = lab.lab.cold.federation(0)
    marks, _ = lab.lab.cold.read_map(board, tape["spectral_tape"])
    shells = lab.lab.cold.distances(board, marks[0])
    rng = np.random.default_rng(480078)
    contact_swap()
    arms, pred = {}, []
    for name, mask in (("blank", []), ("written", marks)):
        arm, coherent = sampled_arm(board, mask, shells, rng)
        arms[name] = arm
        pred.append(coherent)
    out = dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                path,
                HERE / "undisturbed_medium.json",
                HERE / "POISSON_RESPONSE.md",
            )
        },
        trials=READS,
        arms=arms,
    )
    out["decoded"] = decode(out, tape)
    prior = json.loads((HERE / "undisturbed_medium.json").read_bytes())["arms"]
    native = (
        np.array(prior["written_equilibrium"]["acceleration"])[1:-1]
        - np.array(prior["blank_equilibrium"]["acceleration"])[1:-1]
    )
    for idx, row in enumerate(out["decoded"]["contrast"]):
        c = (pred[1] - pred[0])[idx]
        assert row["coherent"][0] <= c <= row["coherent"][1]
        assert row["residual"][0] <= native[idx] - c <= row["residual"][1]
    # Repeat the actual native release at zero contrast with the same temperature.
    # Its stationary acceleration and its finite-time shape are separate readings.
    changed = {}
    for label, mask in (("blank", []), ("written", marks)):
        raw, coherent_zero = sampled_arm(board, mask, shells, rng, stiffness=0)
        hb = lab.lab.cold.motion(board, mask, 1, 0)[0]
        thermal = lab.medium.density_factors(hb)[0]["equilibrium"]
        h, i, j = lab.lab.ordinary.two_loads(board, 0, mask)
        inclusive, laws = lab.medium.fly(h, i, j, shells, thermal)
        np.testing.assert_allclose(inclusive[1:-1], 0, atol=4e-13)
        np.testing.assert_allclose(coherent_zero, 0, atol=4e-13)
        changed[label] = dict(
            swaps=raw["swaps"],
            max_exterior_acceleration=float(max(abs(inclusive[1:-1]))),
            terminal_shell_counts=[
                rng.multinomial(READS, row / row.sum()).tolist() for row in laws
            ],
        )
    out["zero_stiffness"] = changed
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "exchange_geometry.json"
    if args.counts_only:
        out = json.loads(path.read_bytes())
        lab.lab.cold.apparatus.check(
            out["decoded"],
            decode(out, json.loads((HERE / "mediated_release.json").read_bytes())),
        )
    else:
        out = run()
        if args.write:
            path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
        else:
            lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(out["decoded"]["contrast"][:3])
