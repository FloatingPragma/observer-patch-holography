"""Ordinary-calibrated single and double insertion costs on the closed source."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np

from . import closed_energy as energy

HERE = energy.HERE
N = 2**38
# At most twelve new means, each taking values in [0,33].
WIDTH = 33 * np.sqrt(np.log(2400) / (2 * N))


def local_spectra():
    h = energy.pair.stock.kernels()[2][:9, :9]
    spectra = {}
    for n in range(3):
        for t in range(3):
            if n == t:
                assert np.all(h[:, 3 * n + t] == 0)
                spectra[n, t] = [0, 0]
            else:
                ids = [3 * n + t, 3 * t + n]
                diag = 2 + 32 * (n - t) ** 2
                np.testing.assert_array_equal(
                    h[np.ix_(ids, ids)], [[diag, -2], [-2, diag]]
                )
                spectra[n, t] = [(diag - 2) // 4, (diag + 2) // 4]
    return spectra


def changes(w, cell_weights, near):
    board, home, marks, _, _, _ = energy.returned(w)
    assert near not in home and all(near not in m for m in marks)
    h = energy.pair.stock.kernels()[2]
    f = np.zeros((len(w), board.n), int)
    f[:, home] = w[:, :4]
    shift = (cell_weights - 1) / len(w)
    out = {}
    for t in (1, 2):
        total = np.full(len(w), 32 * t**2 * (6 - board.degree[near]), int)
        for v in board.neighbours[near]:
            # All contacts incident to near are blank, including the one
            # with the middle load. The local column is symmetric in port order.
            m = f[:, v]
            total += h.diagonal()[3 * m + t] - h.diagonal()[3 * m]
        out[t] = shift @ total / 4
    return out


def run():
    _, w, cells, _, _, _, _, _ = energy.pair.setup()
    z = energy.least_change(cells)
    n = w[:, 1]
    spectrum = local_spectra()
    calibrated = json.loads((HERE / "closed_energy.json").read_bytes())
    price = changes(w, z, calibrated["protected_echo"]["sender"])
    assert price == {1: Q(3, 56), 2: 0}
    distributions = {
        "native": np.bincount(n) / len(n),
        "repair": np.bincount(n, weights=np.asarray(z, float)) / len(n),
        "empty_center": [1, 0, 0],
        "unit_center": [0, 1, 0],
    }
    rng = np.random.default_rng(671328)
    levels = np.array([0, 8, 9, 32, 33])
    raw = {}
    means = {}
    for tag, dist in distributions.items():
        raw[tag] = {}
        for t in range(3):
            p = np.zeros(34)
            for m, prob in enumerate(dist):
                for value in spectrum[m, t]:
                    p[value] += prob / 2
            bank = rng.multinomial(N, p[levels] / p.sum())
            means[tag, t] = float(bank @ levels / N)
            raw[tag][str(t)] = bank.tolist()

    def compare(a, b, t):
        m = (means[b, t] - means[b, 0]) - (means[a, t] - means[a, 0])
        return np.array([m - 4 * WIDTH, m + 4 * WIDTH])

    gains = [compare("empty_center", "unit_center", t) for t in (1, 2)]
    delta = [compare("native", "repair", t) for t in (1, 2)]
    assert gains[0][0] < -17 < gains[0][1] and gains[1][0] < -32.5 < gains[1][1]
    assert all(x[0] < price[t] < x[1] for t, x in enumerate(delta, 1))
    scalar = [x * y / v for x in delta[0] for y in gains[1] for v in gains[0]]
    assert min(scalar) > delta[1][1]
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "CLOSED_INSERTION.md", HERE / "closed_energy.json")
        },
        trials_per_bank=N,
        levels=levels.tolist(),
        counts=raw,
        ordinary_calibration=[x.tolist() for x in gains],
        from_closed_source={str(t): str(v) for t, v in price.items()},
        measured_increments=[x.tolist() for x in delta],
        common_scalar_double_prediction=[min(scalar), max(scalar)],
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    out = run()
    path = HERE / "closed_insertion.json"
    if opt.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == out
    print(
        out["ordinary_calibration"],
        out["measured_increments"],
        out["common_scalar_double_prediction"],
    )
