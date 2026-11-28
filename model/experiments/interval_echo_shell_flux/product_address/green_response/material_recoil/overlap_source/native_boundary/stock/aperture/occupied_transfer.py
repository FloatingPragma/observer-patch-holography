"""Predict held terminal motion from separate ordinary short-range interference reads."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import occupied_release as flight

HERE = flight.HERE
T = 64
SHOTS = 2**72
BATCH = 2**52


def flag_bank(rng, prob):
    # Return exact integer sums; do not sum an entire arm in a 64-bit counter.
    total = np.zeros(5, dtype=object)
    for _ in range(SHOTS // (2048 * BATCH)):
        draws = rng.multinomial(BATCH, prob / prob.sum(), size=2048)
        total += draws.sum(0, dtype=np.uint64).astype(object)
    assert sum(total) == SHOTS
    return total.tolist()


def flags(meter, rho):
    i, j = np.nonzero(np.triu(meter, 1))
    diag = meter.diagonal()
    degree = abs(meter).sum(axis=1)
    norm = max(degree)
    d = rho.diagonal()
    single = np.array([np.maximum(diag, 0) @ d, np.maximum(-diag, 0) @ d])
    weights = abs(meter[i, j])
    bare = weights @ (d[i] + d[j]) / 2
    signal = meter[i, j] @ rho[i, j]
    counters = np.r_[single, bare + signal, bare - signal] / norm
    assert min(counters) > -1e-13
    np.testing.assert_allclose(
        (counters @ [1, -1, 1, -1]) * norm, np.sum(rho * meter), atol=5e-14
    )
    return float(norm), np.r_[np.maximum(counters, 0), 1 - sum(counters)]


def acquire(twice, marks, rad, r):
    inside = rad <= r
    meter, source = flight.material_ops(twice / 2, marks, inside)
    assert not source.any()
    probability = []
    n = len(rad)
    for a, b in flight.SECTORS:
        h = twice / 2 + np.diag(a * marks[0] + b * marks[1])
        e, v = np.linalg.eigh(h)
        alpha = v.T @ (np.ones(n) / np.sqrt(n))
        diff = e[:, None] - e
        # Real part at |U-V| for independent uniform clocks U,V.
        cos_avg = np.sinc(diff * T / (2 * np.pi)) ** 2
        rho = (v * alpha) @ cos_avg @ (v * alpha).T
        assert abs(np.trace(rho) - 1) < 2e-13
        norm, p = flags(meter, rho)
        probability.append(p)
        # Check the finite spectral identity at every sharp-source calibration.
        wave = v @ (alpha * np.exp(-1j * T * e))
        terminal = abs(wave[inside]) @ abs(wave[inside])
        np.testing.assert_allclose(
            (terminal - inside.mean()) * 2 / T**2,
            np.sum(meter * rho),
            rtol=0,
            atol=2e-13,
        )
    return norm, flight.CONDITIONAL @ np.array(probability)


def run():
    board, old, edges, rad, _ = flight.chart()
    p, q = flight.sector_prior(old)
    delta = flight.lab.probe.log_factor(old["source"])
    release = json.loads((HERE / "occupied_release.json").read_bytes())["results"]
    rng = np.random.default_rng(775920)
    eps = np.sqrt(flight.RISK / (2 * SHOTS))
    result = {}
    for tag, load in (("vacancy", 0), ("excess", 2)):
        twice, marks = flight.apparatus(board, edges, load, -1, 1)
        result[tag] = {}
        for r in (1, 2):
            norm, sharp = acquire(twice, marks, rad, r)
            banks = [
                flag_bank(rng, v) for v in np.vstack([sharp, p @ sharp, q @ sharp])
            ]
            freq = np.array(banks, dtype=float) / SHOTS
            parts = {}
            for label, sign in (
                ("position", [1, -1, 0, 0, 0]),
                ("interference", [0, 0, 1, -1, 0]),
                ("total", [1, -1, 1, -1, 0]),
            ):
                readings = (freq @ sign) * norm * T**2 / 2
                # The signed flag is in [-1,1]; independent native and held
                # differences and the four calibration means use one union bound.
                half = 2 * norm * eps * T**2 / 2
                expected = flight.lab.probe.mix_bands(
                    delta, flight.MIX @ readings[:4] + np.array([-4, 4]) * half
                )
                held = readings[5] - readings[4] + np.array([-2, 2]) * half
                parts[label] = dict(
                    from_source=expected.tolist(), from_material=held.tolist()
                )
            terminal = release[tag]["panels"]["free"]["releases"][str(T)][str(r)][
                "held"
            ]
            assert all(
                max(terminal[0], v[0]) < min(terminal[1], v[1])
                for v in parts["total"].values()
            )
            result[tag][str(r)] = dict(
                norm=norm, counts=banks, terminal=terminal, components=parts
            )
    return dict(
        sha256={
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "OCCUPIED_TRANSFER.md",
                "occupied_release.json",
                "pair_contact.json",
            )
        },
        shots_per_arm=SHOTS,
        time=T,
        flags=[
            "positive_site",
            "negative_site",
            "positive_fringe",
            "negative_fringe",
            "idle",
        ],
        preparations=["00", "01", "10", "11", "native", "product"],
        results=result,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    out = run()
    path = HERE / "occupied_transfer.json"
    if opt.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    for tag, rows in out["results"].items():
        for r, row in rows.items():
            print(tag, r, row["components"]["total"], row["terminal"])
