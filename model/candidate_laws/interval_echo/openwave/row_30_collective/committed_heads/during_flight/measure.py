"""Native strand propagation with two resolved coin records and a distant wall read."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .. import measure as inherited
from . import readback

HERE = Path(__file__).parent
SHOTS = 2**24
SEED = 761322
TIMES = (2, 4, 6)
HANDS = (0, 1, 2, 4, 6)
LABELS = {
    "magnitude": ["2/2"],
    "once": ["0/2", "1/2"],
    "duplicate": ["0/0/2", "0/1/2", "1/0/2", "1/1/2"],
    "spaced": ["0/0", "0/1", "1/0", "1/1"],
}


def generators(kappa, n, meter, sensor):
    """Both coin parities in the conserved even-sign sector; common energy origin."""
    words = np.arange(2**n)
    bits = (words[:, None] >> np.arange(n)) & 1
    # Both material endpoints of a connecting head have load one on this word.
    aligned = (bits[:, :-1] == 1) & (bits[:, 1:] == 0)
    remote = bits[:, sensor - 1] == bits[:, sensor]
    rows = np.tile(words, n)
    cols = np.concatenate([words ^ (1 << i) for i in range(n)])
    exchange = sparse.csr_matrix(
        (np.full(rows.size, -0.5), (rows, cols)), shape=(len(words),) * 2
    )
    _, local_tables = inherited.readings()
    endpoints = (bits[:, meter - 1], 1 - bits[:, meter])
    even = local_tables[kappa, 1, 1][endpoints]
    others = 4 * kappa * (n - 1 - aligned.sum(axis=1)) - even
    # The two paths interfere at the second read. Different phase origins for
    # the two generators would amount to a different local coin experiment.
    origin = np.mean(others + even)
    return [
        (
            -1j
            * (
                exchange
                + sparse.diags(others + local_tables[kappa, 1, c][endpoints] - origin)
            )
        ).tocsr()
        for c in (1, -1)
    ], remote


def evolve(generator, vectors, duration):
    return expm_multiply(
        generator * duration, vectors, traceA=generator.diagonal().sum() * duration
    )


def panels(path, remote):
    """Unconditional joint laws, indexed by record, remote outcome and source."""
    ec, eo, oc, oo = path

    def squared(amplitude):
        return np.stack(
            [np.sum(abs(amplitude[mask]) ** 2, axis=0) for mask in (remote, ~remote)]
        )

    mag = squared(ec)
    one = (mag + squared(oo)) / 4
    # When two different coin labels were observed the odd evolution between
    # them has a minus sign; their unobserved final parities remain orthogonal.
    same = (squared(ec + eo) + squared(oc + oo)) / 16
    flip = (squared(ec - eo) + squared(oc - oo)) / 16
    zero = np.zeros_like(one)
    return {
        "magnitude": np.stack([mag]),
        "once": np.stack([one, one]),
        "duplicate": np.stack([one, zero, zero, one]),
        "spaced": np.stack([same, flip, flip, same]),
    }


def histories(kappa, n, meter, sensor, bonds, times=TIMES):
    (even, odd), remote = generators(kappa, n, meter, sensor)
    source = np.zeros((2**n, len(bonds)))
    source[[inherited.initial_word(n, b) for b in bonds], np.arange(len(bonds))] = 1
    one = evolve(even, source, 1)
    initial = np.column_stack((evolve(even, one, 1), evolve(odd, one, 1)))
    after_even = after_odd = initial
    previous = 2
    for t in times:
        if t > previous:
            after_even = evolve(even, after_even, t - previous)
            after_odd = evolve(odd, after_odd, t - previous)
        path = [*np.hsplit(after_even, 2), *np.hsplit(after_odd, 2)]
        rows = panels(path, remote)
        for probabilities in rows.values():
            np.testing.assert_allclose(
                probabilities.sum(axis=(0, 1)), 1, atol=3e-10, rtol=0
            )
        yield t, rows
        previous = t


def run():
    inherited.strand.native_checks()
    route, src, ref, meter, sensor = inherited.geometry()
    rng = np.random.default_rng(SEED)
    clocks = []
    for t in HANDS:
        probability = np.round((1 + np.array([np.cos(t), np.sin(t)])) / 2, 12)
        clocks.append(
            {"appointment": t, "counts": rng.binomial(SHOTS, probability).tolist()}
        )
    banks = []
    for kappa in (8, 0):
        for t, rows in histories(kappa, len(route), meter, sensor, (src, ref)):
            for name, joint in rows.items():
                for j, bond in enumerate((src, ref)):
                    p = joint[:, :, j].ravel().copy()
                    # Fixed roundoff grid for the multinomial seed, below the
                    # counting resolution. The last bin takes the residual.
                    p = np.maximum(0, np.round(p, 12))
                    p /= p.sum()
                    banks.append(
                        {
                            "kappa": kappa,
                            "appointment": t,
                            "menu": name,
                            "source": bond,
                            "records": LABELS[name],
                            "probabilities": p.tolist(),
                            "counts": rng.multinomial(SHOTS, p).tolist(),
                        }
                    )
    return {
        "shots": SHOTS,
        "seed": SEED,
        "sources": [src, ref],
        "route": route,
        "head": [route[meter - 1][1], route[meter][0]],
        "receiver": [route[sensor - 1][1], route[sensor][0]],
        "clocks": clocks,
        "outcomes_per_record": ["wall", "no wall"],
        "banks": banks,
        "readings": readback.decode(banks, clocks, SHOTS, (src, ref)),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    report = run()
    if parser.parse_args().write:
        (HERE / "records.json").write_text(json.dumps(report, indent=2) + "\n")
    else:
        print(json.dumps(report["readings"], indent=2))
