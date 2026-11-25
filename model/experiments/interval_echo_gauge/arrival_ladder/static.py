"""Release an encoded charge and a material donor under one static native Hamiltonian."""

import argparse
from functools import cache
import json

import numpy as np
from scipy.linalg import expm

from . import measure as m


def geometry():
    bench = m.layout()
    board = m.federation(0)
    a, comparator = bench["a"], bench["donor"]
    b = min(
        p
        for p in board.neighbours[a]
        if p != bench["guard"]
        and p < 12
        and p in board.neighbours[comparator]
        and board.degree[p] == board.degree[a]
        and p not in bench["compensators"]
        and p not in bench["neutral_path"]
    )
    return dict(
        **{k: v for k, v in bench.items() if k not in ("donor", "offline")},
        donor=b,
        comparator=comparator,
        offline=sorted(set(bench["offline"]) | {b}),
    )


@cache
def collar(kappa, parity):
    exchange, contrast = m.LAW.local_term_parts()
    h = exchange + kappa * contrast
    cols = np.kron(m.MEMORIES[:, [parity]], np.eye(9))
    values = np.diag(cols.conj().T @ h @ cols).real
    np.testing.assert_allclose(h @ cols, cols * values, atol=1e-13)
    return values.reshape(3, 3)


@cache
def block(kappa, parity, encoded_x, encoded_s):
    site = geometry()
    board = m.federation(0)
    a, b, r = site["a"], site["donor"], site["comparator"]
    private = 6 - board.degree
    diagonal = []
    for na, nb in m.LAW.SRC:
        energy = kappa * (
            private[a] * na**2
            + private[b] * nb**2
            + private[r]
            + private[site["x"]] * encoded_x**2
            + private[site["guard"]] * encoded_s**2
        )
        energy += collar(kappa, 0)[encoded_x, na] + collar(kappa, parity)[nb, 1]
        diagonal.append(energy)
    return m.free_seam(kappa) + np.diag(diagonal)


def instrument(word, time, kappa, parity):
    x, a, b = word
    echo = m.LAW.echo()
    encoded = echo[:, m.LAW.SRC.index((x, 2))]
    state = np.zeros((99, 9), complex)
    for i in np.flatnonzero(abs(encoded) > 1e-14):
        ex, es = m.LAW.SRC[i % 9]
        state[i] = (
            encoded[i]
            * expm(-1j * time * block(kappa, parity, ex, es))[
                :, m.LAW.SRC.index((a, b))
            ]
        )
    state = echo @ state
    decoded = abs(state[[m.LAW.SRC.index((ox, 2)) for ox in range(3)]]) ** 2
    prob = decoded.ravel()
    # Departure from the restored blank sign collar or the double reference
    # is a resolved rejection, even if the two target loads happen to match.
    reject = float(np.sum(abs(np.delete(state, [2, 5, 8], axis=0)) ** 2))
    np.testing.assert_allclose(sum(prob) + reject, 1, atol=1e-12, rtol=0)
    return np.r_[prob, reject]


def measure():
    site = geometry()
    shots = 2**30
    rng = np.random.default_rng(5405)
    rows = []
    for kappa, parity in [(m.K, 0), (m.K, 1), (m.K / 2, 0)]:
        for time in (0, np.pi / 2, np.pi, 2 * np.pi):
            for source in (
                "negative",
                "neutral",
                "positive",
                "filled",
                "donor_depleted",
            ):
                probs = instrument(m.SOURCES[source], time, kappa, parity)
                rounded = np.round(probs, 14)
                counts = rng.multinomial(shots, rounded / rounded.sum())
                rows.append(
                    dict(
                        kappa=kappa,
                        comparator_parity=parity,
                        source=source,
                        time_over_pi=float(time / np.pi),
                        probabilities=probs.tolist(),
                        counts=counts.tolist(),
                    )
                )
    return dict(
        schema=1,
        site=site,
        shots_per_prepared_input=shots,
        arrival_chances={
            s: m.arrival(site, sign, m.K)
            for s, sign in (("negative", -1), ("positive", 1))
        },
        code_words=m.WORDS,
        restored_then_rejected_banks=rows,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(measure(), indent=2) + "\n"
    if args.write:
        (m.HERE / "static.json").write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
