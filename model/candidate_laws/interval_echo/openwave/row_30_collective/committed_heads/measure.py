"""Resolved connector-head records followed by two-owner wall measurements."""

import json
from pathlib import Path
from itertools import product
from functools import cache
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label
from .. import experiment as strand
from . import readback

HERE = Path(__file__).resolve().parent
TIMES = (0, 2, 4, 6)
SHOTS = 2**25
SEED = 351132
NAMES = ("magnitude", "sign", "coin", "full")
TAIL = readback.TAIL
PARITIES = list(product((1, -1), repeat=2))


@cache
def readings():
    law = strand.candidate.native_law
    ex, contrast = law.local_term_parts()
    initial = label(2, 1, 1)
    prepared = (
        law.echo()
        @ (np.eye(99)[:, law.SRC.index((2, 0))] + np.eye(99)[:, law.SRC.index((0, 2))])
        / np.sqrt(2)
    )
    np.testing.assert_allclose(
        prepared, np.kron(initial, np.eye(9)[:, law.SRC.index((1, 1))]), atol=1e-14
    )
    bases = np.stack([label(2, p, c) for p, c in PARITIES])
    block = {}
    for k in (0, 8):
        for p, c in PARITIES:
            B = np.kron(label(2, p, c)[:, None], np.eye(9))
            h = ex + k * contrast
            expected = np.array(
                [4 * k + (1 - p) / 2 if w != (1, 1) or c == -1 else 0 for w in law.SRC]
            )
            np.testing.assert_allclose(h @ B, B * expected, atol=3e-13, rtol=0)
            block[k, p, c] = expected.reshape(3, 3)
    masks = {name: {} for name in NAMES}
    for m, (coin, sgn) in enumerate(law.LABELS, start=1):
        if abs(sgn) != 2:
            continue
        for name, val in zip(
            NAMES, ("2", str(sgn), str(coin), f"{coin}/{sgn}"), strict=True
        ):
            masks[name].setdefault(val, []).append(m)
    probs = {}
    for name, outcomes in masks.items():
        probs[name] = {}
        for val, idx in outcomes.items():
            state = np.zeros_like(initial)
            state[idx] = initial[idx]
            weights = abs(bases @ state) ** 2
            probs[name][val] = weights
        np.testing.assert_allclose(
            sum(z.sum() for z in probs[name].values()), 1, atol=1e-14
        )
    return probs, block


def geometry():
    board = federation(0)
    route = max((x for cycle, x in strand.grow(board, 312141) if not cycle), key=len)
    n = len(route)
    src, ref = n // 3, 1

    def owner(i):
        return {route[i - 1][1] // 12, route[i][0] // 12}

    sensor = next(
        i
        for i in range(1, n)
        if owner(i).isdisjoint(owner(src) | owner(ref))
        and min(abs(i - src), abs(i - ref)) >= 3
    )
    meter = max(src, ref) + 1
    assert sensor > meter and n == 14
    return route, src, ref, meter, sensor


def initial_word(n, j):
    last = 1
    word = 0
    for i in range(n):
        last = last if i == j else -last
        word += (1 << i) if last == -1 else 0
    return word


def probabilities(k, route, sources, meter, sensor):
    n = len(route)
    words = np.arange(2**n)
    bits = (words[:, None] >> np.arange(n)) & 1  # 1 means vacancy on left
    active = (bits[:, :-1] == 1) & (bits[:, 1:] == 0)
    remote = bits[:, sensor - 1] == bits[:, sensor]
    flips = np.concatenate([words ^ (1 << i) for i in range(n)])
    current = sparse.csr_matrix(
        (-0.5 * np.ones(n * len(words)), (np.tile(words, n), flips)),
        shape=(len(words),) * 2,
    )
    psi = np.zeros((len(words), 2))
    psi[sources, np.arange(2)] = 1
    assert np.all(psi[remote] == 0)
    _, energy = readings()
    result = []
    for p, c in PARITIES:
        if (p, c) == (-1, -1):
            # The two cut Hamiltonians differ by a scalar.
            result.append(result[1])
            continue
        diag = 4 * k * (n - 1 - active.sum(axis=1))
        before = energy[k, 1, 1][bits[:, meter - 1], 1 - bits[:, meter]]
        after = energy[k, p, c][bits[:, meter - 1], 1 - bits[:, meter]]
        diag = diag - before + after
        if k == 0:
            # Only the addressed pair of dimers can interact. The remote
            # dimer pair is disjoint and each exchanges on its own blank seam.
            assert sensor - 1 > meter
            data = (np.sin(TIMES) ** 2 / 2)[:, None] * np.ones((1, 2))
        elif c == -1:
            # Reduce both product sources to the same right-hand state, and
            # propagate the full right factor of the cut physical Hamiltonian.
            right = n - meter
            tags = words[: 2**right] << meter
            assert sources[0] >> meter == sources[1] >> meter
            h = current[tags, :][:, tags] + sparse.diags(diag[tags] - diag[tags].mean())
            right_init = np.eye(1, 2**right, sources[0] >> meter).ravel()
            out = expm_multiply(-1j * h, right_init, start=0, stop=6, num=4)
            data = (abs(out[:, remote[tags]]) ** 2).sum(axis=1)[:, None] * np.ones(
                (1, 2)
            )
        else:
            h = current + sparse.diags(diag - diag.mean())
            out = expm_multiply(-1j * h, psi, start=0, stop=6, num=4)
            np.testing.assert_allclose(
                (abs(out) ** 2).sum(axis=1), 1, atol=5e-11, rtol=0
            )
            data = (abs(out[:, remote, :]) ** 2).sum(axis=1)
        result.append(data)
    return np.stack(result)


def run():
    strand.native_checks()
    route, src, ref, meter, sensor = geometry()
    sources = [initial_word(len(route), i) for i in (src, ref)]
    rng = np.random.default_rng(SEED)
    bank = []
    forecast = {}
    probs, _ = readings()
    for k in (8, 0):
        msg = probabilities(k, route, sources, meter, sensor)
        accuracies = {}
        for name, out in probs.items():
            expectations = []
            for it, t in enumerate(TIMES):
                theory = []
                for a in range(len(sources)):
                    joint = []
                    for val, weights in out.items():
                        click = float(weights @ msg[:, it, a])
                        joint.extend([click, weights.sum() - click])
                    dist = np.array(joint)
                    assert dist.min() > -1e-12 and abs(dist.sum() - 1) < 1e-12
                    dist = np.maximum(dist, 0)
                    dist /= dist.sum()
                    counts = rng.multinomial(SHOTS, dist)
                    # Separate clock cohort; the contact is its own clock. Read
                    # times modulo one turn, unwrapped in appointment order.
                    hand = rng.binomial(
                        SHOTS, (1 + np.array([np.cos(t), np.sin(t)])) / 2
                    )
                    theory.append(float(dist[::2].sum()))
                    bank.append(
                        dict(
                            k=k,
                            instrument=name,
                            t=t,
                            source_bond=(src, ref)[a],
                            record=list(out),
                            counts=counts.tolist(),
                            clock=hand.tolist(),
                            probabilities=dist.tolist(),
                        )
                    )
                expectations.append(theory[0] - theory[1])
            accuracies[name] = dict(source_response=expectations)
        forecast[str(k)] = accuracies
    return dict(
        route=route,
        authority="20-cell established open branch, seed 312141",
        seed=SEED,
        shots=SHOTS,
        sources=[src, ref],
        source_orientations=sources,
        head_index=meter,
        head_addresses=[route[meter - 1][1], route[meter][0]],
        sensor_index=sensor,
        sensor_addresses=[route[sensor - 1][1], route[sensor][0]],
        conditional_counts="each record in indicated order, remote wall then remote no wall",
        tail=TAIL,
        count_readback=readback.evaluate(bank, SHOTS, (src, ref), TIMES),
        summary=forecast,
        banks=bank,
    )


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    d = run()
    if args.write:
        (HERE / "results.json").write_text(json.dumps(d, indent=2) + "\n")
    else:
        print(d["summary"])
