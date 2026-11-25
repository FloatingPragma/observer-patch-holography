"""Joint local clock words and message decisions from two freely evolving loads."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation
from model.experiments.material_inertia.neutral_defect import generators, ordered_index

HERE = Path(__file__).resolve().parent
SHOTS = 2**17
SLOTS = 4
STEP = 0.5
NEAR, FAR, COMMON, OTHER = 1, 82, 78, 74


def native_gates():
    law = clock.law
    echo = law.echo()
    h0, h1 = law.local_term_parts()
    h = h0 + clock.KAPPA * h1
    pi = expm(-1j * np.pi * h)
    pointer = np.tile([p == (1, 1) for p in law.SRC], law.DIM // len(law.SRC))
    for a, b in product(range(3), repeat=2):
        ket = np.eye(law.DIM)[:, law.SRC.index((a, b))]
        np.testing.assert_allclose(
            echo.conj().T @ (pointer * (echo @ ket)),
            int(a + b == 2) * ket,
            rtol=0,
            atol=1e-12,
        )
        swapped = np.eye(law.DIM)[:, law.SRC.index((b, a))]
        np.testing.assert_allclose(pi @ ket, swapped, atol=1e-12, rtol=0)


def flight():
    car = federation(0)
    # x and y are now the double and the single of an otherwise empty screen.
    x, y, adj, k, _, _ = generators(car)
    interaction = -4 * clock.KAPPA * adj.astype(float)
    h = k + sparse.diags(interaction)
    for a, b in [(1, 0), (82, 74), (82, 78), (78, 74), (82, 1)]:
        if b not in car.neighbours[a]:
            raise ValueError("station map changed")
    # Independently assemble representative full blank-memory columns, including
    # reference vacancies; subtract 6*kappa*(2**2+1**2) from the entire sector.
    h0, h1 = clock.law.local_term_parts()
    seam = h0 + clock.KAPPA * h1
    for u, v in [(1, 82), (1, 74), (82, 1), (1, 0)]:
        word = np.zeros(car.n, dtype=int)
        word[[u, v]] = [2, 1]
        i = ordered_index(u, v, car.n)
        col = np.zeros(len(x))
        col[i] = clock.KAPPA * ((6 - car.degree) @ (word**2) - 30)
        for a, b in car.seams:
            m, n = word[a], word[b]
            index = clock.law.SRC.index((m, n))
            col[i] += seam[index, index]
            if m != n:
                swap = {a: b, b: a}
                j = ordered_index(swap.get(u, u), swap.get(v, v), car.n)
                col[j] += seam[clock.law.SRC.index((n, m)), index]
        np.testing.assert_allclose(col, h[:, i].toarray().ravel(), atol=1e-11, rtol=0)
    tick = ((x == NEAR) | (x == 0)) & ~((y == NEAR) | (y == 0))
    return car, x, y, tick, k, sparse.diags(interaction, format="csr")


def branched(h, car, y, tick, starts, *, ranging=False, turn=None):
    """Take one shared evolution per prefix. Each column is one input bit or arm."""
    state = np.zeros((len(y), len(starts)), complex)
    for j, (double, single) in enumerate(starts):
        state[ordered_index(double, single, car.n), j] = 1
    banks = {d: {} for d in (1, 3)}
    terminal = np.minimum(y, 12)
    prop = (-1j * STEP * h).tocsr()
    trace = prop.diagonal().sum()

    def enter(word, psi, mark):
        out = expm_multiply(prop, psi, traceA=trace)
        for digit, select in [("0", ~tick), ("1", tick)]:
            bits = word + digit
            part = select[:, None] * out
            done = len(bits) == SLOTS
            for d, bank in banks.items():
                reached = digit == "1" and bits.count("1") == d
                expired = done and bits.count("1") < d
                if reached or expired:
                    mass = abs(part) ** 2
                    if reached:
                        for j in range(mass.shape[1]):
                            counts = np.bincount(
                                terminal, weights=mass[:, j], minlength=13
                            )
                            for address, p in enumerate(counts):
                                key = f"{mark}:{bits}:{address}"
                                bank.setdefault(key, np.zeros(len(starts)))[j] += p
                    else:
                        bank[f"{mark}:{bits}:timeout"] = mass.sum(axis=0)
            if not done and bits.count("1") < 3:
                if ranging and len(bits) == 1:
                    for tag, mask in [("visited", y == FAR), ("missed", y != FAR)]:
                        branch = part * mask[:, None]
                        if tag == "visited" and turn is not None:
                            branch = branch[turn]
                        enter(bits, branch, tag)
                else:
                    enter(bits, part, mark)

    enter("", state, "unmarked")
    for bank in banks.values():
        np.testing.assert_allclose(sum(bank.values()), 1.0, atol=2e-11, rtol=0)
    return banks


def sample(law, rng):
    p = np.array(law).round(13)
    np.testing.assert_allclose(p.sum(), 1, rtol=0, atol=2e-10)
    return rng.multinomial(SHOTS, p / p.sum())


def readers(keys, train, test):
    outputs = {}
    for name in ("joint", "clock", "address"):
        groups = {}
        for i, key in enumerate(keys):
            _, word, address = key.split(":")
            label = key if name == "joint" else (word if name == "clock" else address)
            groups.setdefault((address == "timeout", label), []).append(i)
        passed = 0
        for (expired, _), indices in groups.items():
            guess = int(
                not expired and train[1, indices].sum() > train[0, indices].sum()
            )
            passed += int(test[guess, indices].sum())
        outputs[name] = passed
    return outputs


def acquire(banks, rng):
    result = {}
    for d, bank in banks.items():
        keys = sorted(bank)
        probs = np.array([bank[key] for key in keys]).T
        train = np.array([sample(p, rng) for p in probs])
        test = np.array([sample(p, rng) for p in probs])
        timeout = np.array([k.endswith(":timeout") for k in keys])
        detect = np.array([not k.endswith((":timeout", ":12")) for k in keys])
        result[str(d)] = dict(
            keys=keys,
            train=train.tolist(),
            test=test.tolist(),
            correct=readers(keys, train, test),
            expired=test[:, timeout].sum(axis=1).tolist(),
            detected=test[:, detect].sum(axis=1).tolist(),
        )
    return result


def run():
    native_gates()
    car, _, y, tick, k, v = flight()
    h = k + v
    # A native pi pulse on the sender is a permutation in this sharp basis.
    # The additional relay permutes gateway addresses after the sender pulse.
    native = branched(h, car, y, tick, [(1, 82), (1, 74), (82, 1), (82, 74)])
    rng = np.random.default_rng(305219)
    result = {}
    for name, columns in [("passive", [0, 1]), ("relay", [2, 3]), ("sham", [0, 0])]:
        result[name] = acquire(
            {
                d: {key: p[columns] for key, p in bank.items()}
                for d, bank in native.items()
            },
            rng,
        )
    result["half_exchange"] = acquire(
        branched(k * 0.5 + v, car, y, tick, [(1, 82), (1, 74)]), rng
    )
    marks = branched(h, car, y, tick, [(1, 0)], ranging=True)
    result["range"] = {}
    for d, bank in marks.items():
        keys = sorted(bank)
        counts = sample([bank[key][0] for key in keys], rng)

        def kept(key):
            tag, _, addr = key.split(":")
            return tag == "visited" and addr not in ("12", "timeout")

        returns = sum(int(n) for key, n in zip(keys, counts) if kept(key))
        result["range"][str(d)] = dict(
            keys=keys, counts=counts.tolist(), returned_via_peer=returns
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = run()
    if args.record:
        with args.record.open("x") as f:
            json.dump(result, f, separators=(",", ":"))
            f.write("\n")
    elif result != json.loads((HERE / "records.json").read_text()):
        raise SystemExit("joint readout differs from retained count banks")
    for arm, data in result.items():
        print(
            arm,
            {
                d: (
                    b.get("correct"),
                    b.get("expired"),
                    b.get("detected"),
                    b.get("returned_via_peer"),
                )
                for d, b in data.items()
            },
        )
