"""Independent local-algebra and all-launch checks for the two deadline records."""

import json
from itertools import combinations

import numpy as np

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation
from .measure import HERE, SHOTS


def terminal():
    car = federation(0)
    echo = clock.law.echo()
    total = np.tile(
        [u + v for u, v in clock.law.SRC], clock.law.DIM // len(clock.law.SRC)
    )
    for m in range(5):
        # A coarse pointer retaining the endpoint sum after the echo and
        # uncomputing exactly preserves all coherence within its total sector.
        np.testing.assert_allclose(
            echo.T.conj() @ ((total == m)[:, None] * echo),
            np.diag((total == m).astype(float)),
            atol=2e-12,
            rtol=0,
        )
    # Every physically distinct local word: the double/single either
    # occupies a local port or is outside (code 12).
    words = np.array(
        [
            [2 * int(i == u) + int(i == v) for i in range(12)]
            for u in range(13)
            for v in range(13)
            if u != v or u == 12
        ]
    )
    decoded = np.empty_like(words)
    for i in range(12):
        neighbours = [j for j in car.neighbours[i] if j < 12]
        j, k = next(
            (j, k) for j, k in combinations(neighbours, 2) if k in car.neighbours[j]
        )
        a = words[:, i] + words[:, j]
        b = words[:, i] + words[:, k]
        c = words[:, j] + words[:, k]
        decoded[:, i] = (a + b - c) // 2
    np.testing.assert_array_equal(words, decoded)
    # Use the actual sender preaddress and certified pi contact to read
    # reset trials before any tick. Relay includes its addressed gateway pulse.
    for hand in [1, 2]:
        for relay in [False, True]:
            receivers = []
            for bit_port in (82, 74):
                word = np.zeros(car.n, dtype=int)
                word[[hand, 78]] = [2, 1]
                word[[78, bit_port]] = word[[bit_port, 78]]
                if relay:
                    word[[1, 82]] = word[[82, 1]]
                hit = np.where(word[:12] == 1)[0]
                receivers.append(int(hit[0]) if len(hit) else 12)
            if receivers != ([1, 12] if relay else [12, 12]):
                raise AssertionError("clock-independent reset law changed")


def counts():
    for name in ["records.json", "separated.json"]:
        data = json.loads((HERE / name).read_text())
        for source in ("passive", "relay", "half_exchange", "sham"):
            for d, row in data[source].items():
                keys = row["keys"]
                test = np.array(row["test"])
                train = np.array(row["train"])
                np.testing.assert_array_equal(test.sum(axis=1), [SHOTS] * 2)
                np.testing.assert_array_equal(train.sum(axis=1), [SHOTS] * 2)
                decisions = []
                for i, key in enumerate(keys):
                    _, bits, addr = key.split(":")
                    if addr == "timeout":
                        if len(bits) != 4 or bits.count("1") >= int(d):
                            raise AssertionError("timeout after deadline or before cap")
                        decisions.append(0)
                    else:
                        if bits[-1] != "1" or bits.count("1") != int(d):
                            raise AssertionError("answer not at enrolled receiver tick")
                        decisions.append(int(train[1, i] > train[0, i]))
                correct = sum(int(test[g, i]) for i, g in enumerate(decisions))
                if correct != row["correct"]["joint"]:
                    raise AssertionError("full decoder score changed")
    print(
        "Local triangle census, reset comparator and complete receiver decisions verified."
    )


if __name__ == "__main__":
    terminal()
    counts()
