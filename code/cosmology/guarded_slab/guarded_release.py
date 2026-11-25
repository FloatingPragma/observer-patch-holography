"""An exact signed readout distribution at a locally triggered release."""

import argparse
import json
from collections import deque
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import numpy as np
from scipy.linalg import expm

from guarded_transfer import HERE, dimer, federation, observations, source

CAP = 16
P = F(1, 4)
LAMBDA = 1 - 2 * P


def within(board, origin, target, used):
    paths = deque([[origin]])
    seen = {origin}
    while paths:
        path = paths.popleft()
        if path[-1] == target:
            return path
        for nxt in board.neighbours[path[-1]]:
            if nxt // 12 == origin // 12 and nxt not in used and nxt not in seen:
                seen.add(nxt)
                paths.append(path + [nxt])
    raise ValueError(
        "Two gateways cannot exchange the tag around the occupied apparatus"
    )


def tag_receipts(start, deliveries, cut):
    # Restricted full-exchange columns move a sharp double and leave a blank pair blank.
    position = start
    seen = []
    for leg, delivery in enumerate(deliveries):
        path = delivery["tag_path"]
        for a, b in zip(path, path[1:]):
            if cut and leg == 0 and b == delivery["arrival"]:
                continue
            if position == a:
                position = b
            elif position == b:
                position = a
        assert position != delivery["vacant_reference"]
        seen.append(position == delivery["arrival"])
    return seen


def route(level):
    board = federation(level)
    lookup = dimer.layouts()
    reserved = set()
    for c in range(board.cells):
        off = 12 * c
        blocked = tuple(np.flatnonzero(board.degree[off : off + 12] == 6))
        reserved.update(off + i for i in lookup[blocked])
    links = {}
    for a, b in board.external:
        links[a // 12, b // 12] = a, b
        links[b // 12, a // 12] = b, a
    chain = [0]
    while len(chain) < 4:
        chain.append(min(v for a, v in links if a == chain[-1] and v not in chain))
    position = links[chain[0], chain[1]][0]
    start = position
    read = []
    for a, b in zip(chain, chain[1:]):
        leaving, arriving = links[a, b]
        path = within(board, position, leaving, reserved) + [arriving]
        assert not (set(path) & reserved)
        ref = next(
            j for j in board.neighbours[arriving] if j // 12 == b and j not in reserved
        )
        read.append(
            dict(receiver=b, tag_path=path, arrival=arriving, vacant_reference=ref)
        )
        position = arriving
    live = tag_receipts(start, read, cut=False)
    cut = tag_receipts(start, read, cut=True)
    assert all(live) and not any(cut)
    return dict(
        cells=board.cells,
        source=chain[0],
        parked_tag=start,
        deliveries=read,
        received=live,
        received_with_first_seam_parked=cut,
    )


def moments(k):
    z = LAMBDA**k
    return P * z / (1 - (1 - P) * z)


def word_probability(word, sign_means):
    result = F(1)
    for bit, m in zip(word, sign_means, strict=True):
        result *= (1 + (-1) ** bit * m) / 2
    return result


def compute(cap=CAP):
    weights = [P * (1 - P) ** (t - 1) for t in range(1, cap + 1)]
    timeout = (1 - P) ** cap
    words = list(product((0, 1), repeat=3))
    timed = []
    for t, w in enumerate(weights, 1):
        means = [LAMBDA ** (t + d) for d in (1, 2, 3)]
        timed.append([w * word_probability(word, means) for word in words])
    sharp = [sum(r[j] for r in timed) for j in range(8)]
    assert sum(sharp) + timeout == 1
    # For independent triggers each receiver retains its own miss digit (2).
    marginals = []
    for d in (1, 2, 3):
        sign = sum(w * LAMBDA ** (t + d) for t, w in enumerate(weights, 1))
        marginals.append(((1 - timeout + sign) / 2, (1 - timeout - sign) / 2, timeout))
    independent = {
        "".join(map(str, word)): np.prod([m[b] for m, b in zip(marginals, word)])
        for word in product(range(3), repeat=3)
    }
    assert sum(independent.values()) == 1
    variance = moments(2) - moments(1) ** 2
    third = moments(3) - 3 * moments(1) * moments(2) + 2 * moments(1) ** 3
    return dict(
        cap=cap,
        timeout=str(timeout),
        words=["".join(map(str, w)) for w in words],
        common_sharp_by_first_change=[[str(x) for x in r] for r in timed],
        common_sharp=[str(x) for x in sharp],
        common_stationary=[str((1 - timeout) / 8)] * 8,
        common_energy_by_first_change=[[str(w / 8)] * 8 for w in weights],
        common_energy=[str((1 - timeout) / 8)] * 8,
        cut_receivers={"222": "1"},
        independent_sharp={key: str(v) for key, v in independent.items()},
        unbounded_sign_mean=[str(LAMBDA**d * moments(1)) for d in (1, 2, 3)],
        unbounded_off_diagonal_cov=[
            str(LAMBDA ** (i + j) * variance) for i, j in ((1, 2), (1, 3), (2, 3))
        ],
        unbounded_connected_third=str(LAMBDA**6 * third),
    )


def energy_control():
    """Use the registered energy decoder on each sharp rail after tag arrival."""
    h, decoder, leak = dimer.gates("native")
    levels, modes = np.linalg.eigh(h[:2, :2])
    np.testing.assert_allclose(levels, [32, 33], atol=1e-13, rtol=0)
    np.testing.assert_allclose(abs(decoder[:2, :2] @ modes) ** 2, np.eye(2), atol=1e-13)
    receipts = {}
    for steps in dimer.STEPS:
        # Every first leak and both terminal labels are on the original denominator.
        expected = [F(1, 2 ** (j + 1)) for j in range(1, steps + 1)] + [
            F(1, 2),
            F(1, 2 ** (steps + 1)),
        ]
        for rail in (0, 1):
            np.testing.assert_allclose(
                dimer.outcomes(decoder, leak, rail, steps),
                list(map(float, expected)),
                atol=1e-13,
                rtol=0,
            )
        receipts[str(steps)] = list(map(str, expected))
    return receipts


def run(out):
    kernel, _ = source()
    np.testing.assert_allclose(kernel, [[0.75, 0.25], [0.25, 0.75]], atol=1e-13)
    law = observations.native_law
    h0, h1 = law.local_term_parts()
    cols = [law.SRC.index(w) for w in ((2, 0), (0, 2))]
    propagator = expm(-1j * np.pi * (h0 + observations.KAPPA * h1))
    pulse = propagator[:, cols]
    np.testing.assert_allclose(
        abs(pulse) ** 2, np.eye(law.DIM)[:, cols[::-1]], atol=1e-13
    )
    blank = law.SRC.index((0, 0))
    np.testing.assert_allclose(
        abs(propagator[:, blank]) ** 2, np.eye(law.DIM)[:, blank], atol=1e-13
    )
    dimer.pointer()
    data = compute()
    data["energy_reader_on_either_rail"] = energy_control()
    data["routes"] = [route(level) for level in range(3)]
    out.mkdir(parents=True, exist_ok=True)
    (out / "release.json").write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE)
    run(p.parse_args().out)
