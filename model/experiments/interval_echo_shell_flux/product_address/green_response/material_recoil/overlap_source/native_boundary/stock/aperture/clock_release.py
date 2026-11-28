"""Initial accelerations of actual full-board material released beside the renewal clock."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np

from . import renewal_clock as clock
from . import thermal_radar

HERE = clock.HERE
lab = clock.lab


def instrument(board, home, marks, h, sites=None, radius=None):
    if radius is None:
        _, distance, _ = thermal_radar.chart(board, marks)
        radius = distance[home[1]]
    if sites is None:
        near = 13
        far = min(
            x
            for x in board.neighbours[near]
            if radius[x] == 2
            and board.degree[x] == board.degree[near]
            and not set(board.neighbours[x]) & set(home)
            and x not in np.ravel(marks)
        )
    else:
        near, far = sites
    mask = np.array([int(i == near) - int(i == far) for i in range(board.n)])
    read = np.column_stack((mask, radius))
    pot = [3**i for i in range(board.n)]
    labels = [3**board.n * 11**j for j in range(len(marks))]
    edges = sorted(
        set(marks)
        | {edge for edge in board.seams if set(edge) & set([*home, near, far])}
    )
    if sites is not None:
        # The first commutator with the two-port contrast only touches these
        # ports. The second reaches the neighbors of their neighbors.
        inner = {near, far, *board.neighbours[near], *board.neighbours[far]}
        edges = [e for e in edges if inner.intersection(e)]
    active = {
        i: [(int(k), int(h[k, i])) for k in np.flatnonzero(h[:, i])] for i in range(99)
    }

    def column(word, site, t):
        flat = dict(zip(home, map(int, word[:4]), strict=True))
        if site is not None:
            flat[site] = t
        heads = dict(zip(marks, map(int, word[4:]), strict=True))
        key = sum(pot[x] * v for x, v in flat.items()) + sum(
            a * b for a, b in zip(labels, word[4:].tolist(), strict=True)
        )
        reading = sum((q * read[x] for x, q in flat.items()), start=np.zeros(2, int))
        diag = sum(32 * q * q * (6 - int(board.degree[x])) for x, q in flat.items())
        out = {}
        for e in edges:
            x, y = e
            n, m = flat.get(x, 0), flat.get(y, 0)
            mark = heads.get(e, 0)
            code = mark * 9 + 3 * n + m
            for k, amp in active[code]:
                if k == code:
                    diag += amp
                    continue
                newhead, ports = divmod(k, 9)
                u, v = divmod(ports, 3)
                dst = key + pot[x] * (u - n) + pot[y] * (v - m)
                if newhead != mark:
                    dst += labels[marks.index(e)] * (newhead - mark)
                val = reading + (u - n) * read[x] + (v - m) * read[y]
                if dst in out:
                    old, other = out[dst]
                    assert np.array_equal(other, val)
                    amp += old
                out[dst] = amp, val
        out[key] = diag, reading
        return reading, out, key

    return near, int(far), column


def acceleration(columns, signs):
    # 2 Re <H psi|O H psi-H O psi>, with H tabulated in quarter units.
    amps = defaultdict(int)
    weighted = defaultdict(lambda: np.zeros(2, int))
    measures = {}
    for (initial, out, _), sign in zip(columns, signs, strict=True):
        for key, (amp, reading) in out.items():
            assert key not in measures or np.array_equal(measures[key], reading)
            measures[key] = reading
            amps[key] += sign * amp
            weighted[key] += sign * amp * initial
    num = sum(amp * (amp * measures[key] - weighted[key]) for key, amp in amps.items())
    return np.array([F(int(n), 8 * len(columns)) for n in num], object)


def run():
    _, w, cells, *_ = lab.energy.pair.setup()
    board, home, marks, _, _, four = lab.energy.returned(w)
    near, far, col = instrument(board, home, marks, lab.energy.pair.stock.kernels()[2])
    z = lab.energy.least_change(cells)
    standards = np.zeros((2, 7), int)
    standards[1, 1] = 1
    rows = {}
    for t in (1, 2):
        arrays = [[], []]
        for j, word in enumerate(np.vstack((w, standards))):
            base = col(word, None, t)
            if j < len(w):
                # Compare the full diagonal with the separately reconstructed
                # whole-board energy, on every original input.
                assert base[1][base[2]][0] == four[j]
            background = acceleration([base], [1])
            placed = [col(word, x, t) for x in (near, far)]
            for seq, sign in zip(arrays, (1, -1), strict=True):
                seq.append((acceleration(placed, [1, sign]) - background) / t)
        row = {}
        for key, seq in zip(("even", "odd"), arrays, strict=True):
            seq = np.array(seq, object)
            native = np.sum(seq[:-2], axis=0) / len(w)
            change = (z - 1) @ seq[:-2] / len(w)
            donor = seq[-1] - seq[-2]
            row[key] = dict(
                original=list(map(str, native)),
                repair_increment=list(map(str, change)),
                ordinary_donor_increment=list(map(str, donor)),
                equivalent_donors=[
                    str(a / b) if b else None
                    for a, b in zip(change, donor, strict=True)
                ],
            )
        rows[str(t)] = row
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "CLOCK_RELEASE.md",
                HERE / "renewal_clock.json",
                HERE / "closed_energy.json",
            )
        },
        near=near,
        far=far,
        observables=["near_minus_far_load", "radial_load_moment"],
        readers=rows,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "clock_release.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.check(json.loads(dst.read_bytes()), result)
    print(result["near"], result["far"], result["readers"])
