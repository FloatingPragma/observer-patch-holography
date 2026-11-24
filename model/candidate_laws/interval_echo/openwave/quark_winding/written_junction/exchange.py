"""Free exchange of two density-wall cores across one native memory contact."""

import argparse
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier
from model.experiments.interval_echo_encounter import memory
from ..belt_walls import belt


def geometry():
    board = carrier.federation(0)
    a, b = min(board.external)
    layouts = []
    for u in (a, b):
        order = next(v for v in belt.belts().values() if u % 12 in v)
        layouts.append([u // 12 * 12 + i for i in order])
    return board, layouts, [a, b]


def words(n):
    return [sum(1 << i for i in sites) for sites in combinations(range(10), n)]


def walls(n):
    return {
        word
        for word in words(n)
        if sum(
            bool(word & 1 << a) and bool(word & 1 << b) for a, b in belt.ring(10).seams
        )
        == (0 if n == 3 else 1)
    }


def sources(tips):
    three, four = walls(3), walls(4)
    matches = [
        (a, b)
        for a, b in product(sorted(four), sorted(three))
        if a & tips[0]
        and not b & tips[1]
        and a ^ tips[0] in three
        and b ^ tips[1] in four
    ]
    return matches


def side(basis, layout, board, kappa):
    index = {w: i for i, w in enumerate(basis)}
    rows, cols, values = [], [], []
    for i, w in enumerate(basis):
        d = kappa * sum(
            (6 - board.degree[p]) for j, p in enumerate(layout) if w & 1 << j
        )
        for a, b in belt.ring(10).seams:
            if bool(w & 1 << a) != bool(w & 1 << b):
                d += kappa + 0.5
                rows.append(i)
                cols.append(index[w ^ (1 << a) ^ (1 << b)])
                values.append(-0.5)
        rows.append(i)
        cols.append(i)
        values.append(d)
    return sparse.csr_matrix((values, (rows, cols)), shape=(len(basis), len(basis)))


def lifted(kappa, marked, board, layouts, tips):
    states = []
    blocks = []
    for nr in (3, 4) if marked else range(8):
        left, right = words(7 - nr), words(nr)
        states.extend(product(left, right))
        blocks.append(
            sparse.kronsum(
                side(right, layouts[1], board, kappa),
                side(left, layouts[0], board, kappa),
                format="csr",
            )
        )
    h = sparse.block_diag(blocks, format="csr")
    idx = {pair: i for i, pair in enumerate(states)}
    nr = np.array([b.bit_count() for _, b in states])
    three, four = walls(3), walls(4)
    first = np.array([a in four and b in three for a, b in states])
    second = np.array([a in three and b in four for a, b in states])
    members, base, _, _ = memory.local_components(kappa)
    matrices = []
    for p in (0, 1) if marked else (0,):
        mat = members[1 + p][1] if marked else base
        rows, cols, values = [], [], []
        for j, (a, b) in enumerate(states):
            coin = 4 - int(nr[j]) if marked else 0
            bits = (int(bool(a & tips[0])), int(bool(b & tips[1])))
            k = 4 * coin + memory.BITS.index(bits)
            for target in np.flatnonzero(abs(mat[:, k]) > 1e-12):
                out_coin, out_bits = divmod(int(target), 4)
                newbits = memory.BITS[out_bits]
                aa = a ^ tips[0] if bits[0] != newbits[0] else a
                bb = b ^ tips[1] if bits[1] != newbits[1] else b
                if marked:
                    assert out_coin == 4 - bb.bit_count()
                rows.append(idx[aa, bb])
                cols.append(j)
                values.append(mat[target, k])
        joint = h + sparse.csr_matrix((values, (rows, cols)), shape=h.shape)
        assert sparse.linalg.norm(joint - joint.T) < 1e-9
        matrices.append(joint)
    return idx, matrices, nr, first, second


def flight(kappa, marked, board, layouts, tips, launch):
    idx, matrices, nr, first, second = lifted(kappa, marked, board, layouts, tips)
    a, b = launch
    # Actual 99-state source, with both parities retained until the final memory read.
    write = memory.local_components(kappa)[0]
    source = []
    for parity, h in enumerate(matrices):
        psi = np.zeros(len(idx), complex)
        if marked:
            local = write[parity + 1][2][:, memory.BITS.index((1, 0))]
            for k in np.flatnonzero(abs(local) > 1e-12):
                _, bits = divmod(int(k), 4)
                pair = (
                    (a, b)
                    if memory.BITS[bits] == (1, 0)
                    else (a ^ tips[0], b ^ tips[1])
                )
                psi[idx[pair]] = local[k]
        else:
            psi[idx[a, b]] = 1
        source.append(psi)
    panels = []
    previous = 0.0
    for t in (0.0, 1 / 16, 0.25, 1.0, 4.0):
        if t > previous:
            source = [
                expm_multiply(-1j * (t - previous) * h, psi)
                for h, psi in zip(matrices, source)
            ]
        previous = t
        prob = sum(abs(psi) ** 2 for psi in source)
        np.testing.assert_allclose(prob.sum(), 1, atol=1e-10)
        signs = {}
        if marked:
            plus, minus = source
            for name, sgn in [("positive", 1), ("negative", -1)]:
                signprob = abs((plus + sgn * minus) / np.sqrt(2)) ** 2
                signs[name] = [
                    [
                        float(signprob[(nr == n) & mask].sum())
                        for mask in (first, second, ~(first | second))
                    ]
                    for n in (3, 4)
                ]
        panels.append(
            dict(
                time=t,
                original=prob[first].sum(),
                exchanged=prob[second].sum(),
                other=prob[~(first | second)].sum(),
                right_count=[prob[nr == i].sum() for i in range(8)],
                sign_by_right_count_and_wall=signs,
            )
        )
    return dict(block_dimensions=[len(idx)] * len(matrices), reads=panels)


def run():
    board, layouts, edge = geometry()
    tips = [1 << layout.index(i) for i, layout in zip(edge, layouts)]
    pairs = sources(tips)
    out = {}
    for k in (candidate.KAPPA, 0):
        for marked in (True, False):
            out[f"k{int(k)}_" + ("written" if marked else "erased")] = flight(
                k, marked, board, layouts, tips, pairs[0]
            )
    return currents._portable(
        dict(
            edge=edge,
            layouts=layouts,
            sharp_pairs_tested=100,
            exchange_pairs=pairs,
            switched_word=[pairs[0][i] ^ tips[i] for i in range(2)],
            runs=out,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(data)
