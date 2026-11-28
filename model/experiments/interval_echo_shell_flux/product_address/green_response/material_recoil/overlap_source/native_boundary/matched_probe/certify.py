"""Integer energy, onset, preparation, and source subalgebra checks on the closed board."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb, factorial

import numpy as np
from scipy.linalg import expm
import sympy as sp

from . import measure as m


def multiply(columns, previous):
    result = defaultdict(int)
    for i, v in previous.items():
        for j, h in columns[i]:
            result[j] += h * v
    return {i: v for i, v in result.items() if v}


def exact_blocks():
    board, seam, dist, origins = m.fp.make_board()
    preps = m.preparations(seam, origins)
    blocks = {}
    witnesses = []
    for total in range(3):
        h, words, bits = m.fp.generator(board, seam, total, 8, True)
        g = (4 * h).astype(int)
        assert (g != 4 * h).nnz == 0
        # Every off-diagonal entry stays in exactly one conserved block.
        q = np.repeat(bits[:, :12].sum(axis=1), 4) - np.tile([0, 0, 1, 1], len(bits))
        a, b = g.nonzero()
        assert np.array_equal(q[a], q[b])
        witnesses.append(
            dict(
                total=total,
                dimension=len(q),
                sectors={str(z): int(sum(q == z)) for z in np.unique(q)},
            )
        )
        columns = [
            list(
                zip(
                    g.indices[g.indptr[i] : g.indptr[i + 1]].tolist(),
                    g.data[g.indptr[i] : g.indptr[i + 1]].tolist(),
                    strict=True,
                )
            )
            for i in range(g.shape[0])
        ]
        radius = np.repeat(bits @ dist, 4).tolist()
        for (n, extra), keys in preps.items():
            if n + int(extra is not None) != total:
                continue
            ids = [4 * words.index(w) + c for w, c in keys]
            order = 99 if extra is None else 12
            powers = [[{i: 1}] for i in ids]
            moments = []
            for k in range(order + 1):
                moments.append([[u[-1].get(i, 0) for i in ids] for u in powers])
                if k == order:
                    break
                for u in powers:
                    u.append(multiply(columns, u[-1]))
                    # Only the first nine full vectors are needed for onset.
                    if k >= 9:
                        u[-2] = {}
            onset = []
            for p in range(0, 9, 2):
                matrix = []
                for u in powers:
                    row = []
                    for v in powers:
                        val = sum(
                            (-1) ** (p // 2 + j)
                            * comb(p, j)
                            * sum(
                                c * v[p - j].get(i, 0) * radius[i]
                                for i, c in u[j].items()
                            )
                            for j in range(p + 1)
                        )
                        row.append(F(val, 4**p * factorial(p)))
                    matrix.append(row)
                onset.append(matrix)
            blocks[n, extra] = dict(
                moments=np.array(moments, dtype=object),
                onset=np.array(onset, dtype=object),
                q=q[ids],
            )
    return blocks, witnesses, origins


def exact_average(arm, blocks, extra, key):
    # Every scalar in the stocks is a multiple of 1/76.
    total = 0
    for n, r in enumerate(arm):
        ints = np.rint(r * 76).astype(int)
        np.testing.assert_allclose(r, ints / 76, atol=1e-15, rtol=0)
        total = total + (blocks[n, extra][key] * ints.T).sum(axis=(1, 2))
    return [F(z, 76) for z in total]


def transfer():
    blocks, sectors, origins = exact_blocks()
    stocks = m.stocks()
    result = {}
    for name, arms in stocks.items():
        moments = []
        onset = []
        for extra in (None, *origins):
            rows = []
            for arm, target in zip(arms, stocks["target"], strict=True):
                x = exact_average(arm, blocks, extra, "moments")
                y = exact_average(target, blocks, extra, "moments")
                first = next((i for i in range(len(x)) if x[i] != y[i]), None)
                rows.append(
                    dict(
                        orders=len(x),
                        first_different_order=first,
                        difference_moment_of_4H=None
                        if first is None
                        else str(x[first] - y[first]),
                    )
                )
            moments.append(rows)
            low, high = [exact_average(arm, blocks, extra, "onset") for arm in arms]
            onset.append([end - start for end, start in zip(high, low, strict=True)])
        result[name] = dict(
            spectral_moments=moments,
            taylor_even=[
                [str(x - y) for x, y in zip(row, onset[0], strict=True)]
                for row in onset[1:]
            ],
        )
    pinched = {}
    for name, ref in [
        ("phase_feedback", "exchange_1_1"),
        ("balanced_coherent", "target"),
    ]:
        q = blocks[1, None]["q"]
        delta = sum(
            (-1) ** k * (stocks[name][k][1] - stocks[ref][k][1]) for k in (0, 1)
        )
        kept = delta * (q == q[:, None])
        assert np.max(np.abs(kept)) < 2e-15
        pinched[name] = dict(reference=ref, visible_delta_difference=0)
    # Q-visible real coordinates: four populations and Re rho[1,2].
    tomography = {}
    for group in ([None], [origins[0]], [origins[1]], origins, [None, *origins]):
        rows = []
        first = []
        for k in range(13):
            for p in group:
                op = blocks[1, p]["moments"][k]
                row = [*op.diagonal(), 2 * op[1, 2]]
                if sp.Matrix([*rows, row]).rank() > len(rows):
                    rows.append(row)
                    first.append(
                        dict(power=k, probe=p, coordinates=list(map(int, row)))
                    )
        tomography[",".join(map(str, group))] = dict(rank=len(rows), independent=first)
    return dict(
        sectors=sectors,
        spectral_boards_and_probes=[None, *origins],
        standards=result,
        pinched_identities=pinched,
        real_source_tomography=tomography,
    )


def preparations():
    law = m.source.LAW
    ids = [
        9 * (1 + law.LABELS.index((c, -1))) + law.SRC.index(w)
        for c in (0, 1)
        for w in ((1, 0), (0, 1))
    ]
    bright = np.zeros(99)
    bright[ids[1]] = bright[ids[2]] = 1 / np.sqrt(2)
    dark = bright.copy()
    dark[ids[2]] *= -1
    blank = np.eye(99)[:, law.SRC.index((0, 1))]
    u = law.echo()
    np.testing.assert_allclose(u @ bright, blank, atol=5e-16)
    np.testing.assert_allclose(u @ dark, dark, atol=5e-16)
    board, seam, _, _ = m.fp.make_board()
    a = seam[0]
    route = min(
        j
        for edge in board.seams
        if a in edge
        for j in edge
        if j != a and board.degree[j] == 5
    )
    # One mobile unit on (A, private-route port); B occupancy, if present,
    # means neither rail is occupied. The parked sector sees identity.
    swap = expm(-1j * np.pi * np.array([[8.5, -0.5], [-0.5, 8.5]]))
    private = np.diag([1.0, np.exp(-8j * np.pi / 8)])
    np.testing.assert_allclose(
        swap.conj().T @ private @ swap, np.diag([-1, 1]), atol=6e-15
    )
    phase = np.tile([(-1) ** a for a, _ in law.SRC], 11)
    np.testing.assert_allclose(phase * bright, dark, atol=0)
    stocks = m.stocks()
    f = stocks["phase_feedback"][1][1]
    e = stocks["target"][1][1]
    np.testing.assert_allclose(
        0.25 * m.pure(bright[ids]) + 0.75 * m.pure(dark[ids]), f, atol=2e-16
    )
    h4 = m.source.formation.emission.four_h(8)
    ex4 = m.source.formation.emission.four_h(0)
    # All three spectral levels are joint eigenlabels of the conjugated
    # exchange and contrast operators on this input's support.
    ex = ex4[np.ix_(ids, ids)]
    contrast = (h4 - ex4)[np.ix_(ids, ids)]
    for rho in (f, e):
        assert (rho * ex.T).sum() / 4 == 3 / 8
        assert (rho * contrast.T).sum() / 32 == 3 / 4
    # Coin and both binary ports as independently traceable physical factors.
    # The signed-label register is a common constant on every prepared word.
    emb = [2, 1, 6, 5]
    delta = np.zeros((8, 8))
    delta[np.ix_(emb, emb)] = f - e
    marginals = {}
    for size in range(1, 4):
        for keep in combinations(range(3), size):
            z = delta.reshape([2] * 6)
            n = 3
            for k in sorted(set(range(3)) - set(keep), reverse=True):
                z = z.trace(axis1=k, axis2=n + k)
                n -= 1
            z = z.reshape(2**n, 2**n)
            marginals["+".join(("coin", "A", "B")[k] for k in keep)] = float(
                abs(np.linalg.eigvalsh(z)).sum() / 2
            )
    # This is a sharp endpoint read followed by a seam-energy read on every
    # occupied outcome; empty outcomes contribute exactly zero.
    ready = np.tile(np.array([int(x == 1) for x, _ in law.SRC]), 11)
    gated = (ready[:, None] * h4 * ready[None, :])[np.ix_(ids, ids)]
    stats = [str(F((rho * gated.T).sum() / 4)) for rho in (e, f)]
    # Exchange and agreement commute on the original 99-state seam. Resolving
    # the five energy bins uniquely determines both labels on their support.
    _, projs = m.source.energy_projectors(8)
    for energy, proj, _ in projs:
        np.testing.assert_array_equal(ex4 @ proj, 4 * (energy % 8) * proj)
        np.testing.assert_array_equal(
            (h4 - ex4) @ proj, 4 * (energy - energy % 8) * proj
        )
    current = (
        1j
        * (h4[np.ix_(ids, ids)] / 4)
        * (np.array([1, 0, 1, 0])[None, :] - np.array([1, 0, 1, 0])[:, None])
    )
    return dict(
        phase_route=[a, route, a],
        energy_bright=float(bright @ h4 @ bright / 4),
        energy_dark=float(dark @ h4 @ dark / 4),
        high_marginal_trace_distances=marginals,
        high_seam_exchange="3/8",
        high_seam_unweighted_contrast="3/4",
        endpoint_gated_seam_energy_erased_feedback=stats,
        current_element_over_i=str(F(current[1, 2].imag)),
        visible_source_order=["c0 A", "c0 B", "c1 A", "c1 B"],
    )


def run():
    return dict(
        protocol_sha256=sha256(
            (m.HERE / "SOURCE_TOMOGRAPHY.md").read_bytes()
        ).hexdigest(),
        **transfer(),
        feedback=preparations(),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = run()
    path = m.HERE / "certificates.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        assert json.loads(path.read_text()) == out
    for name in ("target", "phase_feedback", "balanced_mass"):
        print(name, out["standards"][name])
    print(out["feedback"])
    print({k: v["rank"] for k, v in out["real_source_tomography"].items()})
