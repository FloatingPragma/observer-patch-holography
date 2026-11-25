"""A localized oriented shift charge, its boundary transport, and the force on two marked heads."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .. import bridge
from .chart import Atlas

HERE = Path(__file__).resolve().parent


def map_matrix(board, patch):
    a, b = [], []
    for u, v, perm in zip(patch["cells"], patch["images"], patch["lifts"]):
        a.extend(12 * v + np.asarray(perm))
        b.extend(range(12 * u, 12 * u + 12))
    s = sparse.csr_matrix(
        (np.ones(len(a), dtype=np.int64), (a, b)), shape=(board.n, board.n)
    )
    return sparse.kron(s, sparse.eye(4, dtype=np.int64), format="csr")


def blank_h(board, q):
    h = sparse.diags(
        17 * board.degree + 16 * (2 * q + q * q) * (6 - board.degree), dtype=np.int64
    ) - board.adjacency().astype(np.int64)
    return sparse.kron(h, sparse.eye(4, dtype=np.int64), format="csr")


def entries(mat):
    mat = mat.tocsr()
    mat.eliminate_zeros()
    return mat


def op_table(board, route, patch, q, odd):
    h = bridge.native.generator(board, route, q, odd)
    bg = blank_h(board, q)
    s = map_matrix(board, patch)
    r = s.T @ s
    image = s @ s.T
    ident = sparse.eye(s.shape[0], dtype=np.int64, format="csr")
    assert (r @ r - r).count_nonzero() == 0 and (
        image @ image - image
    ).count_nonzero() == 0
    assert (image @ bg @ image @ s - s @ r @ bg @ r).count_nonzero() == 0
    counted = sorted(set(route))
    v = h - bg
    offset = int(v.diagonal()[4 * next(p for p in range(board.n) if p not in counted)])
    v = v - offset * ident
    sites = {int(i) // 4 for a in v.nonzero() for i in a}
    assert sites == set(route)
    # 2i times the unscaled sine of the partial translation.
    skew = s - s.T
    # For Q=(S-S†)/2i, 4 dot Q = [2H,S] + [2H,S]†.
    out = (ident - image) @ bg @ s
    inc = -s @ bg @ (ident - r)
    bound = out + inc
    assert (bg @ s - s @ bg - bound).count_nonzero() == 0
    marked = v @ s - s @ v
    full = h @ s - s @ h
    assert (full - marked - bound).count_nonzero() == 0
    ops = {
        "shift": (skew, 2j),
        "domain_cut": (inc + inc.T, 4),
        "image_cut": (out + out.T, 4),
        "boundary": (bound + bound.T, 4),
        "marked_force": (marked + marked.T, 4),
        "free_rate": (full + full.T, 4),
    }
    return h, {k: (entries(m), unit) for k, (m, unit) in ops.items()}


def flux_on_edges(board, patch):
    """Integer numerator of each outward-paired source and image cut read."""
    peer = {p: t for a, b in board.external for p, t in ((a, b), (b, a))}
    edges = {}
    eye = sparse.eye(4, dtype=np.int64)
    for u, image, perm in zip(patch["cells"], patch["images"], patch["lifts"]):
        for p in range(12 * u, 12 * u + 12):
            if p not in peer:
                continue
            a = 12 * image + perm[p % 12]
            b, t = peer[a], peer[p]
            # Native off-diagonal 2H=-1. Both port pairs are read in both
            # interference signs and act as identity on all head states.
            sites = sparse.csr_matrix(
                ([-1, -1, 1, 1], ([b, p, a, t], [p, b, t, a])), shape=(board.n, board.n)
            )
            edges[u, t // 12] = entries(sparse.kron(sites, eye, format="csr"))
    for (u, v), op in edges.items():
        if (v, u) in edges:
            assert (op + edges[v, u]).count_nonzero() == 0
    return edges


def regional_energy(board, route, patch, q, odd, h):
    """Active internal and half-bridge energy on exactly the same domain as Q."""
    selected = set(patch["cells"])
    cells = bridge.operators(board, route, q, odd)
    e4 = 2 * bridge.masked(cells, selected)
    interior = []
    boundary = []
    for a, b in board.external:
        ua, ub = a // 12, b // 12
        if ua in selected and ub in selected:
            interior.append((a, b))
        if (ua in selected) != (ub in selected):
            boundary.append((a, b))
    e4 += 2 * bridge.bridges(board, interior) + bridge.bridges(board, boundary)
    # Equal scalar shift of the complete Hamiltonian cancels in the commutator.
    summed = entries(h @ e4 - e4 @ h)
    per = {}
    for a, b in board.external:
        u, v = a // 12, b // 12
        if u not in selected and v not in selected:
            continue
        cur = bridge.current(board, cells, (a, b))
        if u in selected:
            per[u, v] = -cur
        if v in selected:
            per[v, u] = cur
    assert (summed - sum(per.values())).count_nonzero() == 0
    return e4, summed, per


def chosen(atlas):
    # First source cell and first adjacent two-step shift admitting the largest
    # witness in the specified orientation; choice recorded before state reads.
    out = atlas.patch(0, 19, 0)
    assert len(out["cells"]) == 12 and {0, 3} <= set(out["cells"])
    return out


def finite_geometry():
    rows = {}
    for level in range(3):
        atlas = Atlas(level)
        counts, balls = atlas.census()
        best = 0
        for u in sorted(atlas.graph):
            for v, d in atlas.dist[u].items():
                if d == 2:
                    best = max(
                        best, *[len(atlas.patch(u, v, k)["cells"]) for k in range(3)]
                    )
        rows[atlas.board.cells] = dict(
            complete_ball_directions=sorted(counts.items()),
            maximum_witnessed_connected_cells=best,
            largest_complete_ball_cells=len(balls[max(balls)]["cells"]) if balls else 0,
        )
    return rows


def run():
    atlas = Atlas(1)
    patch = chosen(atlas)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    pos = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    states = bridge.tangents()
    rows = {}
    local_rows = {}
    currents = flux_on_edges(board, patch)
    for q in (-1, 1):
        materials = {}
        local = {}
        for odd in (False, True):
            h, ops = op_table(board, route, patch, q, odd)
            start = np.zeros((h.shape[0], len(pos)))
            start[pos] = np.eye(len(pos))
            ht = -2j * h
            w = expm_multiply(ht, start, traceA=ht.diagonal().sum())

            def effect(v, t):
                z = t @ v
                return np.mean(
                    [
                        a.conj().T @ b
                        for a, b in zip(np.hsplit(v, 10), np.hsplit(z, 10))
                    ],
                    axis=0,
                )

            for k, (m, unit) in ops.items():
                materials.setdefault(k, []).append(effect(w, m) / unit)
            net = sum(currents.values())
            assert (net - ops["boundary"][0]).count_nonzero() == 0
            e4, energy_rate, energy_flux = regional_energy(
                board, route, patch, q, odd, h
            )
            for key, mat in energy_flux.items():
                local.setdefault(f"energy_edge_{key[0]}_{key[1]}", []).append(
                    effect(w, mat) * (0.125j)
                )
            materials.setdefault("regional_energy", []).append(effect(w, e4) / 4)
            materials.setdefault("energy_inflow", []).append(
                effect(w, energy_rate) * (0.125j)
            )
            z = blank_h(board, q)
            v = h - z
            for (a, b), op in currents.items():
                local.setdefault(f"edge_{a}_{b}", []).append(effect(w, op) / 4)
            for i, u in enumerate(patch["cells"]):
                single = {k: [patch[k][i]] for k in ("cells", "images", "lifts")}
                part = map_matrix(board, single)
                body = v @ part - part @ v
                total = h @ part - part @ h
                cuts = sum(op for (inner, _), op in currents.items() if inner == u)
                assert (total + total.T - cuts - body - body.T).count_nonzero() == 0
                for key, mat in (("body", body + body.T), ("rate", total + total.T)):
                    local.setdefault(f"{key}_{u}", []).append(effect(w, mat) / 4)
            step = 1 / 2048
            times = []
            energies = []
            for sign in (-1, 1):
                ht = -0.5j * step * sign * h
                v = expm_multiply(ht, w, traceA=ht.diagonal().sum())
                times.append(effect(v, ops["shift"][0]) / ops["shift"][1])
                energies.append(effect(v, e4) / 4)
            changed = [effect(w, ops["shift"][0]) / ops["shift"][1]]
            for ticks in (1, 2):
                ht = -0.5j * step * ticks * z
                after = expm_multiply(ht, w, traceA=ht.diagonal().sum())
                changed.append(effect(after, ops["shift"][0]) / ops["shift"][1])
            materials.setdefault("blank_release_shift_rate", []).append(
                (-3 * changed[0] + 4 * changed[1] - changed[2]) / (2 * step)
            )
            materials.setdefault("independent_free_energy_difference", []).append(
                (energies[1] - energies[0]) / (2 * step)
            )
            materials.setdefault("independent_free_shift_difference", []).append(
                (times[1] - times[0]) / (2 * step)
            )
        val = {k: bridge.evaluate(v, states) for k, v in materials.items()}
        for k in states:
            assert abs(val["blank_release_shift_rate"][k] - val["boundary"][k]) < 1e-8
            np.testing.assert_allclose(
                val["domain_cut"][k] + val["image_cut"][k] + val["marked_force"][k],
                val["free_rate"][k],
                atol=1e-12,
            )
            assert (
                abs(val["free_rate"][k] - val["independent_free_shift_difference"][k])
                < 0.0003
            )
            assert (
                abs(
                    val["energy_inflow"][k]
                    - val["independent_free_energy_difference"][k]
                )
                < 0.0003
            )
        rows[str(q)] = val
        local_rows[str(q)] = {k: bridge.evaluate(v, states) for k, v in local.items()}
    cuts = {}
    for side, selection in (("domain", patch["cells"]), ("image", patch["images"])):
        selected = set(selection)
        cuts[side] = [
            [a, b]
            for a, b in board.external
            if (a // 12 in selected) != (b // 12 in selected)
        ]
    return dict(
        dwell=4,
        charge_definition="(S-S†)/(2i)",
        shift_displacement_in_cell_steps=2,
        chart=patch,
        cuts=cuts,
        finite_survey=finite_geometry(),
        signs=rows,
        local_rates_and_edges=local_rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
