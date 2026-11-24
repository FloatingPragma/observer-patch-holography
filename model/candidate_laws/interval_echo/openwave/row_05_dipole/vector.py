"""Test actual rotation of three separately measurable loop-moment components."""

import argparse
import json
from functools import lru_cache

import numpy as np
from scipy import sparse
from scipy.linalg import expm, null_space
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from .experiment import (
    CHARGE,
    HERE,
    SHOTS,
    defect_generator,
    gates_of,
    geometry,
    portable,
    sample,
)

TIMES = (0, 16, 32, 64, 128, 256)


@lru_cache(maxsize=None)
def operators(gates):
    xyz = geometry()
    z = xyz[list(gates)].sum(axis=0)
    z /= np.linalg.norm(z)
    x = next(
        v - z * (z @ v)
        for v in xyz[list(gates)]
        if np.linalg.norm(v - z * (z @ v)) > 0.0001
    )
    x /= np.linalg.norm(x)
    axes = np.stack([x, np.cross(z, x), z], axis=1)
    basis = null_space(
        np.vstack((np.eye(12) + cell().adjacency().toarray(), np.eye(12)[list(gates)]))
    )
    full = np.zeros((3, 12, 12), complex)
    areas = []
    for a, b in cell().seams:
        w = np.cross(xyz[a], xyz[b]) @ axes / 2
        areas.append(w)
        full[:, b, a] = 0.5j * CHARGE * w
        full[:, a, b] = full[:, b, a].conj()
    m = basis.T @ full @ basis
    s = np.linalg.svd(m.reshape(3, -1), compute_uv=False)
    return basis, m, np.asarray(areas), s


def vector_counts(vec, areas, rng):
    rows = []
    for a, b in cell().seams:
        p = abs(vec[a] - 1j * vec[b]) ** 2 / 2
        m = abs(vec[a] + 1j * vec[b]) ** 2 / 2
        rows.append(sample([p, m, 1 - p - m], rng))
    counts = np.array(rows)
    return dict(
        counts=rows,
        measured=(
            CHARGE * (counts[:, 0] - counts[:, 1]) @ areas / (2 * SHOTS)
        ).tolist(),
    )


def run():
    rng = np.random.default_rng(510243)
    all_layouts = sorted(
        {
            gates_of(federation(i), o)
            for i in range(3)
            for o in range(federation(i).cells)
        }
    )
    layouts = [
        dict(gates=g, singular_values=operators(g)[3].tolist()) for g in all_layouts
    ]
    board = federation(2)
    ranks = [sum(operators(gates_of(board, j))[3] > 1e-10) for j in range(board.cells)]
    owner = int(np.argmax(ranks))
    gates = gates_of(board, owner)
    basis, m, areas, _ = operators(gates)
    _, evec = np.linalg.eigh(m[0])
    raw = evec[:, -1]
    sl = slice(12 * owner, 12 * owner + 12)
    v = np.zeros(board.n, complex)
    v[sl] = basis @ raw
    h = defect_generator(board, CHARGE, candidate.KAPPA, 0).astype(complex) - (
        3 + 4 * candidate.KAPPA
    ) * sparse.eye(board.n)
    a, b = np.array(cell().seams).T + sl.start
    out = []
    for B in (-1 / 32, 0, 1 / 32):
        delta = -0.5 * np.expm1(1j * CHARGE * B * areas[:, 2])
        hB = h + sparse.csr_matrix(
            (np.r_[delta, delta.conj()], (np.r_[b, a], np.r_[a, b])), shape=h.shape
        )
        wave = expm_multiply(-1j * hB, v, start=0, stop=256, num=17)
        for t in TIMES:
            state = wave[t // 16, sl]
            pred = expm(1j * B * t * m[2]) @ raw
            vals = [np.vdot(pred, j @ pred).real for j in m]
            out.append(
                dict(
                    field=B,
                    time=t,
                    projected_forecast=vals,
                    probability_outside_compact=1
                    - np.linalg.norm(basis.T @ state) ** 2,
                    **vector_counts(state, areas, rng),
                )
            )
    rad = np.sqrt(
        np.sum(areas**2, axis=0) * np.log(2 * 3 * len(out) / 0.01) / (2 * SHOTS)
    )
    return dict(
        owner=owner,
        gates=gates,
        layouts=layouts,
        rows=out,
        shots=SHOTS,
        seed=510243,
        simultaneous_radius=rad.tolist(),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    data = portable(run())
    path = HERE / "vector.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(data, indent=1) + "\n")
    elif json.loads(path.read_text()) != data:
        raise ValueError("Precession record differs")
    print("Verified the compact moment ranks and three-axis current precession.")


if __name__ == "__main__":
    main()
